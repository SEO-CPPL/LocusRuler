#!/usr/bin/env python3
"""Intra-genome duplication verification and paralog pruning."""

import itertools
from pathlib import Path
from typing import Optional

# Status constants (mirrored from ruler.py to avoid circular import)
INTACT      = "INTACT"
FRAGMENTED  = "FRAGMENTED"
DUPLICATED  = "DUPLICATED"
PARTIAL_DEL = "PARTIAL_DEL"
ABSENT      = "ABSENT"


def _fna_path_for(acc: str, d: str) -> Optional[Path]:
    """Locate the FNA file for an accession in directory ``d``.
    Tries the common NCBI suffixes (`_genomic.fna.gz`, `.fna`, etc.)
    and falls back to a prefix scan.  Returns None if nothing matches.
    """
    d_path = Path(d)
    if not d_path.exists():
        return None
    for s in (f"{acc}_genomic.fna.gz", f"{acc}_genomic.fna",
              f"{acc}.fna.gz", f"{acc}.fna"):
        if (d_path / s).exists():
            return d_path / s
    for p in d_path.iterdir():
        if p.name.startswith(acc) and ".fna" in p.name:
            return p
    return None


def _load_ruler_genome_meta(db_path: str, accessions: list) -> dict:
    """Fetch the identity columns for ``accessions`` from the target SQLite DB.

    These are the columns every grid opens with, so this has to return the
    same set ``content._load_genome_meta`` does. Duplication verification
    itself no longer needs genus info: the threshold is a single value,
    not split intra/inter.
    """
    import sqlite3
    if not db_path or not Path(db_path).exists():
        return {}
    meta: dict = {}
    try:
        con = sqlite3.connect(db_path)
        con.row_factory = sqlite3.Row
        placeholders = ",".join("?" * len(accessions))
        rows = con.execute(
            f"SELECT accession, species, strain, assembly_level FROM genomes "
            f"WHERE accession IN ({placeholders})",
            accessions,
        ).fetchall()
        con.close()
        for r in rows:
            meta[r["accession"]] = {
                "species": r["species"] or "",
                "strain":  r["strain"]  or "",
                "assembly_level": r["assembly_level"] or "",
            }
    except Exception:
        pass
    return meta


def verify_duplications(results: dict, locus_cfg: dict, settings: dict,
                        target_name: str = "") -> None:
    """Confirm or refute multi-piece duplications and prune paralog pieces."""
    # Deferred import avoids circular ruler.py ↔ blast.py dependency
    from blast import piece_pair_identity

    tgt_name = target_name or locus_cfg.get("reference", {}).get("target", "")
    tgt_cfg  = next((t for t in settings.get("targets", [])
                     if t["name"] == tgt_name), None)
    if not tgt_cfg or not tgt_cfg.get("fna_dir"):
        return

    cbcfg = {**settings.get("cluster_blast", {}),
             **locus_cfg.get("cluster_blast", {})}
    # One threshold: piece-vs-piece is intragenomic, and 90% marks recent duplication.
    dup_min_id = float(cbcfg.get("duplication_min_identity", 90))

    ref_bp      = int(locus_cfg.get("_auto", {}).get("cluster_ref_bp") or 1)
    blastn_exe  = settings.get("tools", {}).get("blastn", "blastn")

    for acc, res in results.items():
        pieces = res.get("_pieces") or []
        if len(pieces) < 2:
            continue

        # ── Build query-space overlap graph
        overlapping_with: dict[int, dict[int, float]] = {}
        exception_tags = {
            a.get("locus_tag")
            for a in locus_cfg.get("_auto", {}).get("anchors", [])
            if a.get("exception")
        }

        for i, j in itertools.combinations(range(len(pieces)), 2):
            pa, pb = pieces[i], pieces[j]
            genes_a = set(pa.get("_genes") or []) - exception_tags
            genes_b = set(pb.get("_genes") or []) - exception_tags
            if not (genes_a & genes_b):
                continue
            qa_lo = min(int(pa["qstart"]), int(pa["qend"]))
            qa_hi = max(int(pa["qstart"]), int(pa["qend"]))
            qb_lo = min(int(pb["qstart"]), int(pb["qend"]))
            qb_hi = max(int(pb["qstart"]), int(pb["qend"]))
            if min(qa_hi, qb_hi) - max(qa_lo, qb_lo) > 50:
                fna = _fna_path_for(acc, tgt_cfg["fna_dir"])
                if not fna:
                    continue
                pid = piece_pair_identity(fna, pa, pb, blastn=blastn_exe)
                overlapping_with.setdefault(i, {})[j] = pid
                overlapping_with.setdefault(j, {})[i] = pid

        if not overlapping_with:
            # All pieces cover different reference regions → genuine fragments
            continue

        max_id = max(pid
                     for d in overlapping_with.values()
                     for pid in d.values())

        if max_id >= dup_min_id:
            # Sequence identity confirms a real intra-genome duplication
            if res["status"] != DUPLICATED:
                res["status"] = DUPLICATED
            continue

        # ── Low identity: drop paralog pieces
        visited: set[int] = set()
        to_drop: set[int] = set()

        for start in sorted(overlapping_with.keys()):
            if start in visited:
                continue
            component: set[int] = set()
            stack = [start]
            while stack:
                cur = stack.pop()
                if cur in visited:
                    continue
                visited.add(cur)
                component.add(cur)
                for neighbor in overlapping_with.get(cur, {}):
                    if neighbor not in visited:
                        stack.append(neighbor)
            best = max(
                component,
                key=lambda k: (len(pieces[k].get("_genes") or []), -k),
            )
            to_drop.update(component - {best})

        if to_drop:
            keep_idx = sorted(i for i in range(len(pieces)) if i not in to_drop)
            breaks_all = res.get("_internal_breaks") or []
            res["_pieces"] = [pieces[i] for i in keep_idx]
            if len(breaks_all) == len(pieces):
                res["_internal_breaks"] = [breaks_all[i] for i in keep_idx]
            res["piece_count"]  = len(res["_pieces"])
            res["contig_count"] = len({p["sseqid"] for p in res["_pieces"]})
            # Recompute coverage from kept pieces' query-union lengths.
            new_cov = (sum(int(p.get("length", 0)) for p in res["_pieces"])
                       / max(1, ref_bp))
            res["coverage"] = min(new_cov, 1.0)
            print(f"[duplication] {acc}: dropped {len(to_drop)} paralog piece(s) "
                  f"(q-overlap, max identity={max_id:.1f}% < {dup_min_id:.0f}%); "
                  f"{len(res['_pieces'])} piece(s) remain, "
                  f"cov={res['coverage']:.0%}")

        # Reclassify status for the (possibly pruned) remaining pieces
        n_remain = len(res.get("_pieces") or pieces)
        cov      = res.get("coverage") or 0.0
        if n_remain == 1:
            if cov >= 0.85:
                res["status"] = INTACT
            elif cov >= 0.30:
                res["status"] = PARTIAL_DEL
            else:
                res["status"] = ABSENT
        else:
            # Multiple non-overlapping genuine fragments
            res["status"] = FRAGMENTED
