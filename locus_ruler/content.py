#!/usr/bin/env python3
"""GFF content extraction and per-gene classification."""

import re
import sqlite3
from pathlib import Path
from typing import Optional

from gff import _find_gff, parse_gff_region
from classify import load_rules, classify_genes
from fragmentation import classify_fragmentation
from flank import (
    get_contig_lengths,
    gff_gene_at,
    adjacent_pseudo_partner,
    adjacent_split_cds_pair,
)
from gene_state import (
    inner_anchor_tblastn_metrics,
    cluster_dna_metrics,
    classify_gene_state,
)
from locus_status import (
    apply_locus_coverage_floor,
    classify_locus_status,
    normalize_status_role,
    summarize_status_detail,
)
from domain_recovery import (
    domain_architecture,
    parse_pfam_domtblout,
    pfam_tokens,
    pfam_accessions_from_anchors,
    resolve_domtblout_path,
    run_hmmsearch_if_configured,
    same_strand_neighbors,
    split_domain_class,
    used_gathering_thresholds,
    write_domain_recovery_diagnostics,
)
from flank_blast import (
    _make_faa_index,
    _build_ref_prot_db,
    _reverse_flank_blast,
    _reverse_cluster_blast,
    _select_aux_genome,
    _build_aux_prot_db,
    _build_anchor_prot_db,
    _product_to_flank_label,
)
from cohort_rescue import (
    build_spans_by_group as build_cohort_spans_by_genus,
    cohort_coverage,
    describe as describe_cohort_spans,
)
from locus_scale import (
    DERIVE,
    describe as describe_window,
    reference_length,
    resolve as resolve_window,
)
from writers import (
    diagnostics_dir,
    tables_dir,
    write_loci_csv,
    write_gene_diagnostics_csv,
    write_hsp_diagnostics_csv,
    write_loci_xlsx,
    write_output_guide,
    write_pieces_csv,
    write_clade_markers_tsv,
    write_marker_matrix_csv,
)


# Statuses from which content extraction makes sense
_CONTENT_STATUSES = {
    "INTACT", "PSEUDOGENIZED", "PARTIAL_DEL", "FRAGMENTED",
    "LARGE_DEL", "DUPLICATED", "ANNOTATION_GAP",
}


def _make_separator_cell(label: str, contig: str, pos: int) -> dict:
    """Pseudo-gene cell used to delimit pieces of a fragmented cluster."""
    return {
        "_zone":        "separator",
        "locus_tag":    "",
        "product":      label,
        "family":       None,
        "is_pseudo":    False,
        "start":        pos,
        "end":          pos,
        "contig":       contig,
        "strand":       "",
        "feature_type": "separator",
        "attrs":        "",
    }


def _rebuild_layout(pl: dict) -> None:
    """Sort in_cluster to match left_n/right_n's orientation; layout is not reversed again."""
    pl["in_cluster"].sort(key=lambda g: int(g["start"]), reverse=(pl["orientation"] < 0))
    pl["layout"] = pl["left_n"] + pl["in_cluster"] + pl["right_n"]


def _load_genome_meta(db_path: str, accessions: list[str]) -> dict[str, dict]:
    """Load species / strain / assembly_level for each accession from the SQLite DB."""
    if not db_path or not Path(db_path).exists():
        return {}
    meta: dict[str, dict] = {}
    try:
        con = sqlite3.connect(db_path)
        con.row_factory = sqlite3.Row
        placeholders = ",".join("?" * len(accessions))
        rows = con.execute(
            f"SELECT accession, species, strain, assembly_level "
            f"FROM genomes WHERE accession IN ({placeholders})",
            accessions,
        ).fetchall()
        con.close()
        for r in rows:
            meta[r["accession"]] = {
                "species":        r["species"] or "",
                "strain":         r["strain"]  or "",
                "assembly_level": r["assembly_level"] or "",
            }
    except Exception as e:
        print(f"[content] WARN: could not load genome metadata: {e}")
    return meta


# ── Main orchestration ─────────────────────────────────────────────────────

def run_content(
    ruler_results: dict[str, dict],
    locus_cfg: dict,
    settings: dict,
    blast_hits: dict[str, dict[str, list[dict]]],
    work_dir: Path,
    output_dir: Path,
    target_name: str = "",
    config_dir: Path = None,
    cluster_hits: Optional[dict[str, list[dict]]] = None,
) -> dict[str, dict]:
    """Extract gene content from GFF files for every genome with a confirmed locus."""
    auto     = locus_cfg.get("_auto", {})
    anchors  = auto.get("anchors", [])
    exp_bp   = auto.get("expected_bp", {})
    locus_id = locus_cfg["locus_id"]

    # target_name identifies the dataset actually scanned; reference.target is a fallback.
    ref_target = target_name or locus_cfg.get("reference", {}).get("target", "")
    tgt_cfg    = next(
        (t for t in settings.get("targets", []) if t["name"] == ref_target),
        None,
    )
    if tgt_cfg is None:
        if settings.get("targets"):
            tgt_cfg = settings["targets"][0]
            print(f"[content] WARN: target '{ref_target}' not found; "
                  f"using '{tgt_cfg['name']}'")
        else:
            print("[content] ERROR: no [[targets]] in settings.toml; skipping GFF step")
            return ruler_results

    gff_dir    = tgt_cfg.get("gff_dir", "")
    rules_file = settings.get("classification", {}).get("rules_file", "")
    if rules_file and not Path(rules_file).exists():
        print(f"[content] WARN: global classification rules declared in "
              f"settings.toml but not found at '{rules_file}' - falling "
              "back to per-locus rules / tblastn override")
        rules_file = ""
    rules = load_rules(rules_file, locus_id=locus_id, config_dir=config_dir)

    # ── Threshold resolution
    ruler_cfg     = {**settings.get("ruler", {}), **locus_cfg.get("ruler", {})}
    anchor_min_id = float(ruler_cfg.get("anchor_min_identity", 40))
    min_decayed_coverage = float(ruler_cfg.get("min_decayed_coverage", 0.30))
    cbcfg         = {**settings.get("cluster_blast", {}),
                     **locus_cfg.get("cluster_blast", {})}
    min_gene_cov          = float(cbcfg.get("min_gene_coverage",          0.30))
    tblastn_intact_cov    = float(cbcfg.get("tblastn_intact_coverage",    0.70))
    tblastn_min_id        = float(cbcfg.get("tblastn_min_identity",       30))
    tblastn_intact_id_intra = float(cbcfg.get(
        "tblastn_intact_identity_intra",
        cbcfg.get("tblastn_intact_identity",
                  cbcfg.get("tblastn_cross_check_identity", 70))))
    tblastn_intact_id_inter = float(cbcfg.get(
        "tblastn_intact_identity_inter",
        cbcfg.get("tblastn_intact_identity",
                  cbcfg.get("tblastn_cross_check_identity", 50))))
    ref_genus         = (auto.get("reference_genus") or "").strip()
    ref_cluster_bp    = reference_length(locus_cfg)
    frag_cfg          = {**settings.get("fragmentation", {}),
                         **locus_cfg.get("fragmentation", {})}
    content_cfg       = settings.get("content", {})
    n_flank           = int(content_cfg.get("flank_genes_shown", 2))
    # Both windows scale with the reference locus unless settings pin them.
    _ref_gap_cfg = content_cfg.get("adjacent_tblastn_rescue_ref_gap_bp", DERIVE)
    _subj_gap_cfg = content_cfg.get("adjacent_tblastn_rescue_subject_gap_bp", DERIVE)
    adjacent_ref_gap = resolve_window(
        _ref_gap_cfg, ref_cluster_bp,
        float(content_cfg.get("adjacent_tblastn_rescue_ref_gap_fraction", 0.20)),
    )
    adjacent_subj_gap = resolve_window(
        _subj_gap_cfg, ref_cluster_bp,
        float(content_cfg.get("adjacent_tblastn_rescue_subject_gap_fraction", 0.50)),
    )
    for _name, _cfg, _val in (
        ("adjacent_tblastn_rescue_ref_gap", _ref_gap_cfg, adjacent_ref_gap),
        ("adjacent_tblastn_rescue_subject_gap", _subj_gap_cfg, adjacent_subj_gap),
    ):
        print(f"[content] window - {describe_window(_name, _cfg, _val, ref_cluster_bp)}")
    flank_label_min_id = float(content_cfg.get("flank_label_min_identity", 0.0))
    aux_anchors       = [a for a in anchors if a.get("role") == "aux"]
    shared_anchors    = [
        a for a in anchors
        if normalize_status_role(
            a.get("status_role", ""),
            a.get("role", ""),
            bool(a.get("exception")),
        ) == "SHARED"
    ]
    n_expected        = len([a for a in anchors if a.get("role") != "aux"])
    exp_flank         = exp_bp.get("flank")
    domain_cfg        = {**settings.get("domain_recovery", {}),
                         **locus_cfg.get("domain_recovery", {})}
    domain_recovery_enabled = bool(domain_cfg.get("enabled", True))
    domain_domtblout = (
        resolve_domtblout_path(settings, target_name, domain_cfg)
        if domain_recovery_enabled else None
    )
    if domain_recovery_enabled and domain_domtblout is None:
        domain_domtblout = run_hmmsearch_if_configured(
            settings,
            tgt_cfg.get("faa", ""),
            work_dir / target_name,
            locus_id,
            cpu=int(settings.get("run", {}).get("cpu", 1)),
            domain_cfg=domain_cfg,
            pfam_accessions=pfam_accessions_from_anchors(anchors),
        )
    if not domain_domtblout:
        domain_hits_by_lt = {}
    elif used_gathering_thresholds(domain_domtblout):
        # Pfam already decided this family by family; re-filtering only undoes it.
        domain_hits_by_lt = parse_pfam_domtblout(domain_domtblout)
    else:
        domain_hits_by_lt = parse_pfam_domtblout(
            domain_domtblout,
            max_seq_evalue=float(domain_cfg.get("evalue", "1e-5")),
            max_domain_cevalue=(
                float(domain_cfg["domain_cevalue"])
                if domain_cfg.get("domain_cevalue") not in (None, "") else None
            ),
            max_domain_ievalue=(
                float(domain_cfg["domain_ievalue"])
                if domain_cfg.get("domain_ievalue") not in (None, "") else None
            ),
        )
    domain_recovery_rows: list[dict] = []
    if domain_hits_by_lt:
        print(f"[content] domain recovery: loaded Pfam domains for "
              f"{len(domain_hits_by_lt):,} proteins from {domain_domtblout}")

    fl_tag = next((a["locus_tag"] for a in anchors if a["role"] == "flank_L"),  None)
    fr_tag = next((a["locus_tag"] for a in anchors if a["role"] == "flank_R"),  None)
    cl_ref_tag = next((a["locus_tag"] for a in anchors if a["role"] == "cluster_L"), None)
    cr_ref_tag = next((a["locus_tag"] for a in anchors if a["role"] == "cluster_R"), None)
    cl_anchor_ref  = next((a for a in anchors if a["role"] == "cluster_L"), None)
    cluster_origin = int(cl_anchor_ref["start"]) if cl_anchor_ref else 0

    genome_meta = _load_genome_meta(
        tgt_cfg.get("db", ""),
        [acc for acc in ruler_results if ruler_results[acc]["status"] != "UNKNOWN"],
    )
    domain_lts_by_acc: dict[str, list[str]] = {}
    domain_pos_by_lt: dict[str, dict] = {}
    if domain_hits_by_lt and tgt_cfg.get("db") and Path(tgt_cfg.get("db")).exists():
        try:
            con = sqlite3.connect(tgt_cfg["db"])
            lts = list(domain_hits_by_lt)
            for i in range(0, len(lts), 500):
                batch = lts[i:i + 500]
                placeholders = ",".join("?" * len(batch))
                rows = con.execute(
                    f"SELECT locus_tag, genome_acc, contig, start, end, strand, product FROM proteins "
                    f"WHERE locus_tag IN ({placeholders})",
                    batch,
                ).fetchall()
                for lt, gacc, ctg, start, end, strand, product in rows:
                    domain_lts_by_acc.setdefault(gacc, []).append(lt)
                    domain_pos_by_lt[lt] = {
                        "locus_tag": lt,
                        "genome_acc": gacc,
                        "contig": ctg,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                        "product": product or "",
                    }
            con.close()
        except Exception as e:
            print(f"[content] WARN: domain recovery locus_tag mapping failed: {e}")

    # ── Reverse flank blast setup
    tools_cfg      = settings.get("tools", {})
    makeblastdb_exe = tools_cfg.get("makeblastdb", "makeblastdb")
    blastp_exe     = tools_cfg.get("blastp", "blastp")
    tgt_faa_path   = tgt_cfg.get("faa", "")

    ref_cfg_name   = locus_cfg.get("reference", {}).get("target", "")
    ref_accession  = locus_cfg.get("reference", {}).get("accession", "")
    ref_tgt_cfg_wb = next(
        (t for t in settings.get("targets", []) if t["name"] == ref_cfg_name),
        None,
    )

    ref_prot_db: Optional[str] = None
    tgt_faa_index: dict[str, int] = {}

    if ref_tgt_cfg_wb and tgt_faa_path and Path(tgt_faa_path).exists():
        ref_prot_db = _build_ref_prot_db(
            ref_tgt_cfg_wb.get("faa", ""),
            ref_accession,
            work_dir / locus_id,
            makeblastdb_exe,
        )
        if ref_prot_db:
            print(f"[content] Indexing target FAA for reverse flank blast -> {tgt_faa_path}")
            tgt_faa_index = _make_faa_index(tgt_faa_path)
            print(f"[content]   indexed {len(tgt_faa_index):,} target proteins")

    # ── Auxiliary proteome for inter-genus reverse flank blast
    aux_prot_db: Optional[str] = None
    if tgt_faa_path and Path(tgt_faa_path).exists():
        _aux_acc = _select_aux_genome(ruler_results)
        if _aux_acc:
            _tgt_db_path = tgt_cfg.get("db", "")
            print(f"[content] Building aux prot DB from {_aux_acc} -> {work_dir / locus_id}")
            aux_prot_db = _build_aux_prot_db(
                _aux_acc,
                _tgt_db_path,
                tgt_faa_path,
                work_dir / locus_id,
                makeblastdb_exe,
            )
            if aux_prot_db:
                print(f"[content]   aux_prot_db ready for inter-genus flank blast")
            if not tgt_faa_index:
                print(f"[content] Indexing target FAA for aux flank blast -> {tgt_faa_path}")
                tgt_faa_index = _make_faa_index(tgt_faa_path)
                print(f"[content]   indexed {len(tgt_faa_index):,} target proteins")

    # ── Anchor proteome for reverse cluster blast
    anchor_prot_db: Optional[str] = None
    anchor_lt_to_family: dict[str, str] = {}
    _anchors_faa = auto.get("_anchors_faa", "")
    if _anchors_faa:
        anchor_prot_db = _build_anchor_prot_db(
            _anchors_faa,
            work_dir / locus_id,
            makeblastdb_exe,
        )
        for _a in anchors:
            if _a.get("role") not in ("cluster_L", "cluster_R", "inner"):
                continue
            _lt  = _a.get("locus_tag", "")
            _fam = _a.get("family") or _lt
            if _lt and _fam:
                anchor_lt_to_family[_lt] = _fam
        if anchor_prot_db:
            print(f"[content] anchor_prot_db ready for reverse cluster blast "
                  f"({len(anchor_lt_to_family)} cluster anchor families)")
        if anchor_prot_db and not tgt_faa_index and tgt_faa_path and Path(tgt_faa_path).exists():
            print(f"[content] Indexing target FAA for reverse cluster blast -> {tgt_faa_path}")
            tgt_faa_index = _make_faa_index(tgt_faa_path)
            print(f"[content]   indexed {len(tgt_faa_index):,} target proteins")

    # Per-accession contig length cache (populated on demand)
    blastdb_root   = work_dir / "blastdbs"
    blastdbcmd_exe = tools_cfg.get("blastdbcmd", "blastdbcmd")
    _contig_lens_cache: dict[str, dict[str, int]] = {}

    aux_lenient_id = float(cbcfg.get("single_gene_min_identity_lenient", 40))
    aux_lenient_cov = float(cbcfg.get("single_gene_min_coverage_lenient", 0.40))
    aux_qlen: dict[str, int] = {}
    for a in aux_anchors:
        try:
            aux_qlen[a["locus_tag"]] = max(
                1, (int(a["end"]) - int(a["start"]) + 1) // 3
            )
        except (KeyError, ValueError, TypeError):
            aux_qlen[a.get("locus_tag", "")] = 1

    # Consensus query window per aux anchor, built once per genus. One window
    # for a cohort spanning several genera would be set by whichever genus
    # brought the most genomes, and the rest measured against an alignment their
    # orthologue never had.
    aux_genus_of = {
        accession: (str(meta.get("species") or "").strip().split() or [""])[0]
        for accession, meta in genome_meta.items()
    }
    aux_cohort_spans = build_cohort_spans_by_genus(
        blast_hits, aux_qlen, aux_genus_of,
        min_coverage=aux_lenient_cov,
        min_identity=aux_lenient_id,
    ) if aux_anchors else {}
    for genus, spans in sorted(aux_cohort_spans.items()):
        for line in describe_cohort_spans(spans, aux_qlen):
            print(f"[content] cohort span - {genus or 'ungrouped'} - {line}")

    def _contig_lens_for(acc: str) -> dict[str, int]:
        if acc not in _contig_lens_cache:
            db_stem = str(blastdb_root / acc / acc)
            _contig_lens_cache[acc] = get_contig_lengths(db_stem, blastdbcmd_exe)
        return _contig_lens_cache[acc]

    # Maps an assigned family back to its reference role label, once per locus.
    _anchor_fam_to_ctx_label: dict[str, str] = {}
    for _a in anchors:
        _fam  = _a.get("family") or ""
        _lt   = _a.get("locus_tag") or ""
        _role = _a["role"]
        _label: Optional[str] = None
        if   _role == "flank_L":  _label = "L1"
        elif _role == "flank_L2": _label = "L2"
        elif _role == "flank_R":  _label = "R1"
        elif _role == "flank_R2": _label = "R2"
        elif _role in ("cluster_L", "cluster_R", "inner") and not _a.get("exception"):
            _label = "INNER"
        if _label:
            # Register both the family string and the locus_tag; reverse-blastp returns the latter.
            if _fam:
                _anchor_fam_to_ctx_label[_fam] = _label
            if _lt and _lt != _fam:
                _anchor_fam_to_ctx_label[_lt] = _label

    # ── Per-genome content extraction loop
    all_locus_genes: dict[str, list[dict]] = {}

    for acc, res in ruler_results.items():
        piece_status = res.get("status", "UNKNOWN")
        if piece_status not in _CONTENT_STATUSES and not res.get("_pieces"):
            all_locus_genes[acc] = []
            if piece_status != "UNKNOWN":
                res["status"], res["_status_detail"] = classify_locus_status({}, anchors)
            continue

        home_contig = res.get("contig") or ""
        if "|" in home_contig:
            home_contig = home_contig.split("|")[0]
        if not home_contig:
            all_locus_genes[acc] = []
            res["status"], res["_status_detail"] = classify_locus_status({}, anchors)
            continue

        gff_path = _find_gff(gff_dir, acc)
        if gff_path is None:
            res["notes"] = (res.get("notes") or "") + "; GFF not found"
            all_locus_genes[acc] = []
            res["status"] = "UNKNOWN"
            continue

        pieces_raw = list(res.get("_pieces") or [])
        if not pieces_raw:
            all_locus_genes[acc] = []
            res["status"], res["_status_detail"] = classify_locus_status({}, anchors)
            continue

        def _ref_q(p):
            try:
                return min(int(p["qstart"]), int(p["qend"]))
            except Exception:
                return 0
        pieces_sorted = sorted(pieces_raw, key=_ref_q)

        # -- Fragmentation annotation (FRAGMENTED genomes only) --------
        if piece_status == "FRAGMENTED":
            _frag = classify_fragmentation(
                pieces        = pieces_sorted,
                assembly_level= genome_meta.get(acc, {}).get("assembly_level", ""),
                coverage      = float(res.get("coverage") or 0.0),
                settings      = frag_cfg,
            )
            res["_fragmentation_type"]  = _frag["fragmentation_type"]
            res["_contig_edge_support"] = _frag["contig_edge_support"]
            res["_q_pair_max_overlap"]  = _frag["q_pair_max_overlap"]
            for _p, _ep, _ed in zip(pieces_sorted,
                                    _frag["piece_edge_proximal"],
                                    _frag["piece_edge_distance"]):
                _p["_edge_proximal"] = _ep
                _p["_edge_distance"] = _ed

        # Parse GFF per unique contig, deduplicate locus_tags
        unique_contigs = {p["sseqid"] for p in pieces_sorted}
        gff_by_contig: dict[str, list[dict]] = {}
        for c in unique_contigs:
            raw = parse_gff_region(gff_path, c, 1, 10**9)
            if rules:
                classify_genes(raw, rules)
            seen_lt: dict[str, dict] = {}
            for g in raw:
                lt = g["locus_tag"]
                if not lt:
                    continue
                if lt not in seen_lt:
                    seen_lt[lt] = g
                else:
                    existing = seen_lt[lt]
                    if g["feature_type"] == "CDS" and existing["feature_type"] != "CDS":
                        seen_lt[lt] = g
                    elif (g["is_pseudo"] and not existing["is_pseudo"]
                          and existing["feature_type"] == "gene"):
                        seen_lt[lt] = g
            gff_by_contig[c] = sorted(seen_lt.values(), key=lambda g: g["start"])

        def _load_contig_genes(contig_id: str) -> list[dict]:
            if contig_id not in gff_by_contig:
                raw = parse_gff_region(gff_path, contig_id, 1, 10**9)
                if rules:
                    classify_genes(raw, rules)
                seen_lt: dict[str, dict] = {}
                for gg in raw:
                    lt = gg["locus_tag"]
                    if not lt:
                        continue
                    if lt not in seen_lt:
                        seen_lt[lt] = gg
                    else:
                        existing = seen_lt[lt]
                        if gg["feature_type"] == "CDS" and existing["feature_type"] != "CDS":
                            seen_lt[lt] = gg
                        elif (gg["is_pseudo"] and not existing["is_pseudo"]
                              and existing["feature_type"] == "gene"):
                            seen_lt[lt] = gg
                gff_by_contig[contig_id] = sorted(
                    seen_lt.values(), key=lambda gg: gg["start"])
            return gff_by_contig.get(contig_id, [])

        acc_hits = blast_hits.get(acc, {})

        def _span_gap(a_lo: int, a_hi: int, b_lo: int, b_hi: int) -> int:
            if a_hi < b_lo:
                return b_lo - a_hi
            if b_hi < a_lo:
                return a_lo - b_hi
            return 0

        def _shared_anchor_hit_state(anchor: dict) -> dict:
            """Genome-wide functional state for SHARED precursor anchors."""
            lt = anchor.get("locus_tag", "")
            hits = acc_hits.get(lt, []) or []
            if not hits:
                return {"state": "ABSENT"}

            state_rank = {"ABSENT": 0, "PSEUDOGENE": 1, "DIVERGENT": 2, "INTACT": 3}
            best = None
            for h in hits:
                try:
                    pid = float(h.get("pident", 0))
                    qstart, qend = int(h["qstart"]), int(h["qend"])
                    qlen = int(h.get("qlen") or 0)
                except (KeyError, ValueError, TypeError):
                    continue
                if qlen <= 0:
                    try:
                        qlen = max(1, (int(anchor["end"]) - int(anchor["start"]) + 1) // 3)
                    except (KeyError, ValueError, TypeError):
                        qlen = 1
                cov = abs(qend - qstart + 1) / max(1, qlen)
                try:
                    contig_id = h.get("sseqid", "")
                    s_lo = min(int(h["sstart"]), int(h["send"]))
                    s_hi = max(int(h["sstart"]), int(h["send"]))
                except (KeyError, ValueError, TypeError):
                    contig_id, s_lo, s_hi = "", 0, 0

                gff_pseudo = False
                if contig_id and s_lo and s_hi:
                    g = gff_gene_at(_load_contig_genes(contig_id), s_lo, s_hi)
                    gff_pseudo = bool(g and g.get("is_pseudo"))

                if anchor.get("lenient"):
                    intact_id = aux_lenient_id
                    intact_cov = aux_lenient_cov
                else:
                    intact_id = float(cbcfg.get(
                        "single_gene_min_identity_intra" if same_genus
                        else "single_gene_min_identity_inter",
                        90 if same_genus else 75,
                    ))
                    intact_cov = float(cbcfg.get("single_gene_min_coverage", 0.70))

                state = classify_gene_state(
                    tblastn_cov=min(1.0, cov),
                    tblastn_pid=pid,
                    cluster_dna_present=False,
                    gff_pseudo=gff_pseudo,
                    intact_coverage=intact_cov,
                    intact_identity=intact_id,
                    min_coverage=min_gene_cov,
                    min_identity=tblastn_min_id,
                    lenient=bool(anchor.get("lenient", False)),
                )
                candidate = {
                    "state": state,
                    "pid": pid,
                    "cov": min(1.0, cov),
                    "source_locus": lt,
                    "hit_contig": contig_id,
                    "hit_lo": s_lo or "",
                    "hit_hi": s_hi or "",
                    "bitscore": h.get("bitscore", ""),
                }
                if best is None:
                    best = candidate
                    continue
                cand_key = (
                    state_rank.get(candidate["state"], 0),
                    float(candidate.get("bitscore") or 0),
                )
                best_key = (
                    state_rank.get(best["state"], 0),
                    float(best.get("bitscore") or 0),
                )
                if cand_key > best_key:
                    best = candidate

            return best or {"state": "ABSENT"}

        def _best_anchor_hit_on(
            tag: Optional[str],
            c: str,
            min_cov: Optional[float] = None,
            min_pid: Optional[float] = None,
        ) -> Optional[dict]:
            if not tag:
                return None
            rows = acc_hits.get(tag, [])
            cands = []
            for r in rows:
                if r.get("sseqid") != c:
                    continue
                try:
                    pid = float(r.get("pident", 0))
                    length = int(r.get("length", 0))
                    qlen = int(r.get("qlen", 0))
                except (TypeError, ValueError):
                    continue
                if min_pid is not None and pid < min_pid:
                    continue
                if min_cov is not None and length / max(1, qlen) < min_cov:
                    continue
                cands.append(r)
            if not cands:
                return None
            best = max(cands, key=lambda r: float(r.get("bitscore", 0)))
            lo = min(int(best["sstart"]), int(best["send"]))
            hi = max(int(best["sstart"]), int(best["send"]))
            return {
                "lo": lo,
                "hi": hi,
                "contig": c,
                "pident": float(best.get("pident", 0)),
                "length": int(best.get("length", 0)),
                "qlen": int(best.get("qlen", 0)),
            }

        def _anchor_q_interval(a: dict) -> tuple[int, int]:
            q_lo = int(a["start"]) - cluster_origin
            q_hi = int(a["end"]) - cluster_origin
            return min(q_lo, q_hi), max(q_lo, q_hi)

        target_species = genome_meta.get(acc, {}).get("species", "") or ""
        target_genus   = target_species.strip().split()[0] if target_species else ""
        same_genus     = bool(ref_genus and target_genus and ref_genus == target_genus)

        # ── tblastn 1:1 family override (INNER ANCHORS ONLY)
        _ref_family_by_lt: dict[str, str] = {}
        for a in anchors:
            if a.get("role") not in ("cluster_L", "cluster_R", "inner"):
                continue
            lt = a.get("locus_tag")
            if not lt:
                continue
            _ref_family_by_lt[lt] = a["family"] if a.get("family") else lt

        # Piece subject intervals indexed by contig, the position filter for the forward override.
        _piece_iv_by_contig: dict[str, list[tuple[int, int]]] = {}
        for _p in pieces_sorted:
            _c = _p.get("sseqid")
            if not _c:
                continue
            ivs = _p.get("subject_intervals")
            if ivs:
                for lo, hi in ivs:
                    _piece_iv_by_contig.setdefault(_c, []).append((int(lo), int(hi)))
            else:
                _plo = min(int(_p["sstart"]), int(_p["send"]))
                _phi = max(int(_p["sstart"]), int(_p["send"]))
                _piece_iv_by_contig.setdefault(_c, []).append((_plo, _phi))

        def _gene_inside_any_piece(contig_id: str,
                                   g_lo: int, g_hi: int) -> bool:
            for plo, phi in _piece_iv_by_contig.get(contig_id, []):
                if g_lo <= phi and g_hi >= plo:
                    return True
            return False

        if _ref_family_by_lt:
            _hits_by_anchor_contig: dict[tuple[str, str], list[tuple[int, int]]] = {}
            for ref_lt in _ref_family_by_lt:
                for h in acc_hits.get(ref_lt, []):
                    cid = h.get("sseqid")
                    if not cid:
                        continue
                    try:
                        lo = min(int(h["sstart"]), int(h["send"]))
                        hi = max(int(h["sstart"]), int(h["send"]))
                    except (KeyError, TypeError, ValueError):
                        continue
                    _hits_by_anchor_contig.setdefault((ref_lt, cid), []).append((lo, hi))

            def _union_overlap(hits: list[tuple[int, int]],
                               g_lo: int, g_hi: int) -> int:
                clipped = [(max(g_lo, lo), min(g_hi, hi))
                           for lo, hi in hits if hi >= g_lo and lo <= g_hi]
                if not clipped:
                    return 0
                clipped.sort()
                merged = [list(clipped[0])]
                for s, e in clipped[1:]:
                    if s <= merged[-1][1] + 1:
                        merged[-1][1] = max(merged[-1][1], e)
                    else:
                        merged.append([s, e])
                return sum(e - s + 1 for s, e in merged)

            for contig_id, gene_list in gff_by_contig.items():
                for g in gene_list:
                    g_lo, g_hi = int(g["start"]), int(g["end"])
                    # Position filter: only piece-internal genes get the forward family override.
                    if not _gene_inside_any_piece(contig_id, g_lo, g_hi):
                        continue
                    g_len = max(1, g_hi - g_lo + 1)
                    best_cov, best_fam = -1.0, None
                    for (ref_lt, c), hits in _hits_by_anchor_contig.items():
                        if c != contig_id:
                            continue
                        cov = _union_overlap(hits, g_lo, g_hi) / g_len
                        if cov >= min_gene_cov and cov > best_cov:
                            best_cov = cov
                            best_fam = _ref_family_by_lt[ref_lt]
                    if best_fam:
                        g["family"] = best_fam

        # ── Per-piece layout 
        piece_layouts: list[dict] = []

        for piece_idx, piece in enumerate(pieces_sorted):
            contig_p = piece["sseqid"]
            s_dir = 1 if int(piece["send"]) > int(piece["sstart"]) else -1

            # Wrap-origin piece: either interval is in-cluster, the contig middle is flank.
            subj_ivs = piece.get("subject_intervals")
            wraps    = bool(piece.get("_wraps_origin"))
            clen     = int(piece.get("_contig_len") or 0)

            if wraps and subj_ivs and len(subj_ivs) >= 2:
                end_lo, end_hi     = int(subj_ivs[0][0]), int(subj_ivs[0][1])
                start_lo, start_hi = int(subj_ivs[1][0]), int(subj_ivs[1][1])
                # p_lo / p_hi span the end side; subject_intervals carries both sides.
                p_lo, p_hi = end_lo, end_hi
                contig_genes = gff_by_contig.get(contig_p, [])
                # in-piece: gene overlaps either interval
                inn  = [g for g in contig_genes
                        if (g["end"]   >= end_lo   and g["start"] <= end_hi) or
                           (g["end"]   >= start_lo and g["start"] <= start_hi)]
                # Origin-crossing GFF feature (end > clen) lands here too
                if clen > 0:
                    for g in contig_genes:
                        if g.get("end", 0) > clen and g not in inn:
                            inn.append(g)
                # Outer flanks live in the contig MIDDLE (between intervals)
                middle = [g for g in contig_genes
                          if g["start"] > start_hi and g["end"] < end_lo]
                # Split the contig middle into each interval's adjacent flanking half.
                pre  = [g for g in middle if g["end"] < end_lo]   # closest to end-side
                post = [g for g in middle if g["start"] > start_hi]  # closest to start-side
                # Take the n genes closest to each piece boundary
                pre_sorted  = sorted(pre,  key=lambda g: -g["end"])    # descending = closest to end_lo first
                post_sorted = sorted(post, key=lambda g: g["start"])   # ascending = closest to start_hi first
                left_n  = list(reversed(pre_sorted[:n_flank]))         # so L1 = left_n[-1] (closest)
                right_n = post_sorted[:n_flank]                        # R1 = right_n[0]   (closest)
            else:
                p_lo = min(int(piece["sstart"]), int(piece["send"]))
                p_hi = max(int(piece["sstart"]), int(piece["send"]))

                # Absorb the full protein footprint when boundary tblastn overlaps the piece.
                for tag in (cl_ref_tag, cr_ref_tag):
                    span = _best_anchor_hit_on(tag, contig_p, min_pid=anchor_min_id)
                    if span and span["hi"] >= p_lo and span["lo"] <= p_hi:
                        p_lo = min(p_lo, span["lo"])
                        p_hi = max(p_hi, span["hi"])

                contig_genes = gff_by_contig.get(contig_p, [])
                pre  = [g for g in contig_genes if g["end"]   <  p_lo]
                post = [g for g in contig_genes if g["start"] >  p_hi]

                inn = [g for g in contig_genes
                       if g["end"] >= p_lo and g["start"] <= p_hi]

                # Flank-anchor hit exclusion: an inner gene overlapping a flank-anchor
                # tblastn hit moves back to pre / post.
                _flank_hit_spans: list[tuple[int, int]] = []
                for _fa in anchors:
                    if _fa.get("role") not in (
                            "flank_L", "flank_L2", "flank_R", "flank_R2"):
                        continue
                    for _fh in acc_hits.get(_fa.get("locus_tag", ""), []):
                        if _fh.get("sseqid") != contig_p:
                            continue
                        try:
                            _fh_lo = min(int(_fh["sstart"]), int(_fh["send"]))
                            _fh_hi = max(int(_fh["sstart"]), int(_fh["send"]))
                        except (KeyError, TypeError, ValueError):
                            continue
                        _flank_hit_spans.append((_fh_lo, _fh_hi))

                if _flank_hit_spans:
                    _piece_mid = (p_lo + p_hi) // 2
                    _new_inn: list[dict] = []
                    for _g in inn:
                        _g_lo, _g_hi = int(_g["start"]), int(_g["end"])
                        # A gene already given a family by the forward override is never evicted.
                        if _g.get("family"):
                            _new_inn.append(_g)
                            continue
                        if any(_g_lo <= _fhi and _g_hi >= _flo
                               for _flo, _fhi in _flank_hit_spans):
                            # Unidentified gene overlapping a flank-anchor hit: re-classify as flank.
                            if _g_hi <= _piece_mid:
                                pre.append(_g)
                            else:
                                post.append(_g)
                        else:
                            _new_inn.append(_g)
                    if len(_new_inn) < len(inn):
                        inn = _new_inn
                        pre.sort(key=lambda x: int(x["start"]))
                        post.sort(key=lambda x: int(x["start"]))

                left_n  = pre[-n_flank:]
                right_n = post[:n_flank]

            for g in left_n:  g["_zone"] = "flank_left"
            for g in right_n: g["_zone"] = "flank_right"
            for g in inn:     g["_zone"] = "cluster"
            for g in (left_n + inn + right_n):
                g["_piece_idx"] = piece_idx
                g["_piece_contig"] = contig_p
                g["_piece_s_lo"] = p_lo
                g["_piece_s_hi"] = p_hi

            layout = left_n + inn + right_n
            if s_dir < 0:
                layout = list(reversed(layout))
                # Reverse-strand piece: swap and reverse the flank buckets to keep index order.
                left_n, right_n = list(reversed(right_n)), list(reversed(left_n))
                for g in left_n:  g["_zone"] = "flank_left"
                for g in right_n: g["_zone"] = "flank_right"

            piece_layouts.append({
                "piece":             piece,
                "contig":            contig_p,
                "p_lo":              p_lo,
                "p_hi":              p_hi,
                "orientation":       s_dir,
                "left_n":            left_n,
                "right_n":           right_n,
                "in_cluster":        inn,
                "layout":            layout,
                "is_remote":         bool(piece.get("_is_remote")),
                # Surface the wrap metadata for gene_state.py's piece-internal check.
                "subject_intervals": ([(int(lo), int(hi))
                                       for lo, hi in (subj_ivs or [])]
                                      if subj_ivs else None),
                "wraps_origin":      wraps,
                "contig_len":        clen,
            })

        res["_piece_count"] = len(piece_layouts)
        res["_piece_layouts_meta"] = [
            {"contig": pl["contig"], "p_lo": pl["p_lo"], "p_hi": pl["p_hi"],
             "orientation": pl["orientation"], "is_remote": pl["is_remote"]}
            for pl in piece_layouts
        ]

        def _boundary_anchor_tblastn_metrics(a: dict) -> dict:
            """Best strong tblastn hit in the reference-derived boundary window."""
            tag = a.get("locus_tag")
            if not tag:
                return {"hit": None}
            aq_lo, aq_hi = _anchor_q_interval(a)
            best = None
            best_score = 0.0

            def _extension_window(pl: dict) -> Optional[tuple[int, int, int, str]]:
                piece = pl["piece"]
                try:
                    pq_lo = min(int(piece["qstart"]), int(piece["qend"]))
                    pq_hi = max(int(piece["qstart"]), int(piece["qend"]))
                    p_lo = int(pl["p_lo"])
                    p_hi = int(pl["p_hi"])
                    orient = int(pl.get("orientation") or 1)
                except (KeyError, TypeError, ValueError):
                    return None

                if aq_hi < pq_lo:
                    ref_extend_bp = max(0, pq_lo - aq_lo)
                    q_side = "left"
                elif aq_lo > pq_hi:
                    ref_extend_bp = max(0, aq_hi - pq_hi)
                    q_side = "right"
                elif aq_lo < pq_lo:
                    ref_extend_bp = max(0, pq_lo - aq_lo)
                    q_side = "left"
                elif aq_hi > pq_hi:
                    ref_extend_bp = max(0, aq_hi - pq_hi)
                    q_side = "right"
                else:
                    return None
                if ref_extend_bp <= 0:
                    return None

                if q_side == "left":
                    if orient >= 0:
                        w_lo, w_hi = p_lo - ref_extend_bp, p_lo - 1
                    else:
                        w_lo, w_hi = p_hi + 1, p_hi + ref_extend_bp
                else:
                    if orient >= 0:
                        w_lo, w_hi = p_hi + 1, p_hi + ref_extend_bp
                    else:
                        w_lo, w_hi = p_lo - ref_extend_bp, p_lo - 1
                return min(w_lo, w_hi), max(w_lo, w_hi), ref_extend_bp, q_side

            for pl_idx, pl in enumerate(piece_layouts):
                ext = _extension_window(pl)
                if not ext:
                    continue
                w_lo, w_hi, ref_extend_bp, q_side = ext
                for h in acc_hits.get(tag, []) or []:
                    if h.get("sseqid") != pl["contig"]:
                        continue
                    try:
                        pid = float(h.get("pident", 0))
                        length = int(h.get("length", 0))
                        qlen = int(h.get("qlen", 0))
                        s_lo = min(int(h["sstart"]), int(h["send"]))
                        s_hi = max(int(h["sstart"]), int(h["send"]))
                    except (KeyError, TypeError, ValueError):
                        continue
                    cov = length / max(1, qlen)
                    if pid < tblastn_min_id or cov < tblastn_intact_cov:
                        continue
                    h_mid = (s_lo + s_hi) / 2.0
                    if not (w_lo <= h_mid <= w_hi):
                        continue
                    score = length * pid
                    if score <= best_score:
                        continue
                    subj_gap = _span_gap(s_lo, s_hi, int(pl["p_lo"]), int(pl["p_hi"]))
                    best_score = score
                    best = {
                        "lo": s_lo,
                        "hi": s_hi,
                        "contig": pl["contig"],
                        "pident": pid,
                        "length": length,
                        "qlen": qlen,
                        "piece_idx": pl_idx,
                        "piece": pl,
                        "outside_piece_gap_bp": subj_gap,
                        "reference_extension_bp": ref_extend_bp,
                        "boundary_side": q_side,
                    }
            if best is None:
                return {"hit": None}
            cov = best["length"] / max(1, best["qlen"])
            return {"cov": min(1.0, cov), "pid": best["pident"], "hit": best}

        # ── Target-centric flank gene identification
        first_pl = piece_layouts[0]
        last_pl  = piece_layouts[-1]
        L1_obj = first_pl["left_n"][-1] if len(first_pl["left_n"]) >= 1 else None
        L2_obj = first_pl["left_n"][-2] if len(first_pl["left_n"]) >= 2 else None
        R1_obj = last_pl["right_n"][0]  if len(last_pl["right_n"]) >= 1 else None
        R2_obj = last_pl["right_n"][1]  if len(last_pl["right_n"]) >= 2 else None

        # ── Reverse flank blast (fallback for unidentified flank genes) --
        _blast_db = (aux_prot_db if (not same_genus and aux_prot_db)
                     else ref_prot_db)
        if _blast_db and tgt_faa_index:
            _seen_lts: set[str] = set()
            _rev_candidates = []
            for _pl in piece_layouts:
                for _g in _pl["left_n"] + _pl["right_n"]:
                    if _g is None or _g.get("family"):
                        continue
                    _lt = _g.get("locus_tag", "")
                    if _lt and _lt in _seen_lts:
                        continue
                    _rev_candidates.append(_g)
                    if _lt:
                        _seen_lts.add(_lt)
            if _rev_candidates:
                _reverse_flank_blast(
                    _rev_candidates,
                    tgt_faa_path, tgt_faa_index, _blast_db,
                    blastp_exe=blastp_exe,
                    min_identity=flank_label_min_id,
                )

        # ── Reverse cluster blastp (fallback for unassigned inner genes) --
        if anchor_prot_db and tgt_faa_index:
            _rev_cluster_cands: list[dict] = []
            _seen_cluster_lts: set[str] = set()
            for _pl in piece_layouts:
                for _g in _pl["in_cluster"]:
                    if _g is None or _g.get("family"):
                        continue
                    _lt = _g.get("locus_tag", "")
                    if not _lt or _lt in _seen_cluster_lts:
                        continue
                    _rev_cluster_cands.append(_g)
                    _seen_cluster_lts.add(_lt)
            if _rev_cluster_cands:
                _reverse_cluster_blast(
                    _rev_cluster_cands,
                    tgt_faa_path, tgt_faa_index, anchor_prot_db,
                    anchor_lt_to_family,
                    blastp_exe=blastp_exe,
                    min_identity=tblastn_min_id,
                )

        # GFF-product fallback for flank genes still unassigned after blastp.
        for _pl in piece_layouts:
            for _g in _pl["left_n"] + _pl["right_n"]:
                if _g is not None and not _g.get("family"):
                    _lbl = _product_to_flank_label(_g.get("product", ""))
                    if _lbl:
                        _g["family"] = _lbl

        def _fmt_flank(g: Optional[dict]) -> str:
            if not g:
                return ""
            fam = g.get("family") or "?"
            ps  = "|P" if g.get("is_pseudo") else ""
            return f"{g.get('locus_tag', '')} [{fam}{ps}] {(g.get('product') or '')[:40]}"

        res["_flank_gene_info"] = {
            "L1": _fmt_flank(L1_obj),
            "L2": _fmt_flank(L2_obj),
            "R1": _fmt_flank(R1_obj),
            "R2": _fmt_flank(R2_obj),
        }

        # ── Per-reference-gene cross-signal state
        ref_inner_genes = [
            a for a in anchors
            if a["role"] in ("cluster_L", "cluster_R", "inner")
        ]
        cl_anchor_for_origin = next(
            (a for a in anchors if a["role"] == "cluster_L"), None)
        origin = int(cl_anchor_for_origin["start"]) if cl_anchor_for_origin else 0
        _cluster_hsps_acc = (cluster_hits or {}).get(acc, []) if cluster_hits else []

        tblastn_intact_id = (tblastn_intact_id_intra if same_genus
                             else tblastn_intact_id_inter)

        gene_states: dict[str, str] = {}
        gene_s_ranges: dict[str, dict] = {}
        gene_signals: dict[str, dict] = {}   # per-anchor (T, D) signal pair
        gene_metrics: dict[str, dict] = {}
        boundary_rescues: dict[str, dict] = {}
        for a in ref_inner_genes:
            lt  = a["locus_tag"]
            rs  = int(a["start"]) - origin
            re_ = int(a["end"])   - origin

            metrics = inner_anchor_tblastn_metrics(
                acc_hits.get(lt, []), piece_layouts,
                min_identity_floor=tblastn_min_id,
            )
            relation = "internal" if metrics.get("hit") else ""
            adj = _boundary_anchor_tblastn_metrics(a)
            if adj.get("hit"):
                cur_hit = metrics.get("hit")
                cur_score = (
                    float(cur_hit.get("length", 0)) * float(cur_hit.get("pident", 0))
                    if cur_hit else 0.0
                )
                adj_hit = adj["hit"]
                adj_score = float(adj_hit.get("length", 0)) * float(adj_hit.get("pident", 0))
                if (not cur_hit) or adj_score > (cur_score * 1.10):
                    metrics = adj
                    relation = "boundary_extension"
                    boundary_rescues[lt] = adj_hit
            dna_metrics = cluster_dna_metrics(
                _cluster_hsps_acc, (rs, re_),
                piece_layouts=piece_layouts,
                threshold=min_gene_cov,
            )
            dna_ok = bool(dna_metrics["present"])
            is_lenient_anchor = bool(a.get("lenient", False))
            eff_intact_cov = aux_lenient_cov if is_lenient_anchor else tblastn_intact_cov
            eff_intact_id = aux_lenient_id if is_lenient_anchor else tblastn_intact_id
            gene_signals[lt] = {"T": bool(metrics["hit"]), "D": bool(dna_ok)}
            gene_metrics[lt] = {
                "cov": metrics["cov"],
                "pid": metrics["pid"],
                "T": bool(metrics["hit"]),
                "D": bool(dna_ok),
                "cluster_dna_cov": dna_metrics["cov"],
                "cluster_dna_intervals": dna_metrics.get("intervals") or [],
                "intact_coverage": eff_intact_cov,
                "intact_identity": eff_intact_id,
                "min_coverage": min_gene_cov,
                "min_identity": tblastn_min_id,
                "lenient": is_lenient_anchor,
                "piece_relation": relation,
                "outside_piece_gap_bp": (
                    metrics["hit"].get("outside_piece_gap_bp", "")
                    if metrics.get("hit") else ""
                ),
                "reference_extension_bp": (
                    metrics["hit"].get("reference_extension_bp", "")
                    if metrics.get("hit") else ""
                ),
                "rescue_reason": (
                    "reference_boundary_extension" if relation == "boundary_extension" else ""
                ),
            }
            gff_pseudo = False
            # Split-CDS evidence demotes the cell to DIVERGENT after classify_gene_state.
            _split_reason = ""
            _force_candidate = False
            if metrics["hit"]:
                contig_genes = gff_by_contig.get(metrics["hit"]["contig"], [])
                g = gff_gene_at(
                    contig_genes,
                    metrics["hit"]["lo"], metrics["hit"]["hi"],
                )
                if g:
                    gff_pseudo = bool(g.get("is_pseudo"))
                    gene_metrics[lt]["target_locus_tag"] = g.get("locus_tag", "")
                    gene_metrics[lt]["target_start"] = g.get("start", "")
                    gene_metrics[lt]["target_end"] = g.get("end", "")
                    gene_metrics[lt]["target_product"] = g.get("product", "")

                # -- Split CDS detection: adjacent pseudo partner, then complementary HSP pair --
                _gene_fragmentary = metrics["cov"] < eff_intact_cov
                if g and not gff_pseudo and _gene_fragmentary:
                    partner_pseudo = adjacent_pseudo_partner(contig_genes, g)
                    if partner_pseudo is not None:
                        _force_candidate = True
                        _split_reason = (
                            f"adjacent_pseudo_partner={partner_pseudo.get('locus_tag','?')}"
                        )
                if g and not gff_pseudo and _gene_fragmentary and not _force_candidate:
                    partner_split = adjacent_split_cds_pair(
                        acc_hits.get(lt, []),
                        contig_genes,
                        metrics["hit"],
                        g,
                    )
                    if partner_split is not None:
                        _force_candidate = True
                        _split_reason = (
                            f"adjacent_split_cds={partner_split.get('locus_tag','?')}"
                        )

                gene_s_ranges[lt] = {
                    "contig": metrics["hit"]["contig"],
                    "s_lo":   metrics["hit"]["lo"],
                    "s_hi":   metrics["hit"]["hi"],
                    "piece_relation": relation,
                    "outside_piece_gap_bp": metrics["hit"].get("outside_piece_gap_bp", ""),
                }
            if _split_reason:
                gene_metrics[lt]["split_cds_reason"] = _split_reason

            gene_states[lt] = classify_gene_state(
                tblastn_cov         = metrics["cov"],
                tblastn_pid         = metrics["pid"],
                cluster_dna_present = dna_ok,
                gff_pseudo          = gff_pseudo,
                intact_coverage     = eff_intact_cov,
                intact_identity     = eff_intact_id,
                min_coverage        = min_gene_cov,
                min_identity        = tblastn_min_id,
                lenient             = is_lenient_anchor,
            )
            if _force_candidate:
                # Split or degraded maps to DIVERGENT, never to PSEUDOGENE on a non-pseudo hit.
                gene_states[lt] = "DIVERGENT"

        # Prevent one target CDS from satisfying multiple local scored anchors.
        local_scored_tags = {
            a.get("locus_tag")
            for a in ref_inner_genes
            if normalize_status_role(
                a.get("status_role", ""),
                a.get("role", ""),
                bool(a.get("exception", False)),
            ) in {"CORE", "ASSOCIATED"}
        }
        by_target_locus: dict[str, list[str]] = {}
        for lt in local_scored_tags:
            if (gene_states.get(lt) or "ABSENT") == "ABSENT":
                continue
            target_lt = str((gene_metrics.get(lt) or {}).get("target_locus_tag") or "")
            if target_lt:
                by_target_locus.setdefault(target_lt, []).append(lt)

        _dup_state_rank = {"INTACT": 3, "DIVERGENT": 2, "PSEUDOGENE": 1, "ABSENT": 0}
        for target_lt, ref_lts in by_target_locus.items():
            if len(ref_lts) < 2:
                continue

            def _dup_score(ref_lt: str) -> tuple[float, float, float, str]:
                m = gene_metrics.get(ref_lt) or {}
                st = gene_states.get(ref_lt, "ABSENT")
                return (
                    float(_dup_state_rank.get(st, 0)),
                    float(m.get("pid") or 0.0),
                    float(m.get("cov") or 0.0),
                    ref_lt,
                )

            winner = max(ref_lts, key=_dup_score)
            for loser in ref_lts:
                if loser == winner:
                    continue
                gene_states[loser] = "ABSENT"
                gene_signals[loser] = {"T": False, "D": False}
                gene_metrics.setdefault(loser, {})["duplicate_suppressed_by"] = winner
                gene_metrics[loser]["duplicate_target_locus_tag"] = target_lt
                gene_metrics[loser]["T"] = False
                gene_metrics[loser]["D"] = False
        res["_tblastn_intact_identity_used"] = tblastn_intact_id
        res["_tblastn_intact_identity_mode"] = (
            "intra-genus" if same_genus else "inter-genus")
        res["_gene_states"]   = gene_states
        res["_gene_s_ranges"] = gene_s_ranges
        res["_gene_metrics"]  = gene_metrics
        res["_boundary_extensions"] = {
            lt: {
                "contig": h.get("contig", ""),
                "s_lo": h.get("lo", ""),
                "s_hi": h.get("hi", ""),
                "piece_idx": h.get("piece_idx", ""),
                "outside_piece_gap_bp": h.get("outside_piece_gap_bp", ""),
                "reference_extension_bp": h.get("reference_extension_bp", ""),
                "boundary_side": h.get("boundary_side", ""),
            }
            for lt, h in boundary_rescues.items()
        }

        # Family-keyed worst-case state (pseudogene signal wins across family)
        _STATE_PRI = {"ABSENT": 0, "INTACT": 1,
                      "DIVERGENT": 2, "PSEUDOGENE": 3}
        ref_family_by_lt = {a["locus_tag"]: a.get("family")
                            for a in ref_inner_genes}
        gene_states_by_family: dict[str, str] = {}
        gene_signals_by_family: dict[str, dict] = {}
        gene_metrics_by_family: dict[str, dict] = {}
        gene_s_ranges_by_family: dict[str, dict] = {}
        gene_relations_by_family: dict[str, dict] = {}
        gene_divergent_reasons_by_family: dict[str, str] = {}
        gene_divergent_details_by_family: dict[str, str] = {}

        def _divergent_reason(m: dict) -> tuple[str, str]:
            cov = float(m.get("cov") or 0.0)
            pid = float(m.get("pid") or 0.0)
            intact_cov = float(m.get("intact_coverage") or 0.70)
            intact_pid = float(m.get("intact_identity") or 50.0)
            min_cov = float(m.get("min_coverage") or 0.30)
            min_pid = float(m.get("min_identity") or 30.0)
            has_d = bool(m.get("D"))
            if cov >= intact_cov and min_pid <= pid < intact_pid:
                return (
                    "pid",
                    f"tblastn identity {pid:.1f}% < intact threshold {intact_pid:.1f}% "
                    f"(coverage {cov:.2f})",
                )
            if min_cov <= cov < intact_cov and pid >= min_pid:
                return (
                    "cov",
                    f"tblastn coverage {cov:.2f} < intact threshold {intact_cov:.2f} "
                    f"(identity {pid:.1f}%)",
                )
            if cov < min_cov and pid >= intact_pid and has_d:
                return (
                    "frag",
                    f"short protein hit coverage {cov:.2f} with DNA support "
                    f"(identity {pid:.1f}%)",
                )
            if cov < min_cov and pid < min_pid and has_d:
                return (
                    "dna_only",
                    f"cluster DNA present but no meaningful protein-level tblastn signal "
                    f"(coverage {cov:.2f}, identity {pid:.1f}%)",
                )
            return (
                "ambig",
                f"divergent signal does not meet intact thresholds "
                f"(coverage {cov:.2f}, identity {pid:.1f}%)",
            )

        for lt, sig in gene_signals.items():
            fam = ref_family_by_lt.get(lt) or lt  # fall back to lt as family key
            fsig = gene_signals_by_family.setdefault(fam, {"T": False, "D": False})
            if sig.get("T"): fsig["T"] = True
            if sig.get("D"): fsig["D"] = True
            m = gene_metrics.get(lt, {})
            prev_m = gene_metrics_by_family.get(fam)
            if prev_m is None or float(m.get("cov") or 0.0) > float(prev_m.get("cov") or 0.0):
                gene_metrics_by_family[fam] = {**m, "ref_locus_tag": lt}
                if lt in gene_s_ranges:
                    gene_s_ranges_by_family[fam] = gene_s_ranges[lt]
                if m.get("piece_relation"):
                    gene_relations_by_family[fam] = {
                        "piece_relation": m.get("piece_relation", ""),
                        "outside_piece_gap_bp": m.get("outside_piece_gap_bp", ""),
                        "reference_extension_bp": m.get("reference_extension_bp", ""),
                        "rescue_reason": m.get("rescue_reason", ""),
                    }
        for lt, st in gene_states.items():
            fam = ref_family_by_lt.get(lt)
            if not fam:
                continue
            prev = gene_states_by_family.get(fam)
            if prev is None or _STATE_PRI[st] > _STATE_PRI[prev]:
                gene_states_by_family[fam] = st
            if st == "DIVERGENT":
                reason, detail = _divergent_reason(gene_metrics.get(lt, {}))
                prev_reason = gene_divergent_reasons_by_family.get(fam)
                if prev_reason is None or _STATE_PRI[st] >= _STATE_PRI.get(prev or "", 0):
                    gene_divergent_reasons_by_family[fam] = reason
                    gene_divergent_details_by_family[fam] = detail
        res["_gene_states_by_family"]  = gene_states_by_family
        res["_gene_signals_by_family"] = gene_signals_by_family
        res["_gene_metrics_by_family"] = gene_metrics_by_family
        res["_gene_s_ranges_by_family"] = gene_s_ranges_by_family
        res["_gene_relations_by_family"] = gene_relations_by_family
        res["_gene_divergent_reasons_by_family"] = gene_divergent_reasons_by_family
        res["_gene_divergent_details_by_family"] = gene_divergent_details_by_family

        # ── Reference-boundary extension cells 
        for a in ref_inner_genes:
            lt = a["locus_tag"]
            adj_hit = boundary_rescues.get(lt)
            if not adj_hit:
                continue
            target_family = a.get("family") or lt
            already_visible = any(
                g.get("family") == target_family
                for pl in piece_layouts
                for g in pl["in_cluster"]
            )
            if already_visible:
                continue
            placed_in = adj_hit["piece"]
            g = gff_gene_at(
                gff_by_contig.get(adj_hit["contig"], []),
                adj_hit["lo"], adj_hit["hi"],
            )
            if g:
                rescue_cell = dict(g)
                rescue_cell["family"] = target_family
            else:
                rescue_cell = {
                    "locus_tag": "",
                    "product": a.get("product", "") or "boundary tblastn signal",
                    "family": target_family,
                    "is_pseudo": False,
                    "start": adj_hit["lo"],
                    "end": adj_hit["hi"],
                    "contig": adj_hit["contig"],
                    "strand": "",
                    "feature_type": "boundary_tblastn_rescue",
                    "attrs": "",
                }
            rescue_cell["_zone"] = "cluster"
            rescue_cell["_piece_idx"] = adj_hit["piece_idx"]
            rescue_cell["_piece_contig"] = placed_in["contig"]
            rescue_cell["_piece_s_lo"] = placed_in["p_lo"]
            rescue_cell["_piece_s_hi"] = placed_in["p_hi"]
            rescue_cell["_piece_relation"] = "boundary_extension"
            rescue_cell["_outside_piece_gap_bp"] = adj_hit.get("outside_piece_gap_bp", "")
            rescue_cell["_reference_extension_bp"] = adj_hit.get("reference_extension_bp", "")
            rescue_cell["_rescue_reason"] = "reference_boundary_extension"

            rescue_lt = rescue_cell.get("locus_tag")
            if rescue_lt:
                for bucket in ("left_n", "right_n", "in_cluster"):
                    placed_in[bucket] = [
                        x for x in placed_in[bucket]
                        if x.get("locus_tag") != rescue_lt
                    ]
            placed_in["in_cluster"].append(rescue_cell)

        # ── DNA-only ABSENT synthetic cells 
        for a in ref_inner_genes:
            lt = a["locus_tag"]
            sig = gene_signals.get(lt, {})
            if gene_states.get(lt) != "ABSENT" or sig.get("T") or not sig.get("D"):
                continue
            # Skip only when the family already has a piece-internal GFF cell.
            target_family = a.get("family") or lt
            already_visible = any(
                g.get("family") == target_family
                for pl in piece_layouts
                for g in pl["in_cluster"]
            )
            if already_visible:
                continue

            # Best signal location in a piece: the tblastn hit footprint, else the
            # cluster blastn HSP overlapping the anchor's reference q-range.
            s_contig, s_lo, s_hi = "", 0, 0
            placed_in: Optional[dict] = None
            tb_hit = gene_s_ranges.get(lt)
            if tb_hit:
                for pl in piece_layouts:
                    if (pl["contig"] == tb_hit["contig"]
                            and pl["p_lo"] <= tb_hit["s_hi"]
                            and pl["p_hi"] >= tb_hit["s_lo"]):
                        s_contig, s_lo, s_hi = (
                            tb_hit["contig"], tb_hit["s_lo"], tb_hit["s_hi"])
                        placed_in = pl
                        break
            if placed_in is None:
                rs  = int(a["start"]) - origin
                re_ = int(a["end"])   - origin
                best_ov = 0
                best_hsp = None
                for h in _cluster_hsps_acc:
                    try:
                        q_lo = min(int(h["qstart"]), int(h["qend"]))
                        q_hi = max(int(h["qstart"]), int(h["qend"]))
                        hs_lo = min(int(h["sstart"]), int(h["send"]))
                        hs_hi = max(int(h["sstart"]), int(h["send"]))
                    except (KeyError, ValueError):
                        continue
                    ov = min(q_hi, re_) - max(q_lo, rs)
                    if ov <= best_ov:
                        continue
                    h_contig = h.get("sseqid", "")
                    for pl in piece_layouts:
                        if (pl["contig"] == h_contig
                                and pl["p_lo"] <= hs_hi
                                and pl["p_hi"] >= hs_lo):
                            best_ov  = ov
                            best_hsp = (h_contig, hs_lo, hs_hi, pl)
                            break
                if best_hsp is None:
                    continue
                s_contig, s_lo, s_hi, placed_in = best_hsp

            placed_in["in_cluster"].append({
                "_zone":        "cluster",
                "locus_tag":    "",
                "product":      a.get("product", "") or "DNA-only signal",
                "family":       a.get("family") or lt,
                "is_pseudo":    False,
                "is_dna_only":  True,
                "start":        s_lo,
                "end":          s_hi,
                "contig":       s_contig,
                "strand":       "",
                "feature_type": "dna_only_absent",
                "attrs":        "",
                "_piece_idx":    piece_layouts.index(placed_in),
                "_piece_contig": placed_in["contig"],
                "_piece_s_lo":   placed_in["p_lo"],
                "_piece_s_hi":   placed_in["p_hi"],
            })

        # ── Domain-guided orphan/split display pieces 
        domain_anchors = [
            a for a in ref_inner_genes
            if pfam_tokens(a.get("pfam", ""))
        ]
        if domain_anchors and domain_hits_by_lt:
            for domain_anchor in domain_anchors:
                target_family = domain_anchor.get("family") or domain_anchor["locus_tag"]
                already_visible = any(
                    g.get("family") == target_family
                    for pl in piece_layouts
                    for g in pl["in_cluster"]
                )
                if already_visible:
                    continue
                required_pfams_ordered = pfam_tokens(domain_anchor.get("pfam", ""))
                required_pfams = set(required_pfams_ordered)
                if not required_pfams:
                    continue
                architecture_label = (
                    domain_anchor.get("family") or domain_anchor["locus_tag"]
                )
                split_enabled = bool(domain_anchor.get("pfam_split"))
                candidates = []
                domains_by_lt = {
                    lt: domain_hits_by_lt.get(lt, set())
                    for lt in domain_lts_by_acc.get(acc, [])
                }
                for lt, lt_domains in domains_by_lt.items():
                    arch_name = domain_architecture(
                        lt_domains, required_pfams, architecture_label)
                    if not arch_name:
                        continue
                    pos = domain_pos_by_lt.get(lt)
                    if not pos:
                        continue
                    contig_genes = _load_contig_genes(pos["contig"])
                    g_obj = next((g for g in contig_genes
                                  if g.get("locus_tag") == lt), None)
                    if not g_obj:
                        continue
                    final_class = architecture_label
                    partner_lt = ""
                    state = "INTACT"
                    reason = "domain_architecture"
                    if arch_name != architecture_label:
                        if not split_enabled:
                            continue
                        prev_g, next_g = same_strand_neighbors(contig_genes, lt)
                        split = split_domain_class(
                            lt_domains, required_pfams, required_pfams_ordered,
                            architecture_label, prev_g, next_g, domains_by_lt,
                        )
                        if not split:
                            continue
                        final_class = f"{split[0]}_{split[2]}"
                        partner_lt = split[1]
                        state = "PSEUDOGENE"
                        reason = split[2]
                    candidates.append((lt, g_obj, final_class, arch_name, partner_lt, state, reason))

                if candidates:
                    # Prefer split-pseudo / split-cds over canonical orphan, then the highest domain count.
                    def _cand_score(item):
                        lt, _g, final_class, _arch, _partner, _state, _reason = item
                        return (
                            1 if "split" in final_class else 0,
                            len(required_pfams & domain_hits_by_lt.get(lt, set())),
                        )

                    lt, g_obj, final_class, arch_name, partner_lt, state, reason = max(
                        candidates, key=_cand_score)
                    contig_id = g_obj["contig"]
                    contig_genes = _load_contig_genes(contig_id)
                    s_lo, s_hi = int(g_obj["start"]), int(g_obj["end"])
                    pre = [g for g in contig_genes if int(g["end"]) < s_lo]
                    post = [g for g in contig_genes if int(g["start"]) > s_hi]
                    left_n = pre[-n_flank:]
                    right_n = post[:n_flank]
                    display_cell = dict(g_obj)
                    display_cell["family"] = target_family
                    display_cell["_zone"] = "cluster"
                    display_cell["_domain_recovered"] = True
                    display_cell["_domain_final_class"] = final_class
                    display_cell["_domain_subtype"] = architecture_label
                    display_cell["_state"] = state
                    display_cell["_decayed_reason"] = "domain_split"
                    display_cell["_decayed_detail"] = (
                        f"{final_class}; partner={partner_lt or 'none'}; "
                        f"domains={','.join(sorted(domain_hits_by_lt.get(lt, set())))}"
                    )
                    piece_idx = len(piece_layouts)
                    for g in left_n:
                        g["_zone"] = "flank_left"
                    for g in right_n:
                        g["_zone"] = "flank_right"
                    for g in left_n + [display_cell] + right_n:
                        g["_piece_idx"] = piece_idx
                        g["_piece_contig"] = contig_id
                        g["_piece_s_lo"] = s_lo
                        g["_piece_s_hi"] = s_hi
                    piece_layouts.append({
                        "piece": None,
                        "contig": contig_id,
                        "p_lo": s_lo,
                        "p_hi": s_hi,
                        "orientation": 1,
                        "left_n": left_n,
                        "right_n": right_n,
                        "in_cluster": [display_cell],
                        "layout": left_n + [display_cell] + right_n,
                        "is_remote": True,
                        "subject_intervals": None,
                        "wraps_origin": False,
                        "contig_len": 0,
                        "display_only_domain": True,
                    })
                    fam_sig = gene_signals_by_family.setdefault(
                        target_family, {"T": False, "D": False})
                    fam_sig["T"] = True
                    gene_states_by_family[target_family] = state
                    gene_states[domain_anchor["locus_tag"]] = state
                    gene_metrics_by_family[target_family] = {
                        "cov": 0.0,
                        "pid": 0.0,
                        "T": True,
                        "D": False,
                        "cluster_dna_cov": 0.0,
                        "piece_relation": "domain_recovered",
                        "rescue_reason": reason,
                        "ref_locus_tag": domain_anchor["locus_tag"],
                        "piece_idx": piece_idx,
                    }
                    gene_relations_by_family[target_family] = {
                        "piece_relation": "domain_recovered",
                        "outside_piece_gap_bp": "",
                        "rescue_reason": reason,
                        "piece_idx": piece_idx,
                    }
                    gene_s_ranges_by_family[target_family] = {
                        "contig": contig_id,
                        "s_lo": s_lo,
                        "s_hi": s_hi,
                        "piece_relation": "domain_recovered",
                        "outside_piece_gap_bp": "",
                        "piece_idx": piece_idx,
                    }
                    if state == "PSEUDOGENE":
                        gene_divergent_reasons_by_family[target_family] = "domain_split"
                        gene_divergent_details_by_family[target_family] = display_cell["_decayed_detail"]
                    domain_recovery_rows.append({
                        "genome_acc": acc,
                        "locus_tag": lt,
                        "family": target_family,
                        "subtype": architecture_label,
                        "state": state,
                        "final_class": final_class,
                        "pfam_arch": arch_name,
                        "domains": ",".join(sorted(domain_hits_by_lt.get(lt, set()))),
                        "partner_locus_tag": partner_lt,
                        "contig": contig_id,
                        "start": s_lo,
                        "end": s_hi,
                        "reason": reason,
                    })

        for pl in piece_layouts:
            _rebuild_layout(pl)

        res["_piece_count"] = len(piece_layouts)
        res["_piece_layouts_meta"] = [
            {"contig": pl["contig"], "p_lo": pl["p_lo"], "p_hi": pl["p_hi"],
             "orientation": pl["orientation"], "is_remote": pl["is_remote"],
             "display_only_domain": bool(pl.get("display_only_domain"))}
            for pl in piece_layouts
        ]

        # Assemble the final deduped list, including synthetic DNA-only cells.
        _SEP_LABELS = {
            "ASSEMBLY_SPLIT_CANDIDATE": "[contig?]",
            "BIOLOGICAL_SPLIT":         "[split]",
            "DIVERGENT_FRAGMENTED":     "[?]",
            "UNCERTAIN":                "[?]",
        }
        _frag_type = res.get("_fragmentation_type", "")
        _multi_contig = int(res.get("contig_count") or 1) > 1

        deduped: list[dict] = []
        for i, pl in enumerate(piece_layouts):
            if i > 0:
                if _multi_contig and _frag_type in _SEP_LABELS:
                    sep_bits = [_SEP_LABELS[_frag_type],
                                f" -> P{i+1}", f"{pl['contig']}:{pl['p_lo']:,}"]
                else:
                    sep_bits = [f"||| -> P{i+1}", f"{pl['contig']}:{pl['p_lo']:,}"]
                if pl["orientation"] < 0:
                    sep_bits.append("(inv)")
                if pl["is_remote"]:
                    sep_bits.append("[remote]")
                deduped.append(_make_separator_cell(
                    " ".join(sep_bits), pl["contig"], pl["p_lo"]))
            deduped.extend(pl["layout"])

        # Correct display-only family labels using the per-anchor best-hit ranges.
        def _range_overlap_len(a_lo: int, a_hi: int, b_lo: int, b_hi: int) -> int:
            return max(0, min(a_hi, b_hi) - max(a_lo, b_lo) + 1)

        def _family_hit_overlap(fam: str, g_obj: dict) -> float:
            r = gene_s_ranges_by_family.get(fam) or {}
            if not r or r.get("contig") != g_obj.get("contig"):
                return 0.0
            try:
                g_lo, g_hi = int(g_obj["start"]), int(g_obj["end"])
                r_lo, r_hi = int(r["s_lo"]), int(r["s_hi"])
            except (KeyError, TypeError, ValueError):
                return 0.0
            hit_len = max(1, r_hi - r_lo + 1)
            return _range_overlap_len(g_lo, g_hi, r_lo, r_hi) / hit_len

        for g in deduped:
            if g.get("_separator") or g.get("_zone") != "cluster":
                continue
            current_fam = g.get("family") or ""
            current_overlap = (
                _family_hit_overlap(current_fam, g)
                if current_fam in gene_s_ranges_by_family else 0.0
            )
            if current_overlap >= 0.50:
                continue
            best_fam = ""
            best_score = 0.0
            for fam, r in gene_s_ranges_by_family.items():
                if r.get("contig") != g.get("contig"):
                    continue
                score = _family_hit_overlap(fam, g)
                if score > best_score:
                    best_fam, best_score = fam, score
            if best_fam and best_score >= 0.50 and best_fam != current_fam:
                g["_family_corrected_from"] = current_fam
                g["family"] = best_fam

        # ── Synteny fingerprints (GFF-family-based) 
        from synteny import row_fingerprints, render_fingerprint

        def _ctx_from_genes(gene_list: list[dict],
                            add_edge: bool = False) -> frozenset:
            # Collect every family token in left_n / right_n.
            labels: set[str] = set()
            for g in gene_list:
                fam = g.get("family")
                if fam:
                    labels.add(_anchor_fam_to_ctx_label.get(fam, fam))
            if add_edge and not gene_list:
                labels.add("EDGE")
            return frozenset(labels) if labels else frozenset({"?"})

        _n_pl = len(piece_layouts)
        _piece_ctx = [
            # EDGE propagated only to the outermost boundary of the first / last piece.
            (_ctx_from_genes(pl["left_n"],
                             i == 0        and not pl.get("wraps_origin")),
             _ctx_from_genes(pl["right_n"],
                             i == _n_pl-1  and not pl.get("wraps_origin")))
            for i, pl in enumerate(piece_layouts)
        ]
        L_FP, R_FP = row_fingerprints(_piece_ctx)
        res["_piece_side_contexts"] = [[sorted(l), sorted(r)] for (l, r) in _piece_ctx]
        res["_L_fingerprint"] = render_fingerprint(L_FP)
        res["_R_fingerprint"] = render_fingerprint(R_FP)

        # ── Per-side flank state

        def _identified(g_obj: Optional[dict]) -> bool:
            return g_obj is not None and bool(g_obj.get("family"))

        def _side_edge_touch(side: str) -> bool:
            # The contig ends on that side when the outermost piece has no GFF genes further out.
            if side == "left":
                return (not piece_layouts[0].get("wraps_origin", False)
                        and not bool(piece_layouts[0]["left_n"]))
            return (not piece_layouts[-1].get("wraps_origin", False)
                    and not bool(piece_layouts[-1]["right_n"]))

        def _side_state(inner_obj: Optional[dict], outer_obj: Optional[dict],
                        side: str) -> str:
            if inner_obj is None:
                return "EDGE" if _side_edge_touch(side) else "MISSING"
            if _identified(inner_obj) and _identified(outer_obj):
                return "FULL"
            if _identified(inner_obj):
                return "HALF"
            return "PRESENT"

        L_state = _side_state(L1_obj, L2_obj, "left")
        R_state = _side_state(R1_obj, R2_obj, "right")
        res["_flank_states"] = {"left": L_state, "right": R_state}

        # ── Synteny classification
        L_full = (L_state == "FULL")
        R_full = (R_state == "FULL")
        if   L_full and R_full: synteny_class = "FULL"
        elif L_full:            synteny_class = "HALF_L"
        elif R_full:            synteny_class = "HALF_R"
        else:                   synteny_class = "NONE"
        res["_synteny_class"] = synteny_class

        # bp measurements from GFF gene positions (intergenic gap between genes)
        def _gene_gap(ga: Optional[dict], gb: Optional[dict]) -> Optional[int]:
            if not ga or not gb:
                return None
            if ga.get("contig") != gb.get("contig"):
                return None
            lo = min(int(ga["end"]), int(gb["end"]))
            hi = max(int(ga["start"]), int(gb["start"]))
            return max(0, hi - lo - 1)

        res["L2_L1_bp"] = _gene_gap(L2_obj, L1_obj) if L_full else None
        res["R1_R2_bp"] = _gene_gap(R1_obj, R2_obj) if R_full else None
        res["flank_bp"]  = _gene_gap(L1_obj, R1_obj) if (L_full and R_full) else None

        cpos_min_c = res.get("_cluster_pos_min")
        cpos_max_c = res.get("_cluster_pos_max")
        if synteny_class == "FULL" and cpos_min_c is not None:
            if (L1_obj and R1_obj
                    and L1_obj.get("contig") == R1_obj.get("contig")):
                l1_end, r1_start = int(L1_obj["end"]), int(R1_obj["start"])
                lo, hi = min(l1_end, r1_start), max(l1_end, r1_start)
                blast_between = (lo <= cpos_min_c and cpos_max_c <= hi)
                res["_synteny_case"] = "CASE_1_AGREE" if blast_between else "CASE_2_DISJOINT"
            else:
                res["_synteny_case"] = "CASE_1_AGREE"
        elif synteny_class in ("HALF_L", "HALF_R"):
            res["_synteny_case"] = "CASE_4_HALF"
        else:
            res["_synteny_case"] = "CASE_3_BLAST_ONLY"

        all_locus_genes[acc] = deduped

        # functional_required feeds the notes column; status comes from status_role.
        func_req_tags = {
            fr["locus_tag"]
            for fr in auto.get("functional_required_resolved", [])
        }
        func_pseudo_lts = [
            g["locus_tag"] for g in deduped
            if g["is_pseudo"] and g.get("locus_tag") in func_req_tags
        ]

        # Strip any previously-appended content notes (idempotent reruns)
        prior = res.get("notes") or ""
        for marker in ("; bp OK but no GFF genes found",
                       "; critical gene(s) pseudo:",
                       "; non-critical gene(s) have pseudo annotation"):
            idx = prior.find(marker)
            if idx >= 0:
                prior = prior[:idx]
        prior = re.sub(r";\s*\d+\s+non-critical gene\(s\) have pseudo annotation",
                       "", prior)
        res["notes"] = prior

        if piece_status == "INTACT":
            if func_pseudo_lts:
                res["notes"]  = ((res.get("notes") or "") +
                                 f"; critical gene(s) pseudo: {func_pseudo_lts}")

        # ── Aux display integration 
        if aux_anchors:
            aux_display_layouts: list[dict] = []

            def _piece_overlaps(pl: dict, contig_id: str, lo: int, hi: int) -> bool:
                if pl["contig"] != contig_id:
                    return False
                intervals = pl.get("subject_intervals") or [(pl["p_lo"], pl["p_hi"])]
                return any(lo <= int(phi) and hi >= int(plo) for plo, phi in intervals)

            def _load_contig_genes(contig_id: str) -> list[dict]:
                if contig_id not in gff_by_contig:
                    raw = parse_gff_region(gff_path, contig_id, 1, 10**9)
                    if rules:
                        classify_genes(raw, rules)
                    seen_lt: dict[str, dict] = {}
                    for gg in raw:
                        lt = gg["locus_tag"]
                        if not lt:
                            continue
                        if lt not in seen_lt:
                            seen_lt[lt] = gg
                        else:
                            existing = seen_lt[lt]
                            if gg["feature_type"] == "CDS" and existing["feature_type"] != "CDS":
                                seen_lt[lt] = gg
                            elif (gg["is_pseudo"] and not existing["is_pseudo"]
                                  and existing["feature_type"] == "gene"):
                                seen_lt[lt] = gg
                    gff_by_contig[contig_id] = sorted(
                        seen_lt.values(), key=lambda gg: gg["start"])
                return gff_by_contig.get(contig_id, [])

            def _piece_gap(pl: dict, contig_id: str, lo: int, hi: int) -> Optional[int]:
                if pl["contig"] != contig_id:
                    return None
                intervals = pl.get("subject_intervals") or [(pl["p_lo"], pl["p_hi"])]
                return min(
                    _span_gap(int(plo), int(phi), lo, hi)
                    for plo, phi in intervals
                )

            def _best_aux_hit(aux_a: dict) -> Optional[dict]:
                aux_lt = aux_a.get("locus_tag", "")
                hits = acc_hits.get(aux_lt, []) or []
                if not hits:
                    return None
                candidates = []
                near_candidates = []
                qlen = aux_qlen.get(aux_lt, 1)
                if aux_a.get("lenient"):
                    remote_min_pid = aux_lenient_id
                    remote_min_cov = aux_lenient_cov
                    near_min_pid = aux_lenient_id
                    near_min_cov = aux_lenient_cov
                else:
                    remote_min_pid = float(cbcfg.get(
                        "single_gene_min_identity_intra" if same_genus
                        else "single_gene_min_identity_inter",
                        90 if same_genus else 75,
                    ))
                    remote_min_cov = float(cbcfg.get("single_gene_min_coverage", 0.70))
                    # Aux hits beside a main piece use the inter-lineage single-gene floor.
                    near_min_pid = float(cbcfg.get("single_gene_min_identity_inter", 75))
                    near_min_cov = remote_min_cov
                span = aux_cohort_spans.get(target_genus, {}).get(aux_lt)

                def _clears(cov_value: float, floor: float, hit: dict) -> bool:
                    """Ordinary coverage, or coverage of what the cohort aligns."""
                    if cov_value >= floor:
                        return True
                    if span is None:
                        return False
                    return cohort_coverage(hit, span) >= floor

                for h in hits:
                    try:
                        pid = float(h.get("pident", 0))
                        qstart, qend = int(h["qstart"]), int(h["qend"])
                        sseqid = h.get("sseqid", "")
                        s_lo = min(int(h["sstart"]), int(h["send"]))
                        s_hi = max(int(h["sstart"]), int(h["send"]))
                    except (KeyError, ValueError, TypeError):
                        continue
                    cov = abs(qend - qstart + 1) / max(1, qlen)
                    hh = dict(h)
                    hh["_aux_cov"] = cov
                    hh["_aux_pid"] = pid
                    if span is not None and cov < min(near_min_cov, remote_min_cov):
                        hh["_aux_cohort_cov"] = cohort_coverage(h, span)
                    piece_gaps = [
                        gap for pl in piece_layouts
                        for gap in [_piece_gap(pl, sseqid, s_lo, s_hi)]
                        if gap is not None
                    ]
                    nearest_gap = min(piece_gaps) if piece_gaps else None
                    # The cohort fallback applies only to hits that sit beside an accepted piece.
                    if nearest_gap is not None and nearest_gap <= adjacent_subj_gap:
                        if pid >= near_min_pid and _clears(cov, near_min_cov, h):
                            hh["_aux_nearest_piece_gap"] = nearest_gap
                            hh["_aux_state"] = "INTACT"
                            near_candidates.append(hh)
                    elif pid >= remote_min_pid and cov >= remote_min_cov:
                        hh["_aux_state"] = "INTACT"
                        candidates.append(hh)
                if near_candidates:
                    return max(
                        near_candidates,
                        key=lambda h: (
                            -int(h.get("_aux_nearest_piece_gap", 10**12)),
                            float(h.get("bitscore", 0)),
                        ),
                    )
                if not candidates:
                    return None
                return max(candidates, key=lambda h: float(h.get("bitscore", 0)))

            def _aux_state(aux_a: dict, pid: float, cov: float) -> str:
                if aux_a.get("lenient"):
                    if pid >= aux_lenient_id and cov >= aux_lenient_cov:
                        return "INTACT"
                else:
                    min_pid = float(cbcfg.get(
                        "single_gene_min_identity_intra" if same_genus
                        else "single_gene_min_identity_inter",
                        90 if same_genus else 75,
                    ))
                    min_cov = float(cbcfg.get("single_gene_min_coverage", 0.70))
                    if pid >= min_pid and cov >= min_cov:
                        return "INTACT"
                if pid >= tblastn_min_id:
                    return "DIVERGENT"
                return "ABSENT"

            def _make_aux_cell(aux_a: dict, hit: dict, piece_idx, p_lo, p_hi) -> dict:
                aux_fam = aux_a.get("family") or aux_a.get("locus_tag", "aux")
                s_lo = min(int(hit["sstart"]), int(hit["send"]))
                s_hi = max(int(hit["sstart"]), int(hit["send"]))
                sseqid = hit.get("sseqid", "")
                g = gff_gene_at(_load_contig_genes(sseqid), s_lo, s_hi)
                if g:
                    cell = dict(g)
                    cell["family"] = aux_fam
                else:
                    cell = {
                        "locus_tag": f"aux:{aux_a.get('locus_tag', '')}",
                        "product": aux_a.get("product", ""),
                        "family": aux_fam,
                        "is_pseudo": False,
                        "is_dna_only": False,
                        "start": s_lo,
                        "end": s_hi,
                        "contig": sseqid,
                        "strand": "+" if int(hit["send"]) >= int(hit["sstart"]) else "-",
                        "feature_type": "aux_tblastn_hit",
                        "attrs": "",
                    }
                cell["_zone"] = "aux"
                cell["_piece_idx"] = piece_idx
                cell["_piece_contig"] = sseqid
                cell["_piece_s_lo"] = p_lo
                cell["_piece_s_hi"] = p_hi
                cell["_piece_relation"] = "aux_internal" if piece_idx != "aux" else "aux_piece"
                cell["_tblastn_cov"] = hit.get("_aux_cov", 0.0)
                cell["_tblastn_pid"] = hit.get("_aux_pid", 0.0)
                cell["_state"] = hit.get("_aux_state") or _aux_state(
                    aux_a, cell["_tblastn_pid"], cell["_tblastn_cov"])
                # Flag calls that only cleared the gate against the cohort consensus.
                if hit.get("_aux_cohort_cov") is not None and \
                        cell["_tblastn_cov"] < aux_lenient_cov:
                    cell["_rescue_reason"] = (
                        f"cohort_span:{hit['_aux_cohort_cov']:.2f}"
                    )
                if cell.get("is_pseudo"):
                    cell["_state"] = "PSEUDOGENE"
                cell["_aux_source"] = aux_a.get("locus_tag", "")
                return cell

            def _integrate_aux_near_piece(aux_a: dict, hit: dict) -> bool:
                """Integrate any aux single-gene recovery into a nearby main piece."""
                aux_fam = aux_a.get("family") or aux_a.get("locus_tag", "aux")
                try:
                    sseqid = hit.get("sseqid", "")
                    s_lo = min(int(hit["sstart"]), int(hit["send"]))
                    s_hi = max(int(hit["sstart"]), int(hit["send"]))
                except (KeyError, ValueError, TypeError):
                    return False
                best_pl_idx = None
                best_gap = None
                for pl_idx, pl in enumerate(piece_layouts):
                    gap = _piece_gap(pl, sseqid, s_lo, s_hi)
                    if gap is None:
                        continue
                    if best_gap is None or gap < best_gap:
                        best_gap = gap
                        best_pl_idx = pl_idx
                if best_pl_idx is None or best_gap is None or best_gap > adjacent_subj_gap:
                    return False

                pl = piece_layouts[best_pl_idx]
                new_lo = min(int(pl["p_lo"]), s_lo)
                new_hi = max(int(pl["p_hi"]), s_hi)
                contig_genes = _load_contig_genes(sseqid)

                # Carry over family labels already assigned to the pre-extension layout.
                old_by_lt = {
                    g.get("locus_tag"): g
                    for bucket in ("left_n", "right_n", "in_cluster")
                    for g in pl[bucket]
                    if g.get("locus_tag")
                }
                inn = [
                    dict(g) for g in contig_genes
                    if g["end"] >= new_lo and g["start"] <= new_hi
                ]
                aux_cell = _make_aux_cell(aux_a, hit, best_pl_idx, new_lo, new_hi)
                aux_lt = aux_cell.get("locus_tag")
                replaced = False
                new_inn: list[dict] = []
                for g in inn:
                    lt = g.get("locus_tag")
                    old = old_by_lt.get(lt)
                    if old and old.get("family"):
                        g["family"] = old.get("family")
                    try:
                        ov = min(int(g["end"]), s_hi) - max(int(g["start"]), s_lo) + 1
                    except (KeyError, ValueError, TypeError):
                        ov = 0
                    if ov > 0 and (lt == aux_lt or not replaced):
                        g.update({
                            "family": aux_fam,
                            "_aux_source": aux_a.get("locus_tag", ""),
                            "_tblastn_cov": hit.get("_aux_cov", 0.0),
                            "_tblastn_pid": hit.get("_aux_pid", 0.0),
                            "_state": hit.get("_aux_state") or _aux_state(
                                aux_a, hit.get("_aux_pid", 0.0), hit.get("_aux_cov", 0.0)),
                            "_piece_relation": "aux_integrated",
                        })
                        # writers.py reads _rescue_reason off the gene, so tag it here.
                        if hit.get("_aux_cohort_cov") is not None and \
                                hit.get("_aux_cov", 0.0) < aux_lenient_cov:
                            g["_rescue_reason"] = (
                                f"cohort_span:{hit['_aux_cohort_cov']:.2f}"
                            )
                        replaced = True
                    new_inn.append(g)
                synthetic_cells = [
                    dict(g) for g in pl.get("in_cluster", [])
                    if not g.get("locus_tag")
                    and (
                        g.get("feature_type") in {
                            "dna_only_decayed",
                            "dna_only_pseudo",
                            "adjacent_tblastn_rescue",
                        }
                        or g.get("is_dna_only")
                    )
                ]
                present_families = {g.get("family") for g in new_inn}
                for g in synthetic_cells:
                    if g.get("family") in present_families:
                        continue
                    try:
                        g_start, g_end = int(g["start"]), int(g["end"])
                    except (KeyError, ValueError, TypeError):
                        continue
                    if g.get("contig") == sseqid and g_end >= new_lo and g_start <= new_hi:
                        new_inn.append(g)
                        present_families.add(g.get("family"))
                if not replaced:
                    aux_cell["_zone"] = "cluster"
                    aux_cell["_piece_relation"] = "aux_integrated"
                    new_inn.append(aux_cell)

                pre = [dict(g) for g in contig_genes if g["end"] < new_lo]
                post = [dict(g) for g in contig_genes if g["start"] > new_hi]
                left_n = pre[-n_flank:]
                right_n = post[:n_flank]
                for g in left_n:
                    g["_zone"] = "flank_left"
                for g in right_n:
                    g["_zone"] = "flank_right"
                for g in new_inn:
                    g["_zone"] = "cluster"
                for g in left_n + new_inn + right_n:
                    g["_piece_idx"] = best_pl_idx
                    g["_piece_contig"] = sseqid
                    g["_piece_s_lo"] = new_lo
                    g["_piece_s_hi"] = new_hi

                pl["p_lo"] = new_lo
                pl["p_hi"] = new_hi
                pl["left_n"] = left_n
                pl["right_n"] = right_n
                pl["in_cluster"] = new_inn
                return True

            for aux_a in aux_anchors:
                hit = _best_aux_hit(aux_a)
                if not hit:
                    continue
                try:
                    sseqid = hit.get("sseqid", "")
                    s_lo = min(int(hit["sstart"]), int(hit["send"]))
                    s_hi = max(int(hit["sstart"]), int(hit["send"]))
                except (KeyError, ValueError, TypeError):
                    continue
                aux_fam = aux_a.get("family") or aux_a.get("locus_tag", "aux")
                state = hit.get("_aux_state") or _aux_state(
                    aux_a, hit.get("_aux_pid", 0.0), hit.get("_aux_cov", 0.0))
                aux_gff = gff_gene_at(_load_contig_genes(sseqid), s_lo, s_hi)
                if aux_gff and aux_gff.get("is_pseudo"):
                    state = "PSEUDOGENE"
                hit["_aux_state"] = state
                if state == "ABSENT":
                    continue

                integrated = _integrate_aux_near_piece(aux_a, hit)
                if integrated:
                    aux_fam = aux_a.get("family") or aux_a.get("locus_tag", "aux")
                    res.setdefault("_aux_hits", {})[aux_fam] = {
                        "state": state,
                        "pid": hit.get("_aux_pid", 0.0),
                        "cov": hit.get("_aux_cov", 0.0),
                        "source_locus": aux_a.get("locus_tag", ""),
                        "hit_contig": sseqid,
                        "hit_lo": s_lo,
                        "hit_hi": s_hi,
                    }
                    continue

                placed = False
                for pl_idx, pl in enumerate(piece_layouts):
                    if not _piece_overlaps(pl, sseqid, s_lo, s_hi):
                        continue
                    cell = _make_aux_cell(aux_a, hit, pl_idx, pl["p_lo"], pl["p_hi"])
                    for bucket in ("left_n", "right_n", "in_cluster"):
                        pl[bucket] = [
                            g for g in pl[bucket]
                            if g.get("locus_tag") != cell.get("locus_tag")
                            and g.get("family") != aux_fam
                        ]
                    pl["in_cluster"].append(cell)
                    placed = True
                    break

                if not placed:
                    contig_genes = _load_contig_genes(sseqid)
                    pre = [g for g in contig_genes if g["end"] < s_lo]
                    post = [g for g in contig_genes if g["start"] > s_hi]
                    left_n = pre[-n_flank:]
                    right_n = post[:n_flank]
                    for g in left_n:
                        g["_zone"] = "flank_left"
                    for g in right_n:
                        g["_zone"] = "flank_right"
                    aux_cell = _make_aux_cell(aux_a, hit, "aux", s_lo, s_hi)
                    p_idx = len(piece_layouts) + len(aux_display_layouts)
                    aux_cell["_piece_idx"] = p_idx
                    for g in left_n + right_n + [aux_cell]:
                        g["_piece_idx"] = p_idx
                        g["_piece_contig"] = sseqid
                        g["_piece_s_lo"] = s_lo
                        g["_piece_s_hi"] = s_hi
                    aux_display_layouts.append({
                        "piece": None,
                        "contig": sseqid,
                        "p_lo": s_lo,
                        "p_hi": s_hi,
                        "orientation": 1,
                        "left_n": left_n,
                        "right_n": right_n,
                        "in_cluster": [aux_cell],
                        "layout": left_n + [aux_cell] + right_n,
                        "is_remote": True,
                        "display_only_aux": True,
                    })

                res.setdefault("_aux_hits", {})[aux_fam] = {
                    "state": state,
                    "pid": hit.get("_aux_pid", 0.0),
                    "cov": hit.get("_aux_cov", 0.0),
                    "source_locus": aux_a.get("locus_tag", ""),
                    "hit_contig": sseqid,
                    "hit_lo": s_lo,
                    "hit_hi": s_hi,
                }

            if any(pl["in_cluster"] for pl in piece_layouts if pl.get("in_cluster")) or aux_display_layouts:
                display_layouts = piece_layouts + aux_display_layouts
                for pl in display_layouts:
                    _rebuild_layout(pl)

                display_deduped: list[dict] = []
                _frag_type = res.get("_fragmentation_type", "")
                _multi_contig = int(res.get("contig_count") or 1) > 1
                for i, pl in enumerate(display_layouts):
                    if i > 0:
                        if pl.get("display_only_aux"):
                            sep_bits = [f"[aux] -> P{i+1}", f"{pl['contig']}:{pl['p_lo']:,}"]
                        elif _multi_contig and _frag_type in _SEP_LABELS:
                            sep_bits = [_SEP_LABELS[_frag_type],
                                        f" -> P{i+1}", f"{pl['contig']}:{pl['p_lo']:,}"]
                        else:
                            sep_bits = [f"||| -> P{i+1}", f"{pl['contig']}:{pl['p_lo']:,}"]
                        if pl["orientation"] < 0:
                            sep_bits.append("(inv)")
                        if pl["is_remote"] and not pl.get("display_only_aux"):
                            sep_bits.append("[remote]")
                        display_deduped.append(_make_separator_cell(
                            " ".join(sep_bits), pl["contig"], pl["p_lo"]))
                    display_deduped.extend(pl["layout"])
                all_locus_genes[acc] = display_deduped

        anchor_states = dict(gene_states)
        for aux_a in aux_anchors:
            aux_fam = aux_a.get("family") or aux_a.get("locus_tag", "aux")
            aux_hit = (res.get("_aux_hits") or {}).get(aux_fam, {})
            anchor_states[aux_a["locus_tag"]] = aux_hit.get("state", "ABSENT")
        if shared_anchors:
            shared_rank = {"ABSENT": 0, "PSEUDOGENE": 1, "DIVERGENT": 2, "INTACT": 3}
            for shared_a in shared_anchors:
                lt = shared_a.get("locus_tag", "")
                fam = shared_a.get("family") or lt or "shared"
                local_state = (anchor_states.get(lt) or "ABSENT").upper()
                remote_hit = _shared_anchor_hit_state(shared_a)
                remote_state = (remote_hit.get("state") or "ABSENT").upper()
                if shared_rank.get(remote_state, 0) >= shared_rank.get(local_state, 0):
                    anchor_states[lt] = remote_state
                    res.setdefault("_shared_hits", {})[fam] = remote_hit
        res["status"], res["_status_detail"] = classify_locus_status(
            anchor_states, anchors)
        res["status"] = apply_locus_coverage_floor(
            res["status"],
            res.get("coverage"),
            min_decayed_coverage,
        )

    if aux_anchors:
        n_aux_total = sum(
            len((r.get("_aux_hits") or {})) for r in ruler_results.values()
        )
        print(f"[content] aux display hits integrated: {n_aux_total}")

    # Output counts must describe the same curated cells used by the final status.
    for res in ruler_results.values():
        detail = res.get("_status_detail")
        if detail:
            res.update(summarize_status_detail(detail))

    # ── Write outputs 
    output_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir = tables_dir(output_dir)
    diag_dir = diagnostics_dir(output_dir)

    write_loci_csv(all_locus_genes, ruler_results, locus_cfg, genome_meta,
                   tbl_dir / "loci.csv")
    write_gene_diagnostics_csv(all_locus_genes, ruler_results, locus_cfg,
                               genome_meta, diag_dir / "gene_diagnostics.csv",
                               tgt_faa_path=tgt_faa_path,
                               tgt_faa_index=tgt_faa_index if tgt_faa_index else None,
                               fna_dir=tgt_cfg.get("fna_dir", ""))
    write_hsp_diagnostics_csv(ruler_results, locus_cfg, genome_meta,
                              diag_dir / "hsp_diagnostics.csv")
    write_domain_recovery_diagnostics(
        domain_recovery_rows,
        diag_dir / "domain_recovery_diagnostics.csv",
    )
    try:
        import openpyxl  # noqa: F401
        write_loci_xlsx(all_locus_genes, ruler_results, locus_cfg, genome_meta,
                        tbl_dir / "loci.xlsx")
    except ImportError:
        print("[content] openpyxl not installed; skipping loci.xlsx "
              "(pip install openpyxl to enable)")

    write_pieces_csv(ruler_results, locus_cfg, genome_meta,
                     tbl_dir / "pieces.csv")
    write_clade_markers_tsv(ruler_results, locus_cfg, genome_meta,
                            tbl_dir / "clade_markers.tsv")
    write_marker_matrix_csv(ruler_results, locus_cfg, genome_meta,
                            tbl_dir / "marker_matrix.csv")
    write_output_guide(output_dir / "OUTPUTS.md")

    return ruler_results
