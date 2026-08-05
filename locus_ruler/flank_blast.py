#!/usr/bin/env python3
"""Reverse flank-gene blastp utilities."""

import re
import sqlite3
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

# ── GFF-product → normalized flank label ──────────────────────────────────────────────────────
_STRIP_PREFIXES = re.compile(
    r"^(?:putative|bifunctional|multifunctional|probable|predicted|type\s+\S+)\s+",
    re.IGNORECASE,
)
_STRIP_SUFFIXES = [
    re.compile(r"\s+domain[- ]containing\s+protein\s*$",  re.IGNORECASE),
    re.compile(r"\s+(?:family|related|associated)\s+protein\s*$", re.IGNORECASE),
    re.compile(r"\s+protein\s*$",                          re.IGNORECASE),
    re.compile(r"\s+\d+/\d+[^\s]*.*$",                    re.IGNORECASE),
]
_STRIP_LEADING_PROTEIN = re.compile(r"^protein\s+", re.IGNORECASE)


def _product_to_flank_label(product: str) -> Optional[str]:
    """Derive a short, normalized label from a GFF product description."""
    p = (product or "").strip()
    if not p or p.lower() in {"hypothetical protein", "unknown", "uncharacterized protein"}:
        return None

    # 1. DUF numbers: always keep verbatim (e.g. "DUF3060")
    duf = re.match(r"(DUF\d+)", p, re.IGNORECASE)
    if duf:
        return duf.group(1)

    # 2. Strip leading "protein YacL" → "YacL"
    p = _STRIP_LEADING_PROTEIN.sub("", p)

    # 3. Strip leading qualifier words
    p = _STRIP_PREFIXES.sub("", p)

    # 4. Strip trailing generic phrases
    for pat in _STRIP_SUFFIXES:
        p = pat.sub("", p).strip()

    # 5. Take only the first slash-separated token ("DeoR/GlpR…" → "DeoR")
    p = p.split("/")[0].split(",")[0].strip()

    # 6. Collapse whitespace → underscore, cap at 25 chars
    p = re.sub(r"\s+", "_", p)[:25].strip("_")

    return p or None


def _make_faa_index(faa_path: str) -> dict[str, int]:
    """Build ``{locus_tag: byte_offset_of_header}`` for a FAA file."""
    index: dict[str, int] = {}
    try:
        with open(faa_path, "rb") as f:
            pos = 0
            for line in f:
                if line.startswith(b">"):
                    lt = line[1:].split()[0].decode(errors="replace")
                    index[lt] = pos
                pos += len(line)
    except Exception as e:
        print(f"[flank_blast] WARN: could not index FAA {faa_path}: {e}")
    return index


def _get_seq_from_faa(faa_path: str, locus_tag: str,
                      index: dict[str, int]) -> Optional[str]:
    """Retrieve protein sequence for ``locus_tag`` using the byte-offset
    index from :func:`_make_faa_index`.  Returns None if the tag is not
    indexed or the file cannot be read."""
    pos = index.get(locus_tag)
    if pos is None:
        return None
    seq_parts: list[str] = []
    try:
        with open(faa_path, "rb") as f:
            f.seek(pos)
            f.readline()          # skip the header line
            for line in f:
                if line.startswith(b">"):
                    break
                seq_parts.append(line.decode(errors="replace").strip())
    except Exception:
        return None
    return "".join(seq_parts) if seq_parts else None


def _build_ref_prot_db(
    ref_faa: str,
    ref_accession: str,
    db_dir: Path,
    makeblastdb_exe: str = "makeblastdb",
) -> Optional[str]:
    """Build a blastp database from the reference target's FAA."""
    if not ref_faa or not Path(ref_faa).exists():
        return None

    db_dir.mkdir(parents=True, exist_ok=True)
    db_stem = str(db_dir / "ref_prot_db" / ref_accession)
    stamp    = Path(db_stem + ".done")

    if stamp.exists():
        return db_stem

    Path(db_stem).parent.mkdir(parents=True, exist_ok=True)
    faa_out = Path(db_stem + ".faa")

    written = 0
    with open(faa_out, "w") as out_f:
        current_lt: Optional[str] = None
        seq_buf: list[str] = []

        def _flush():
            nonlocal written
            if current_lt is not None and seq_buf:
                out_f.write(f">{current_lt}\n")
                out_f.write("".join(seq_buf) + "\n")
                written += 1

        # Try to filter to locus_tags belonging to ref_accession.
        ref_lts: set[str] = set()
        db_candidates = list(Path(ref_faa).parent.glob("*.db"))
        for db_cand in db_candidates:
            try:
                con = sqlite3.connect(str(db_cand))
                rows = con.execute(
                    "SELECT locus_tag FROM proteins WHERE genome_acc = ?",
                    (ref_accession,),
                ).fetchall()
                con.close()
                if rows:
                    ref_lts = {r[0] for r in rows}
                    break
            except Exception:
                continue

        use_filter = bool(ref_lts)
        try:
            with open(ref_faa) as f:
                for line in f:
                    if line.startswith(">"):
                        _flush()
                        lt = line[1:].strip().split()[0]
                        current_lt = lt if (not use_filter or lt in ref_lts) else None
                        seq_buf = []
                    elif current_lt is not None:
                        seq_buf.append(line)
                _flush()
        except Exception as e:
            print(f"[flank_blast] WARN: error reading reference FAA: {e}")
            return None

    if written == 0:
        # Filter found nothing → use the full FAA directly
        print(f"[flank_blast]   ref_prot_db: {ref_accession} not in DB, "
              f"using full reference FAA ({Path(ref_faa).name})")
        faa_out.unlink(missing_ok=True)
        faa_out = Path(ref_faa)
        db_stem = str(db_dir / "ref_prot_db" / "all_ref")
        stamp   = Path(db_stem + ".done")
        if stamp.exists():
            return db_stem
        Path(db_stem).parent.mkdir(parents=True, exist_ok=True)
    else:
        print(f"[flank_blast]   ref_prot_db: {written} proteins for {ref_accession}")

    result = subprocess.run(
        [makeblastdb_exe, "-in", str(faa_out), "-dbtype", "prot",
         "-out", db_stem, "-parse_seqids"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"[flank_blast] WARN: makeblastdb (ref prot db) failed: "
              f"{result.stderr[:300]}")
        return None

    stamp.touch()
    return db_stem


def _build_anchor_prot_db(
    anchors_faa: str,
    db_dir: Path,
    makeblastdb_exe: str = "makeblastdb",
) -> Optional[str]:
    """Build a blastp database from the locus anchor proteins (_anchors_faa)."""
    if not anchors_faa or not Path(anchors_faa).exists():
        return None

    db_dir.mkdir(parents=True, exist_ok=True)
    db_stem = str(db_dir / "anchor_prot_db" / Path(anchors_faa).stem)
    stamp   = Path(db_stem + ".done")
    if stamp.exists():
        return db_stem

    Path(db_stem).parent.mkdir(parents=True, exist_ok=True)

    # Write a simplified FAA: only ">locus_tag" as header (avoids 50-char limit)
    simplified_faa = Path(db_stem + "_simplified.faa")
    written = 0
    try:
        with open(anchors_faa) as fin, open(simplified_faa, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    # Header format: role|locus_tag|family|product
                    header = line[1:].strip()
                    parts  = header.split("|")
                    lt     = parts[1] if len(parts) >= 2 else header.split()[0]
                    fout.write(f">{lt}\n")
                    written += 1
                else:
                    fout.write(line)
    except Exception as e:
        print(f"[flank_blast] WARN: could not simplify anchor FAA: {e}")
        return None

    if written == 0:
        print("[flank_blast] WARN: anchor FAA is empty")
        simplified_faa.unlink(missing_ok=True)
        return None

    result = subprocess.run(
        [makeblastdb_exe, "-in", str(simplified_faa), "-dbtype", "prot",
         "-out", db_stem],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"[flank_blast] WARN: makeblastdb (anchor prot db) failed: "
              f"{result.stderr[:300]}")
        simplified_faa.unlink(missing_ok=True)
        return None

    stamp.touch()
    return db_stem


def _reverse_cluster_blast(
    genes: list[dict],
    tgt_faa: str,
    tgt_faa_index: dict[str, int],
    anchor_prot_db: str,
    anchor_lt_to_family: dict[str, str],
    blastp_exe: str = "blastp",
    min_identity: float = 30.0,
    evalue: float = 1e-5,
) -> None:
    """For each unassigned cluster gene, blastp its protein against the
    reference anchor proteins and assign the family of the best-hit anchor.

    Targets: genes in the cluster zone whose family is None or empty.
    Mutates genes in place and sets ``_reverse_cluster_blast = True``.

    This handles cases where the forward tblastn hit fell outside the piece
    boundary — the target protein is the query, so piece coordinates are
    irrelevant.  Only runs for genes already placed inside a cluster piece
    by cluster blastn, so the requirement for cluster-membership evidence is
    already satisfied.
    """
    targets = [g for g in genes
               if g is not None and not g.get("family") and g.get("locus_tag")]
    if not targets or not anchor_prot_db:
        return

    seqs: dict[str, str] = {}
    for g in targets:
        lt  = g["locus_tag"]
        seq = _get_seq_from_faa(tgt_faa, lt, tgt_faa_index)
        if seq:
            seqs[lt] = seq

    if not seqs:
        return

    result_map: dict[str, str] = {}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa",
                                     delete=False, dir="/tmp") as tf:
        for lt, seq in seqs.items():
            tf.write(f">{lt}\n{seq}\n")
        tmp_faa = tf.name

    try:
        proc = subprocess.run(
            [blastp_exe,
             "-query",          tmp_faa,
             "-db",             anchor_prot_db,
             "-outfmt",         "6 qseqid sseqid pident",
             "-evalue",         str(evalue),
             "-num_alignments", "1",
             "-num_threads",    "1"],
            capture_output=True, text=True, timeout=120,
        )
        if proc.returncode == 0:
            for line in proc.stdout.strip().splitlines():
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                qlt, ref_lt, pident = parts[0], parts[1], float(parts[2])
                fam = anchor_lt_to_family.get(ref_lt)
                if fam and pident >= min_identity and qlt not in result_map:
                    result_map[qlt] = fam
        elif proc.stderr:
            print(f"[flank_blast] WARN: blastp (reverse cluster) stderr: "
                  f"{proc.stderr[:200]}")
    except subprocess.TimeoutExpired:
        print("[flank_blast] WARN: blastp (reverse cluster) timed out")
    except Exception as e:
        print(f"[flank_blast] WARN: blastp (reverse cluster) error: {e}")
    finally:
        try:
            Path(tmp_faa).unlink()
        except Exception:
            pass

    for g in targets:
        lt = g["locus_tag"]
        if lt in result_map:
            fam = result_map[lt]
            g["family"] = fam
            g["_reverse_cluster_blast"] = True
            print(f"[flank_blast]   reverse cluster blastp: {lt} → {fam}")


def _reverse_flank_blast(
    genes: list[dict],
    tgt_faa: str,
    tgt_faa_index: dict[str, int],
    ref_prot_db: str,
    blastp_exe: str = "blastp",
    min_identity: float = 30.0,
    evalue: float = 1e-5,
) -> None:
    """For each gene with no family, blastp its protein against the
    reference proteome and set ``gene["family"]`` to the best-hit
    reference locus_tag.

    Mutates ``genes`` in place.  Genes that already have a family, or
    whose protein cannot be retrieved from ``tgt_faa``, are skipped.

    A single blastp invocation handles all candidate genes at once
    (-num_alignments 1 returns one best hit per query).  Hits below
    ``min_identity`` are ignored — this is the same 30 % twilight-zone
    floor used elsewhere in the pipeline.
    """
    targets = [g for g in genes
               if g is not None and not g.get("family") and g.get("locus_tag")]
    if not targets or not ref_prot_db:
        return

    seqs: dict[str, str] = {}
    for g in targets:
        lt  = g["locus_tag"]
        seq = _get_seq_from_faa(tgt_faa, lt, tgt_faa_index)
        if seq:
            seqs[lt] = seq

    if not seqs:
        return

    result_map: dict[str, str] = {}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa",
                                     delete=False, dir="/tmp") as tf:
        for lt, seq in seqs.items():
            tf.write(f">{lt}\n{seq}\n")
        tmp_faa = tf.name

    try:
        proc = subprocess.run(
            [blastp_exe,
             "-query",          tmp_faa,
             "-db",             ref_prot_db,
             "-outfmt",         "6 qseqid sseqid pident",
             "-evalue",         str(evalue),
             "-num_alignments", "1",
             "-num_threads",    "1"],
            capture_output=True, text=True, timeout=60,
        )
        if proc.returncode == 0:
            for line in proc.stdout.strip().splitlines():
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                qlt, ref_lt, pident = parts[0], parts[1], float(parts[2])
                if pident >= min_identity and qlt not in result_map:
                    result_map[qlt] = ref_lt   # first hit = best-scoring
        elif proc.stderr:
            print(f"[flank_blast] WARN: blastp stderr: {proc.stderr[:200]}")
    except subprocess.TimeoutExpired:
        print("[flank_blast] WARN: blastp (reverse flank) timed out")
    except Exception as e:
        print(f"[flank_blast] WARN: blastp (reverse flank) error: {e}")
    finally:
        try:
            Path(tmp_faa).unlink()
        except Exception:
            pass

    for g in targets:
        lt = g["locus_tag"]
        if lt in result_map:
            g["family"] = result_map[lt]


def _select_aux_genome(ruler_results: dict) -> Optional[str]:
    """Pick the median-coverage target genome for use as an auxiliary proteome."""
    GOOD = {"INTACT", "FRAGMENTED", "PSEUDOGENIZED", "ANNOTATION_GAP",
            "PARTIAL_DEL", "LARGE_DEL", "DUPLICATED"}
    cands: list[tuple[float, str]] = []
    for acc, res in ruler_results.items():
        status = res.get("status", "")
        if status not in GOOD:
            continue
        cov = float(res.get("coverage") or 0.0)
        if cov < 0.5:
            continue
        cands.append((cov, acc))
    if not cands:
        return None
    cands.sort()                        # ascending by coverage
    median_idx = len(cands) // 2       # 50th-percentile index (lower bound)
    return cands[median_idx][1]


def _build_aux_prot_db(
    aux_acc: str,
    db_path: str,
    faa_path: str,
    work_dir: Path,
    makeblastdb_exe: str = "makeblastdb",
) -> Optional[str]:
    """Build a blastp DB from the proteins of a single target-genus genome."""
    if not faa_path or not Path(faa_path).exists():
        return None
    if not db_path or not Path(db_path).exists():
        return None

    out_dir = work_dir / "aux_prot_db" / aux_acc
    out_dir.mkdir(parents=True, exist_ok=True)
    db_stem = str(out_dir / "aux_prot")
    stamp   = Path(db_stem + ".done")
    if stamp.exists():
        return db_stem

    # Get locus_tags for aux_acc from the target DB
    aux_lts: set[str] = set()
    try:
        con  = sqlite3.connect(db_path)
        rows = con.execute(
            "SELECT locus_tag FROM proteins WHERE genome_acc = ?",
            (aux_acc,),
        ).fetchall()
        con.close()
        aux_lts = {r[0] for r in rows}
    except Exception as e:
        print(f"[flank_blast] WARN: could not query aux genome {aux_acc}: {e}")
        return None

    if not aux_lts:
        print(f"[flank_blast] WARN: no proteins found for aux genome {aux_acc}")
        return None

    faa_out = Path(db_stem + ".faa")
    written = 0
    try:
        with open(faa_path) as fin, open(faa_out, "w") as fout:
            current_lt: Optional[str] = None
            buf: list[str] = []

            def _flush_aux() -> None:
                nonlocal written
                if current_lt and buf:
                    fout.write(f">{current_lt}\n{''.join(buf)}\n")
                    written += 1

            for line in fin:
                if line.startswith(">"):
                    _flush_aux()
                    lt = line[1:].strip().split()[0]
                    current_lt = lt if lt in aux_lts else None
                    buf = []
                elif current_lt is not None:
                    buf.append(line)
            _flush_aux()
    except Exception as e:
        print(f"[flank_blast] WARN: error writing aux FAA for {aux_acc}: {e}")
        return None

    if written == 0:
        print(f"[flank_blast] WARN: 0 proteins written for aux genome {aux_acc}")
        faa_out.unlink(missing_ok=True)
        return None

    print(f"[flank_blast]   aux_prot_db: {written} proteins for {aux_acc}")
    result = subprocess.run(
        [makeblastdb_exe, "-in", str(faa_out), "-dbtype", "prot",
         "-out", db_stem, "-parse_seqids"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"[flank_blast] WARN: makeblastdb (aux prot db) failed: "
              f"{result.stderr[:300]}")
        return None

    stamp.touch()
    return db_stem
