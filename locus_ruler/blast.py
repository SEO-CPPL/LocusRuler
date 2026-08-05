#!/usr/bin/env python3
"""BLAST database construction and anchor searches."""

import csv, gzip, json, os, re, sqlite3, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ── BLAST output columns ──────────────────────────────────────────────────────
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length",
    "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "qlen", "slen",
]
BLAST_FMT = f"6 {' '.join(BLAST_COLS)}"

# Cluster blastn uses the same column set (nucl vs nucl, same outfmt 6 fields)
CLUSTER_COLS = BLAST_COLS
CLUSTER_FMT  = BLAST_FMT


# ── utilities ──────────────────────────────────────────────────────
def _resolve_anchors_faa(locus_cfg: dict) -> Path:
    """Locate the anchors.faa file from _auto or by conventional naming."""
    stored = locus_cfg.get("_auto", {}).get("_anchors_faa")
    if stored:
        return Path(stored)
    # Fallback: look next to the config file (if _cfg_path was set by the caller)
    cfg_path = locus_cfg.get("_cfg_path")
    if cfg_path:
        p = Path(cfg_path).with_name(
            Path(cfg_path).stem + "_anchors.faa"
        )
        if p.exists():
            return p
    raise FileNotFoundError(
        "anchors.faa not found. Run build_config.py first, or set "
        "_auto._anchors_faa in the locus config."
    )


def _hit_csv_path(work_dir: Path, locus_id: str, accession: str, tag: str) -> Path:
    return work_dir / locus_id / "hits" / f"{accession}__{tag}.tsv"


def _load_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _save_csv(rows: list[dict], path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\t".join(BLAST_COLS) + "\n")
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=BLAST_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _accessions_from_db(db_path: str) -> list[str]:
    con = sqlite3.connect(db_path)
    accs = [r[0] for r in con.execute("SELECT accession FROM genomes").fetchall()]
    con.close()
    return accs


def _fna_path_for(accession: str, fna_dir: str) -> Path | None:
    """Find the genome FASTA (*.fna or *.fna.gz) for an accession."""
    d = Path(fna_dir)
    for suffix in (f"{accession}_genomic.fna.gz",
                   f"{accession}_genomic.fna",
                   f"{accession}.fna.gz",
                   f"{accession}.fna"):
        p = d / suffix
        if p.exists():
            return p
    # Fallback: any file that starts with the accession
    for p in d.iterdir():
        if p.name.startswith(accession) and ".fna" in p.name:
            return p
    return None


# ── 1. Build per-genome BLAST databases ──────────────────────────────────────────────────────
def _make_one_db(args):
    """Worker: makeblastdb for a single genome. Returns (accession, ok, msg)."""
    accession, fna, db_dir, makeblastdb = args
    db_dir = Path(db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = str(db_dir / accession)
    nsq = Path(db_prefix + ".nsq")
    if nsq.exists():
        return accession, True, "exists"

    # Decompress on-the-fly if gzipped
    if str(fna).endswith(".gz"):
        cat_cmd = ["zcat", str(fna)]
        cat = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
        inp = cat.stdout
    else:
        cat = None
        inp = open(fna, "rb")

    cmd = [makeblastdb, "-in", "-", "-dbtype", "nucl",
           "-title", accession, "-out", db_prefix]
    res = subprocess.run(cmd, stdin=inp, capture_output=True)
    if cat:
        cat.wait()
    else:
        inp.close()

    ok = res.returncode == 0
    msg = res.stderr.decode()[:120] if not ok else "ok"
    return accession, ok, msg


def make_dbs(target_cfg: dict, work_dir: Path, cpu: int = 4,
             makeblastdb: str = "makeblastdb") -> dict[str, str]:
    """Build BLAST nucleotide DBs for all genomes in the target."""
    accessions = _accessions_from_db(target_cfg["db"])
    fna_dir    = target_cfg["fna_dir"]
    db_base    = work_dir / "blastdbs"
    db_base.mkdir(parents=True, exist_ok=True)

    tasks = []
    for acc in accessions:
        fna = _fna_path_for(acc, fna_dir)
        if fna is None:
            print(f"  [WARN] no FASTA for {acc} in {fna_dir}, skipping")
            continue
        db_dir = db_base / acc
        tasks.append((acc, str(fna), str(db_dir), makeblastdb))

    results = {}
    done = failed = 0
    with ProcessPoolExecutor(max_workers=cpu) as pool:
        futs = {pool.submit(_make_one_db, t): t[0] for t in tasks}
        for fut in as_completed(futs):
            acc, ok, msg = fut.result()
            if ok:
                results[acc] = str(db_base / acc / acc)
                done += 1
            else:
                print(f"  [ERROR] makeblastdb failed for {acc}: {msg}")
                failed += 1

    print(f"[blast] DBs: {done} built/existing, {failed} failed")
    return results


# ── 2. Anchor tblastn ──────────────────────────────────────────────────────
def _tblastn_one(args) -> tuple[str, str, list[dict]]:
    """Worker: tblastn a single anchor protein against one genome DB."""
    accession, db_prefix, query_faa, tag, evalue, threads, tblastn = args

    cmd = [tblastn,
           "-db",     db_prefix,
           "-query",  query_faa,
           "-out",    "-",
           "-outfmt", BLAST_FMT,
           "-evalue", str(evalue),
           "-num_threads", str(threads),
           "-max_target_seqs", "10"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return accession, tag, []

    _CACHE_MIN_ID = 25.0
    rows = []
    for line in res.stdout.splitlines():
        if not line.strip():
            continue
        vals = line.split("\t")
        if len(vals) < len(BLAST_COLS):
            continue
        row = dict(zip(BLAST_COLS, vals))
        if float(row["pident"]) >= _CACHE_MIN_ID:
            rows.append(row)
    return accession, tag, rows


def run_anchor_blast(locus_cfg: dict, settings: dict,
                     db_map: dict[str, str],
                     work_dir: Path,
                     cpu: int = 4,
                     force: bool = False) -> dict[str, dict[str, list[dict]]]:
    """Run tblastn for every anchor against every genome in db_map."""
    locus_id  = locus_cfg["locus_id"]
    anchors   = locus_cfg["_auto"]["anchors"]
    blast_cfg = {**settings.get("blast", {}), **locus_cfg.get("blast", {})}
    evalue    = blast_cfg.get("evalue", "1e-10")
    tblastn   = settings.get("tools", {}).get("tblastn", "tblastn")
    threads   = int(settings.get("blast", {}).get("threads_per_job", 2))

    anchors_faa = _resolve_anchors_faa(locus_cfg)

    # Write per-anchor temporary FASTAs
    anchor_fastas: dict[str, Path] = {}
    for a in anchors:
        tag = a["locus_tag"]
        role = a["role"]
        tmp = work_dir / locus_id / "query_fastas" / f"{tag}.faa"
        tmp.parent.mkdir(parents=True, exist_ok=True)
        if not tmp.exists() or force:
            # Extract this anchor's sequence from the combined anchors.faa
            seq = _extract_one_seq(anchors_faa, tag)
            if seq:
                tmp.write_text(f">{tag}|{role}\n{seq}\n")
            else:
                print(f"  [WARN] no sequence for {tag} in {anchors_faa}")
                continue
        anchor_fastas[tag] = tmp

    # Build task list (skip if cache exists)
    tasks = []
    for acc, db_prefix in db_map.items():
        for tag, faa in anchor_fastas.items():
            out = _hit_csv_path(work_dir, locus_id, acc, tag)
            if not force and out.exists():
                continue
            tasks.append((acc, db_prefix, str(faa), tag,
                          evalue, threads, tblastn))

    print(f"[blast] {len(tasks)} tblastn jobs queued "
          f"({len(db_map)} genomes × {len(anchor_fastas)} anchors, "
          f"cached={len(db_map)*len(anchor_fastas)-len(tasks)})")

    with ProcessPoolExecutor(max_workers=cpu) as pool:
        futs = {pool.submit(_tblastn_one, t): t for t in tasks}
        for fut in as_completed(futs):
            acc, tag, rows = fut.result()
            _save_csv(rows, _hit_csv_path(work_dir, locus_id, acc, tag))

    # Collect results from cache
    results: dict[str, dict[str, list[dict]]] = {}
    for acc in db_map:
        results[acc] = {}
        for tag in anchor_fastas:
            results[acc][tag] = _load_csv(_hit_csv_path(work_dir, locus_id, acc, tag))

    return results


# ── 3.
# ── 3a. ──────────────────────────────────────────────────────
def _extract_subject_seq(fna_path: Path, contig: str,
                         lo: int, hi: int) -> str | None:
    """Return the nucleotide subsequence [lo, hi] (1-based, inclusive) on
    `contig` from a genome FASTA (.fna / .fna.gz).  Returns None when the
    contig is missing or coordinates are out of range."""
    if not fna_path or not Path(fna_path).exists():
        return None
    opener = gzip.open if str(fna_path).endswith(".gz") else open
    in_target = False
    buf: list[str] = []
    with opener(fna_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if in_target:
                    break
                if line[1:].split()[0] == contig:
                    in_target = True
            elif in_target:
                buf.append(line.strip().upper())
    if not buf:
        return None
    seq = "".join(buf)
    if lo < 1 or hi > len(seq) or lo >= hi:
        return None
    return seq[lo - 1: hi]   # GFF/BLAST are 1-based, inclusive


def piece_pair_identity(
    fna_path: Path,
    piece_a: dict,
    piece_b: dict,
    blastn: str = "blastn",
    min_aligned_bp: int = 200,
) -> float:
    """Compute the best HSP identity between two pieces' subject sequences."""
    def _se(piece):
        s, e = int(piece["sstart"]), int(piece["send"])
        return piece["sseqid"], min(s, e), max(s, e)

    ca, lo_a, hi_a = _se(piece_a)
    cb, lo_b, hi_b = _se(piece_b)
    seq_a = _extract_subject_seq(fna_path, ca, lo_a, hi_a)
    seq_b = _extract_subject_seq(fna_path, cb, lo_b, hi_b)
    if not seq_a or not seq_b:
        return 0.0

    with tempfile.TemporaryDirectory() as td:
        fa_a = Path(td) / "a.fna"
        fa_b = Path(td) / "b.fna"
        fa_a.write_text(f">pieceA\n{seq_a}\n")
        fa_b.write_text(f">pieceB\n{seq_b}\n")
        try:
            res = subprocess.run(
                [blastn, "-query", str(fa_a), "-subject", str(fa_b),
                 "-task", "blastn", "-word_size", "11",
                 "-outfmt", "6 pident length"],
                capture_output=True, text=True, timeout=120,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return 0.0
        if res.returncode != 0:
            return 0.0

        best = 0.0
        for line in res.stdout.splitlines():
            try:
                pid_s, len_s = line.split("\t")[:2]
                pid, alen = float(pid_s), int(len_s)
            except ValueError:
                continue
            if alen >= min_aligned_bp:
                best = max(best, pid)
        return best


# ── 4. Cluster blastn (nucleotide chunk vs target genomes) ──────────────────────────────────────────────────────
def _blastn_one_cluster(args) -> tuple[str, list[dict]]:
    """
    Worker: blastn reference cluster FASTA chunk against one genome DB.
    Uses -task blastn (not megablast) + word_size 7 for sensitivity ~60%+.
    """
    accession, db_prefix, cluster_fna, word_size, perc_id, evalue, threads, blastn = args
    cmd = [
        blastn,
        "-db",            db_prefix,
        "-query",         cluster_fna,
        "-out",           "-",
        "-outfmt",        CLUSTER_FMT,
        "-task",          "blastn",
        "-word_size",     str(word_size),
        "-perc_identity", str(perc_id),
        "-evalue",        str(evalue),
        "-num_threads",   str(threads),
        "-max_target_seqs", "500",   # collect all HSPs across the genome
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return accession, []
    rows = []
    for line in res.stdout.splitlines():
        if not line.strip():
            continue
        vals = line.split("\t")
        if len(vals) < len(CLUSTER_COLS):
            continue
        rows.append(dict(zip(CLUSTER_COLS, vals)))
    return accession, rows


def run_cluster_blastn(
    locus_cfg: dict,
    settings: dict,
    db_map: dict[str, str],
    work_dir: Path,
    cpu: int = 4,
    force: bool = False,
    locus_dir: "Path | None" = None,
) -> dict[str, list[dict]]:
    """
    Run blastn of the reference cluster nucleotide chunk against all target genomes.

    Parameters are read from settings [cluster_blast] block (or locus-level overrides):
      word_size      default 7      (blastn standard; catches ~60%+ identity)
      perc_identity  default 50     (minimum identity filter; noise floor)
      evalue         default 1e-5

    Caches per-genome TSVs at:
      {work_dir}/{locus_id}/cluster_blast/{accession}.tsv

    Returns {accession: [hsp_rows]} for all genomes in db_map.
    """
    locus_id    = locus_cfg["locus_id"]
    cluster_fna = locus_cfg.get("_auto", {}).get("_cluster_fna")
    if not cluster_fna or not Path(cluster_fna).exists():
        print(f"[blast] cluster blastn: _cluster_fna not found; skipping "
              f"(run build_config first, or define cluster_L/cluster_R in config)")
        return {}

    blast_cfg   = {**settings.get("blast", {}), **locus_cfg.get("blast", {})}
    cbcfg       = {**settings.get("cluster_blast", {}),
                   **locus_cfg.get("cluster_blast", {})}
    word_size   = int(cbcfg.get("word_size", 7))
    perc_id     = float(cbcfg.get("perc_identity", 50))
    evalue      = cbcfg.get("evalue", "1e-5")
    blastn      = settings.get("tools", {}).get("blastn", "blastn")
    threads     = int(blast_cfg.get("threads_per_job", 2))

    out_dir = (locus_dir if locus_dir is not None else work_dir / locus_id) / "cluster_blast"
    out_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for acc, db_prefix in db_map.items():
        out = out_dir / f"{acc}.tsv"
        if not force and out.exists():
            continue
        tasks.append((acc, db_prefix, cluster_fna,
                       word_size, perc_id, evalue, threads, blastn))

    cached = len(db_map) - len(tasks)
    print(f"[blast] cluster blastn: {len(tasks)} jobs queued, {cached} cached  "
          f"(word_size={word_size}, perc_id={perc_id}%)")

    with ProcessPoolExecutor(max_workers=cpu) as pool:
        futs = {pool.submit(_blastn_one_cluster, t): t[0] for t in tasks}
        for fut in as_completed(futs):
            acc, rows = fut.result()
            out = out_dir / f"{acc}.tsv"
            if rows:
                with open(out, "w", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=CLUSTER_COLS, delimiter="\t")
                    w.writeheader()
                    w.writerows(rows)
            else:
                out.write_text("\t".join(CLUSTER_COLS) + "\n")

    # Load all results from cache (including previously-cached genomes)
    results: dict[str, list[dict]] = {}
    for acc in db_map:
        out = out_dir / f"{acc}.tsv"
        results[acc] = _load_csv(out)  # reuse existing _load_csv helper
    return results


# ── sequence helper ──────────────────────────────────────────────────────
def _extract_one_seq(faa_path: Path, target_tag: str) -> str | None:
    """Return the protein sequence for target_tag from an anchor .faa file."""
    found = False
    buf = []
    try:
        with open(faa_path) as f:
            for line in f:
                if line.startswith(">"):
                    if found:
                        break
                    parts = line[1:].split("|")
                    # parts[1] is locus_tag in the anchors.faa header format
                    if len(parts) > 1 and parts[1] == target_tag:
                        found = True
                elif found:
                    buf.append(line.strip())
    except FileNotFoundError:
        return None
    return "".join(buf) if buf else None
