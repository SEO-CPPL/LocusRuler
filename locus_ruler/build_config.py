#!/usr/bin/env python3
"""Reference cluster extraction into the locus config."""

import argparse, csv, gzip, json, re, sqlite3, sys
from datetime import datetime
from pathlib import Path

_PKG_DIR = Path(__file__).resolve().parent
if str(_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(_PKG_DIR))

from config_utils import load_settings, load_locus_cfg, save_locus_cfg, get_target_cfg
from classify import load_rules, classify_product
from domain_recovery import pfam_tokens
from locus_status import VALID_STATUS_ROLES, normalize_status_role


def check_anchors_csv(path: Path) -> list[str]:
    """Problems that would make an edited anchor table quietly do nothing."""
    try:
        with open(path, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            rows = list(reader)
    except UnicodeDecodeError:
        return [f"{path.name} is not UTF-8; save it as UTF-8"]
    except OSError as exc:
        return [f"cannot read {path.name}: {exc}"]

    missing = [name for name in ("locus_tag", "role", "status_role")
               if name not in fields]
    if missing:
        return [f"missing column {name!r}" for name in missing]
    if not rows:
        return ["the table has no rows left"]

    allowed = ", ".join(sorted(VALID_STATUS_ROLES))
    problems: list[str] = []
    for number, row in enumerate(rows, start=2):
        if not (row.get("locus_tag") or "").strip():
            problems.append(f"line {number}: no locus_tag, so the row is ignored")
            continue
        value = (row.get("status_role") or "").strip().upper()
        if value and value not in VALID_STATUS_ROLES:
            problems.append(
                f"line {number}: status_role {value!r} is not one of {allowed}")
        for flag in ("exception", "lenient", "pfam_split"):
            raw = (row.get(flag) or "").strip().upper()
            if raw and raw not in ("TRUE", "FALSE"):
                # Anything other than TRUE reads as FALSE, so a typo silently disables it.
                problems.append(
                    f"line {number}: {flag} is {raw!r}, which counts as FALSE; "
                    "write TRUE or FALSE")
        pfam = (row.get("pfam") or "").strip()
        if pfam and not pfam_tokens(pfam):
            problems.append(
                f"line {number}: pfam {pfam!r} holds no PF accession, so no "
                "profile is searched for; they read PF and five digits")
        if (row.get("pfam_split") or "").strip().upper() == "TRUE" and not pfam:
            problems.append(
                f"line {number}: pfam_split is TRUE but pfam is empty, "
                "so there is nothing to split")
    return problems


def apply_existing_family_overrides(genes: list[dict], cfg_path: Path) -> int:
    """Read existing _anchors.csv and merge curated family/exception labels into genes."""
    anchors_csv = cfg_path.with_name(cfg_path.stem + "_anchors.csv")
    if not anchors_csv.exists():
        return 0
    existing: dict[str, dict] = {}
    try:
        with open(anchors_csv, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                lt = row.get("locus_tag")
                if lt:
                    d = {
                        "family":    row.get("family", ""),
                        "exception": row.get("exception", "").upper() == "TRUE",
                        "lenient":   row.get("lenient", "").upper() == "TRUE",
                        "status_role": normalize_status_role(
                            row.get("status_role", ""),
                            row.get("role", ""),
                            row.get("exception", "").upper() == "TRUE",
                        ),
                    }
                    if "pfam" in row:
                        d["pfam"] = row.get("pfam", "")
                    if "pfam_split" in row:
                        d["pfam_split"] = row.get("pfam_split", "").upper() == "TRUE"
                    existing[lt] = d
    except ValueError:
        raise
    except Exception as exc:
        # The table exists but is unreadable, so its curated values are discarded.
        print(f"[build_config] WARNING: could not read {anchors_csv.name} "
              f"({exc}); curated edits in it are being ignored")
        return 0

    if not existing:
        print(f"[build_config] WARNING: {anchors_csv.name} has no readable "
              f"locus_tag column; curated edits in it are being ignored")

    n = 0
    for g in genes:
        lt = g.get("locus_tag", "")
        if lt in existing:
            d = existing[lt]
            g["family"]    = d["family"]
            g["exception"] = d["exception"]
            g["lenient"]   = d["lenient"]
            g["status_role"] = d["status_role"]
            for col in ("pfam", "pfam_split"):
                if col in d:
                    g[col] = d[col]
            n += 1
    if n:
        print(f"[build_config] preserved {n} family/status_role/exception/lenient edit(s) "
              f"from existing {anchors_csv.name}")
    return n


def read_aux_rows(cfg_path: Path) -> list[dict]:
    """Read role='aux' rows from existing _anchors.csv."""
    anchors_csv = cfg_path.with_name(cfg_path.stem + "_anchors.csv")
    if not anchors_csv.exists():
        return []
    out: list[dict] = []
    try:
        with open(anchors_csv, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if (row.get("role", "").strip().lower() == "aux"
                        and row.get("locus_tag", "").strip()):
                    d = {
                        "role":      "aux",
                        "locus_tag": row["locus_tag"].strip(),
                        "family":    row.get("family", "").strip(),
                        "exception": row.get("exception", "").upper() == "TRUE",
                        "lenient":   row.get("lenient", "").upper() == "TRUE",
                        "status_role": normalize_status_role(
                            row.get("status_role", ""),
                            "aux",
                            row.get("exception", "").upper() == "TRUE",
                        ),
                    }
                    if "pfam" in row:
                        d["pfam"] = row.get("pfam", "").strip()
                    if "pfam_split" in row:
                        d["pfam_split"] = row.get("pfam_split", "").upper() == "TRUE"
                    out.append(d)
    except ValueError:
        raise
    except Exception as e:
        print(f"[WARN] failed to read aux rows from {anchors_csv.name}: {e}")
        return []
    return out


def resolve_aux_rows(con: sqlite3.Connection, aux_rows: list[dict]) -> list[dict]:
    """Look up each aux locus_tag across all genomes and fill in coordinates."""
    resolved: list[dict] = []
    for r in aux_rows:
        info = fetch_gene_anywhere(con, r["locus_tag"])
        if not info:
            print(f"[build_config] [WARN] aux locus_tag '{r['locus_tag']}' "
                  f"not found in DB — skipping")
            continue
        g = dict(info)            # locus_tag, genome_acc, contig, start, end, strand, product
        g.update({
            "role":      "aux",
            "family":    r["family"],
            "exception": r["exception"],
            "lenient":   r["lenient"],
            "status_role": r["status_role"],
            "pfam":       r.get("pfam", ""),
            "pfam_split": r.get("pfam_split", False),
            "inner_idx": None,    # aux has no progressive-round index
        })
        resolved.append(g)
        print(f"[build_config] aux: {r['locus_tag']} → {info['genome_acc']} "
              f"({info['contig']}:{info['start']}-{info['end']}) "
              f"family={r['family'] or '—'}")
    return resolved


def write_anchors_csv(genes: list[dict], cfg_path: Path) -> Path:
    """Write the editable ``<config_stem>_anchors.csv``."""
    anchors_csv = cfg_path.with_name(cfg_path.stem + "_anchors.csv")
    fields = ["role", "locus_tag", "family", "status_role", "exception", "lenient",
              "pfam", "pfam_split", "product",
              "contig", "start", "end", "strand"]
    with open(anchors_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        # Reference-genome rows first (preserve relative order), aux rows last.
        cluster_rows = [g for g in genes if g.get("role") != "aux"]
        aux_rows     = [g for g in genes if g.get("role") == "aux"]
        for g in cluster_rows + aux_rows:
            row = {k: g.get(k, "") for k in fields}
            # Use uppercase string for TRUE/FALSE in CSV for legibility
            row["exception"] = "TRUE" if g.get("exception") else "FALSE"
            row["lenient"]   = "TRUE" if g.get("lenient") else "FALSE"
            row["pfam_split"] = "TRUE" if g.get("pfam_split") else "FALSE"
            g["status_role"] = normalize_status_role(
                g.get("status_role", ""), g.get("role", ""), bool(g.get("exception")))
            row["status_role"] = g["status_role"]
            writer.writerow(row)
    return anchors_csv


# ── DB queries ──────────────────────────────────────────────────────
def fetch_gene_pos(con: sqlite3.Connection, accession: str, locus_tag: str) -> dict | None:
    row = con.execute(
        "SELECT locus_tag, contig, start, end, strand, product "
        "FROM proteins WHERE genome_acc=? AND locus_tag=?",
        (accession, locus_tag)
    ).fetchone()
    return dict(row) if row else None


def fetch_gene_anywhere(con: sqlite3.Connection, locus_tag: str) -> dict | None:
    """Look up a locus_tag across ALL genomes in the DB."""
    row = con.execute(
        "SELECT locus_tag, genome_acc, contig, start, end, strand, product "
        "FROM proteins WHERE locus_tag=? LIMIT 1",
        (locus_tag,)
    ).fetchone()
    return dict(row) if row else None


def fetch_genes_between(con: sqlite3.Connection, accession: str,
                        contig: str, start: int, end: int) -> list[dict]:
    """All CDS whose entire body lies within [start, end], ordered by start."""
    rows = con.execute(
        "SELECT locus_tag, contig, start, end, strand, product "
        "FROM proteins "
        "WHERE genome_acc=? AND contig=? AND start>=? AND end<=? "
        "ORDER BY start",
        (accession, contig, start, end)
    ).fetchall()
    return [dict(r) for r in rows]


def fetch_adjacent_gene(con: sqlite3.Connection, accession: str,
                        contig: str, position: int, direction: str) -> dict | None:
    """Fetch the CDS immediately adjacent to `position` on the given contig."""
    if direction == "left":
        row = con.execute(
            "SELECT locus_tag, contig, start, end, strand, product "
            "FROM proteins WHERE genome_acc=? AND contig=? AND end < ? "
            "ORDER BY end DESC LIMIT 1",
            (accession, contig, position),
        ).fetchone()
    else:
        row = con.execute(
            "SELECT locus_tag, contig, start, end, strand, product "
            "FROM proteins WHERE genome_acc=? AND contig=? AND start > ? "
            "ORDER BY start ASC LIMIT 1",
            (accession, contig, position),
        ).fetchone()
    return dict(row) if row else None


# ── genome FASTA helpers ──────────────────────────────────────────────────────
def _find_fna(accession: str, fna_dir: str) -> Path | None:
    """Locate genome FASTA (.fna or .fna.gz) for a given accession."""
    d = Path(fna_dir)
    for suffix in (f"{accession}_genomic.fna.gz",
                   f"{accession}_genomic.fna",
                   f"{accession}.fna.gz",
                   f"{accession}.fna"):
        p = d / suffix
        if p.exists():
            return p
    for p in d.iterdir():
        if p.name.startswith(accession) and ".fna" in p.name:
            return p
    return None


def extract_cluster_fna(fna_dir: str, accession: str,
                         contig: str, start: int, end: int) -> str | None:
    """Extract the nucleotide subsequence [start:end] from accession's genome FASTA."""
    fna = _find_fna(accession, fna_dir)
    if fna is None:
        return None
    open_fn = gzip.open if str(fna).endswith(".gz") else open
    in_target = False
    buf: list[str] = []
    with open_fn(fna, "rt") as f:
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
    return seq[start:end]


# ── protein sequence extraction ──────────────────────────────────────────────────────
def extract_sequences(faa_path: str, wanted: set[str]) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur_id, buf = None, []
    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id and cur_id in wanted:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].split()[0]
                buf = []
            elif cur_id in wanted:
                buf.append(line.strip())
    if cur_id and cur_id in wanted:
        seqs[cur_id] = "".join(buf)
    return seqs


# ── anchor role assignment ──────────────────────────────────────────────────────
def assign_roles(genes: list[dict],
                 flank_L_tag: str | None,
                 flank_R_tag: str | None,
                 flank_L2_tag: str | None = None,
                 flank_R2_tag: str | None = None) -> list[dict]:
    """Assign 'role' to each gene in the ordered list."""
    # ── Tier 1: no flanks 
    if flank_L_tag is None or flank_R_tag is None:
        inner_idx = 0
        for i, g in enumerate(genes):
            if i == 0:
                g["role"] = "cluster_L"
            elif i == len(genes) - 1:
                g["role"] = "cluster_R"
            else:
                g["role"] = "inner"
                g["inner_idx"] = inner_idx
                inner_idx += 1
        return genes

    # ── Standard: flanks present 
    tags = [g["locus_tag"] for g in genes]

    def find_idx(tag: str | None, label: str) -> int | None:
        if tag is None:
            return None
        try:
            return tags.index(tag)
        except ValueError:
            sys.exit(f"ERROR: {label} '{tag}' not found in extracted gene list.")

    i_L2 = find_idx(flank_L2_tag, "flank_L2")
    i_L  = find_idx(flank_L_tag,  "flank_L")
    i_R  = find_idx(flank_R_tag,  "flank_R")
    i_R2 = find_idx(flank_R2_tag, "flank_R2")

    inner_idx = 0
    for i, g in enumerate(genes):
        if i == i_L2:
            g["role"] = "flank_L2"
        elif i == i_L:
            g["role"] = "flank_L"
        elif i == i_R:
            g["role"] = "flank_R"
        elif i == i_R2:
            g["role"] = "flank_R2"
        elif i_L is not None and i == i_L + 1:
            g["role"] = "cluster_L"
        elif i_R is not None and i == i_R - 1:
            g["role"] = "cluster_R"
        else:
            g["role"] = "inner"
            g["inner_idx"] = inner_idx
            inner_idx += 1
    return genes


def build_round_pairs(genes: list[dict]) -> list[dict]:
    """Generate the ordered list of anchor pairs for progressive bp measurement."""
    by_role = {g["role"]: g for g in genes if g["role"] != "inner"}
    inners  = [g for g in genes if g["role"] == "inner"]

    pairs = []
    # Round 0: flank pair — only when both flanks are present
    if "flank_L" in by_role and "flank_R" in by_role:
        pairs.append({
            "round":     0,
            "left_tag":  by_role["flank_L"]["locus_tag"],
            "right_tag": by_role["flank_R"]["locus_tag"],
        })
    if "cluster_L" in by_role and "cluster_R" in by_role:
        pairs.append({
            "round": 1,
            "left_tag":  by_role["cluster_L"]["locus_tag"],
            "right_tag": by_role["cluster_R"]["locus_tag"],
        })

    if inners:
        li, ri = 0, len(inners) - 1
        r = 2
        left_tag  = by_role["cluster_L"]["locus_tag"]
        right_tag = by_role["cluster_R"]["locus_tag"]
        move_left = True
        while li <= ri:
            if move_left:
                left_tag = inners[li]["locus_tag"]
                li += 1
            else:
                right_tag = inners[ri]["locus_tag"]
                ri -= 1
            pairs.append({"round": r, "left_tag": left_tag, "right_tag": right_tag})
            r += 1
            move_left = not move_left

    return pairs


# ── main ──────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Auto-populate _auto block of a LocusRuler locus config.")
    ap.add_argument("--config",   required=True,  help="Path to locus JSON config")
    ap.add_argument("--settings", required=True,  help="Path to settings.toml")
    ap.add_argument("--dry-run",  action="store_true",
                    help="Print extracted table; do not write files")
    args = ap.parse_args()

    cfg_path  = Path(args.config).resolve()
    settings  = load_settings(Path(args.settings).resolve())
    locus_cfg = load_locus_cfg(cfg_path)

    ref       = locus_cfg["reference"]
    accession = ref["accession"]
    target    = ref["target"]
    flank_L   = ref.get("flank_L")
    flank_R   = ref.get("flank_R")
    func_req  = locus_cfg.get("functional_required", [])

    # Validate flank specification: both or neither
    if bool(flank_L) != bool(flank_R):
        sys.exit("ERROR: specify both flank_L and flank_R, or neither "
                 "(Tier 1 / flank-less mode).")
    has_flanks = bool(flank_L and flank_R)

    # Tier 1: user must supply cluster terminals directly
    cl_direct = ref.get("cluster_L")
    cr_direct = ref.get("cluster_R")
    if not has_flanks and (not cl_direct or not cr_direct):
        sys.exit("ERROR: flank-less mode (Tier 1) requires reference.cluster_L "
                 "and reference.cluster_R — the first and last biosynthetic genes "
                 "of the cluster.")

    tgt_cfg   = get_target_cfg(settings, target)
    rules_file= settings.get("classification", {}).get("rules_file", "")
    rules     = load_rules(rules_file,
                           locus_id   = locus_cfg.get("locus_id"),
                           config_dir = cfg_path.parent)

    print(f"[build_config] accession : {accession}")
    print(f"[build_config] target    : {target}  →  {tgt_cfg['db']}")
    if has_flanks:
        print(f"[build_config] flank_L   : {flank_L}")
        print(f"[build_config] flank_R   : {flank_R}")
    else:
        print(f"[build_config] mode      : Tier 1 (flank-less)")
        print(f"[build_config] cluster_L : {cl_direct}")
        print(f"[build_config] cluster_R : {cr_direct}")

    # ── DB lookup
    con = sqlite3.connect(tgt_cfg["db"])
    con.row_factory = sqlite3.Row

    if has_flanks:
        g_L = fetch_gene_pos(con, accession, flank_L)
        g_R = fetch_gene_pos(con, accession, flank_R)
        if not g_L:
            sys.exit(f"ERROR: flank_L '{flank_L}' not found in DB for {accession}.")
        if not g_R:
            sys.exit(f"ERROR: flank_R '{flank_R}' not found in DB for {accession}.")
        if g_L["contig"] != g_R["contig"]:
            sys.exit(f"ERROR: flank_L ({g_L['contig']}) and flank_R ({g_R['contig']}) "
                     f"are on different contigs. Cannot define a single locus.")
        contig    = g_L["contig"]
        coord_min = min(g_L["start"], g_R["start"])
        coord_max = max(g_L["end"],   g_R["end"])
    else:
        g_CL = fetch_gene_pos(con, accession, cl_direct)
        g_CR = fetch_gene_pos(con, accession, cr_direct)
        if not g_CL:
            sys.exit(f"ERROR: cluster_L '{cl_direct}' not found in DB for {accession}.")
        if not g_CR:
            sys.exit(f"ERROR: cluster_R '{cr_direct}' not found in DB for {accession}.")
        if g_CL["contig"] != g_CR["contig"]:
            sys.exit(f"ERROR: cluster_L ({g_CL['contig']}) and cluster_R "
                     f"({g_CR['contig']}) are on different contigs.")
        contig    = g_CL["contig"]
        coord_min = min(g_CL["start"], g_CR["start"])
        coord_max = max(g_CL["end"],   g_CR["end"])

    genes_raw = fetch_genes_between(con, accession, contig, coord_min, coord_max)

    # ── Fetch L2 / R2 (outer 2nd flanks)
    g_L2 = g_R2 = None
    if has_flanks:
        g_L2 = fetch_adjacent_gene(con, accession, contig, g_L["start"], "left")
        g_R2 = fetch_adjacent_gene(con, accession, contig, g_R["end"],   "right")

    con.close()

    if len(genes_raw) < 2:
        sys.exit("ERROR: fewer than 2 genes found between the specified locus_tags.")

    # ── Classify families
    for g in genes_raw:
        g["family"] = classify_product(g["product"], rules)
    for g in [g_L2, g_R2]:
        if g:
            g["family"] = classify_product(g["product"], rules)

    # ── Build full gene list: [L2] + [L1..R1] + [R2]
    genes_full = []
    if g_L2:
        genes_full.append(g_L2)
    genes_full.extend(genes_raw)
    if g_R2:
        genes_full.append(g_R2)

    genes = assign_roles(
        genes_full,
        flank_L_tag  = flank_L,
        flank_R_tag  = flank_R,
        flank_L2_tag = g_L2["locus_tag"] if g_L2 else None,
        flank_R2_tag = g_R2["locus_tag"] if g_R2 else None,
    )

    # Curated CSV family/exception overrides take precedence over fresh rule-based classification.
    apply_existing_family_overrides(genes, cfg_path)

    # ── Aux locus rows: auxiliary search probes that may come from any genome.
    con_aux = sqlite3.connect(tgt_cfg["db"])
    con_aux.row_factory = sqlite3.Row
    aux_rows_raw = read_aux_rows(cfg_path)
    aux_genes = resolve_aux_rows(con_aux, aux_rows_raw) if aux_rows_raw else []
    con_aux.close()
    if aux_genes:
        genes.extend(aux_genes)
        print(f"[build_config] added {len(aux_genes)} aux locus probe(s)")

    # ── expected bp
    by_role = {g["role"]: g for g in genes if g["role"] not in ("inner", "flank_L2", "flank_R2")}
    flank_bp = (by_role["flank_R"]["start"] - by_role["flank_L"]["end"]
                if has_flanks else None)
    cluster_bp = None
    if "cluster_L" in by_role and "cluster_R" in by_role:
        cluster_bp = by_role["cluster_R"]["end"] - by_role["cluster_L"]["start"]

    # L2↔L1 and R1↔R2 gaps (intergenic between outer and inner flanks)
    L2_L1_bp = (g_L["start"] - g_L2["end"]) if g_L2 else None
    R1_R2_bp = (g_R2["start"] - g_R["end"]) if g_R2 else None

    # ── functional_required resolution
    tag_to_gene = {g["locus_tag"]: g for g in genes}
    func_resolved = []
    for lt in func_req:
        g = tag_to_gene.get(lt)
        if not g:
            print(f"[WARN] functional_required '{lt}' not found in extracted genes.")
            func_resolved.append({"locus_tag": lt, "family": None, "product": "NOT FOUND"})
        else:
            func_resolved.append({
                "locus_tag": lt,
                "family":    g.get("family"),
                "product":   g["product"],
            })

    round_pairs = build_round_pairs(genes)

    # ── protein sequence extraction
    all_tags = {g["locus_tag"] for g in genes}
    seqs = extract_sequences(tgt_cfg["faa"], all_tags)
    missing_seqs = all_tags - set(seqs)
    if missing_seqs:
        print(f"[WARN] {len(missing_seqs)} locus_tag(s) not found in FAA: {missing_seqs}")

    # ── print summary
    print(f"\n  contig     : {contig}")
    if g_L2: print(f"  L2↔L1 bp   : {L2_L1_bp:,} bp")
    print(f"  flank_bp   : {flank_bp:,} bp" if flank_bp is not None else "  flank_bp   : N/A (flank-less)")
    if g_R2: print(f"  R1↔R2 bp   : {R1_R2_bp:,} bp")
    print(f"  cluster_bp : {cluster_bp:,} bp" if cluster_bp else "  cluster_bp : N/A")
    print(f"  n_genes_ref: {len(genes)}")
    print()
    print(f"  {'role':12s} {'locus_tag':22s} {'family':12s} {'exc':7s} product")
    print(f"  {'-'*12} {'-'*22} {'-'*12} {'-'*7} {'-'*40}")
    for g in genes:
        ri = f"  inner[{g['inner_idx']}]" if g["role"] == "inner" else f"  {g['role']}"
        exc_s = "YES" if g.get("exception") else ""
        print(f"  {g['role']:12s} {g['locus_tag']:22s} {str(g.get('family','—')):12s} {exc_s:7s} {g['product'][:60]}")
    print()
    print(f"  Progressive round pairs:")
    for rp in round_pairs:
        print(f"    Round {rp['round']}: {rp['left_tag']}  ↔  {rp['right_tag']}")

    if args.dry_run:
        print("\n[dry-run] No files written.")
        return

    # ── write anchors.faa
    anchors_faa = cfg_path.with_name(cfg_path.stem + "_anchors.faa")
    with open(anchors_faa, "w") as f:
        for g in genes:
            seq = seqs.get(g["locus_tag"], "")
            if seq:
                header = f">{g['role']}|{g['locus_tag']}|{g.get('family','unknown')}|{g['product'][:60]}"
                f.write(f"{header}\n{seq}\n")
    print(f"[build_config] wrote → {anchors_faa}  ({sum(1 for g in genes if g['locus_tag'] in seqs)} sequences)")

    # ── write anchors.csv (EDITABLE)
    anchors_csv = write_anchors_csv(genes, cfg_path)
    print(f"[build_config] wrote → {anchors_csv} (EDIT THIS to change family labels and exceptions)")

    # ── extract cluster nucleotide FASTA (for cluster blastn)
    cluster_fna_path = None
    cluster_ref_bp   = None
    if "cluster_L" in by_role and "cluster_R" in by_role:
        cl_g = by_role["cluster_L"]
        cr_g = by_role["cluster_R"]
        cl_s = int(cl_g["start"])
        cr_e = int(cr_g["end"])
        cseq = extract_cluster_fna(tgt_cfg["fna_dir"], accession, contig, cl_s, cr_e)
        if cseq:
            cluster_fna_path = cfg_path.with_name(cfg_path.stem + "_cluster.fna")
            with open(cluster_fna_path, "w") as f:
                f.write(f">{accession}|{contig}|{cl_s}:{cr_e}|cluster\n")
                for i in range(0, len(cseq), 80):
                    f.write(cseq[i:i+80] + "\n")
            cluster_ref_bp = len(cseq)
            print(f"[build_config] cluster FASTA → {cluster_fna_path}  ({cluster_ref_bp:,} bp)")
        else:
            print(f"[WARN] could not extract cluster FASTA for {accession} "
                  f"(contig={contig} {cl_s}-{cr_e}); check fna_dir and coordinates")
    else:
        print("[build_config] no cluster_L/cluster_R defined; cluster FASTA skipped")

    # Reference genus
    ref_genus = ""
    try:
        con = sqlite3.connect(tgt_cfg["db"])
        row = con.execute(
            "SELECT species FROM genomes WHERE accession=?",
            (accession,),
        ).fetchone()
        con.close()
        if row and row[0]:
            ref_genus = row[0].strip().split()[0]
    except Exception as e:
        print(f"[build_config] WARN: could not look up reference genus: {e}")

    # ── update _auto block
    locus_cfg["_auto"] = {
        "built_at":    datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        "built_from":  {"db": tgt_cfg["db"], "faa": tgt_cfg["faa"]},
        "contig":      contig,
        "reference_genus": ref_genus,
        "_anchors_faa": str(anchors_faa),
        **({"_cluster_fna":   str(cluster_fna_path),
            "cluster_ref_bp": cluster_ref_bp} if cluster_fna_path else {}),
        "expected_bp": {
            **({"L2_L1":   L2_L1_bp}  if L2_L1_bp  is not None else {}),
            **({"flank":   flank_bp}  if flank_bp   is not None else {}),
            **({"R1_R2":   R1_R2_bp}  if R1_R2_bp   is not None else {}),
            **({"cluster": cluster_bp} if cluster_bp is not None else {}),
        },
        "n_genes_ref":    len(genes),
        "anchors":        [
            {k: v for k, v in g.items() if k != "inner_idx"}
            for g in genes
        ],
        "functional_required_resolved": func_resolved,
        "round_pairs": round_pairs,
    }

    save_locus_cfg(cfg_path, locus_cfg)
    print(f"[build_config] wrote → {cfg_path}  (_auto populated)")
    print_locus_preview(genes, accession, locus_cfg.get("locus_id", ""),
                        cluster_ref_bp)


def print_locus_preview(genes, accession, locus_id, cluster_ref_bp=None):
    """Show what the chosen flanks actually bracket, before any BLAST runs."""
    def _role(gene):
        return (gene.get("role") or "").strip()

    left = [g for g in genes if _role(g).startswith("flank_L")]
    right = [g for g in genes if _role(g).startswith("flank_R")]
    inside = [g for g in genes
              if _role(g) in ("cluster_L", "cluster_R", "inner")]

    def _line(gene, index=None):
        tag = gene.get("locus_tag", "?")
        product = (gene.get("product") or "").strip() or "(no product)"
        marker = f"{index:>3}" if index is not None else "   "
        return f"    {marker}  {tag:<18s} {product[:56]}"

    print()
    print(f"[build_config] locus {locus_id!r} resolved from {accession}")
    for gene in sorted(left, key=lambda g: _role(g), reverse=True)[-1:]:
        print(f"    flank {_line(gene).strip()}")
    span = f", {cluster_ref_bp:,} bp" if cluster_ref_bp else ""
    print(f"    {'-' * 12} {len(inside)} genes{span} {'-' * 12}")
    for index, gene in enumerate(inside, start=1):
        print(_line(gene, index))
    for gene in sorted(right, key=lambda g: _role(g))[:1]:
        print(f"    flank {_line(gene).strip()}")
    print()
    print("    Is this the locus you meant? If not, rerun make-config with")
    print("    different flanks -- everything after this step is the slow part.")


if __name__ == "__main__":
    main()
