#!/usr/bin/env python3
"""Pipeline entry point."""

import argparse
import json
import sys
from pathlib import Path

# Support both `python locus_ruler/run.py` and `python -m locus_ruler.run`, keeping
# the script-style intra-package imports.
_PKG_DIR = Path(__file__).resolve().parent
if str(_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(_PKG_DIR))

from config_utils import (load_settings, load_locus_cfg, save_locus_cfg,
                          get_target_cfg, resolve_target_name)


# ── Step helpers ──────────────────────────────────────────────────────
def step_build_config(
    cfg_path: Path,
    settings_path: Path,
    dry_run: bool,
    force_rebuild: bool = False,
) -> dict:
    """Populate the _auto block of the locus config."""
    from build_config import (
        load_rules, fetch_gene_pos, fetch_genes_between, fetch_adjacent_gene,
        extract_sequences, extract_cluster_fna,
        classify_product, assign_roles, build_round_pairs,
        write_anchors_csv, apply_existing_family_overrides,
        read_aux_rows, resolve_aux_rows,
    )
    import sqlite3
    from datetime import datetime

    settings  = load_settings(settings_path)
    locus_cfg = load_locus_cfg(cfg_path)

    if not force_rebuild and locus_cfg.get("_auto") and locus_cfg["_auto"].get("built_at"):
        print(f"[run] Step 1: _auto already populated "
              f"(built_at={locus_cfg['_auto']['built_at']}); skipping build_config "
              f"(use --rebuild to force)")
        return locus_cfg

    print("[run] Step 1: running build_config …")
    if dry_run:
        print("  [dry-run] would rebuild _auto block")
        return locus_cfg

    ref       = locus_cfg["reference"]
    accession = ref["accession"]
    target    = ref["target"]
    flank_L   = ref.get("flank_L")
    flank_R   = ref.get("flank_R")
    has_flanks = bool(flank_L and flank_R)
    cl_direct  = ref.get("cluster_L")
    cr_direct  = ref.get("cluster_R")

    if not has_flanks and (not cl_direct or not cr_direct):
        sys.exit("ERROR: flank-less mode (Tier 1) requires reference.cluster_L "
                 "and reference.cluster_R.")

    tgt_cfg    = get_target_cfg(settings, target)
    rules_file = settings.get("classification", {}).get("rules_file", "")
    rules      = load_rules(rules_file,
                            locus_id   = locus_cfg.get("locus_id"),
                            config_dir = cfg_path.parent)

    con = sqlite3.connect(tgt_cfg["db"])
    con.row_factory = sqlite3.Row

    if has_flanks:
        g_L = fetch_gene_pos(con, accession, flank_L)
        g_R = fetch_gene_pos(con, accession, flank_R)
        if not g_L:
            sys.exit(f"ERROR: flank_L '{flank_L}' not found in DB for {accession}")
        if not g_R:
            sys.exit(f"ERROR: flank_R '{flank_R}' not found in DB for {accession}")
        if g_L["contig"] != g_R["contig"]:
            sys.exit(f"ERROR: flank_L ({g_L['contig']}) and flank_R ({g_R['contig']}) "
                     "are on different contigs")
        contig = g_L["contig"]
        cmin   = min(g_L["start"], g_R["start"])
        cmax   = max(g_L["end"],   g_R["end"])
    else:
        g_CL = fetch_gene_pos(con, accession, cl_direct)
        g_CR = fetch_gene_pos(con, accession, cr_direct)
        if not g_CL:
            sys.exit(f"ERROR: cluster_L '{cl_direct}' not found in DB for {accession}")
        if not g_CR:
            sys.exit(f"ERROR: cluster_R '{cr_direct}' not found in DB for {accession}")
        if g_CL["contig"] != g_CR["contig"]:
            sys.exit(f"ERROR: cluster_L ({g_CL['contig']}) and cluster_R "
                     f"({g_CR['contig']}) are on different contigs")
        contig = g_CL["contig"]
        cmin   = min(g_CL["start"], g_CR["start"])
        cmax   = max(g_CL["end"],   g_CR["end"])

    genes_raw = fetch_genes_between(con, accession, contig, cmin, cmax)

    # L2 / R2 outer flanks (Tier 2+ only)
    g_L2 = g_R2 = None
    if has_flanks:
        g_L2 = fetch_adjacent_gene(con, accession, contig, g_L["start"], "left")
        g_R2 = fetch_adjacent_gene(con, accession, contig, g_R["end"],   "right")

    con.close()

    for g in genes_raw:
        g["family"] = classify_product(g["product"], rules)
    for g in [g_L2, g_R2]:
        if g:
            g["family"] = classify_product(g["product"], rules)

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

    # Curated CSV family overrides take precedence over fresh rule-based classification.
    apply_existing_family_overrides(genes, cfg_path)

    # Aux rows live only in the CSV, so read them back before write_anchors_csv drops them
    con_aux = sqlite3.connect(tgt_cfg["db"])
    con_aux.row_factory = sqlite3.Row
    aux_rows_raw = read_aux_rows(cfg_path)
    aux_genes = resolve_aux_rows(con_aux, aux_rows_raw) if aux_rows_raw else []
    con_aux.close()
    if aux_genes:
        genes.extend(aux_genes)
        print(f"[run]   added {len(aux_genes)} aux locus probe(s)")

    by_role = {g["role"]: g for g in genes
               if g["role"] not in ("inner", "flank_L2", "flank_R2")}
    flank_bp   = (by_role["flank_R"]["start"] - by_role["flank_L"]["end"]
                  if has_flanks else None)
    cluster_bp = None
    if "cluster_L" in by_role and "cluster_R" in by_role:
        cluster_bp = by_role["cluster_R"]["end"] - by_role["cluster_L"]["start"]

    # L2↔L1 and R1↔R2 intergenic gaps
    g_L_info = next((g for g in genes if g["role"] == "flank_L"),  None)
    g_R_info = next((g for g in genes if g["role"] == "flank_R"),  None)
    L2_L1_bp = (g_L_info["start"] - g_L2["end"]) if g_L2 and g_L_info else None
    R1_R2_bp = (g_R2["start"] - g_R_info["end"]) if g_R2 and g_R_info else None

    all_tags = {g["locus_tag"] for g in genes}
    seqs     = extract_sequences(tgt_cfg["faa"], all_tags)

    anchors_faa = cfg_path.with_name(cfg_path.stem + "_anchors.faa")
    with open(anchors_faa, "w") as f:
        for g in genes:
            seq = seqs.get(g["locus_tag"], "")
            if seq:
                header = (f">{g['role']}|{g['locus_tag']}|"
                          f"{g.get('family','unknown')}|{g['product'][:60]}")
                f.write(f"{header}\n{seq}\n")
    n_seqs = sum(1 for g in genes if g["locus_tag"] in seqs)
    print(f"[run]   anchors.faa → {anchors_faa}  ({n_seqs} sequences)")

    # Editable per-locus family CSV (preserves manual edits — see build_config.write_anchors_csv).
    anchors_csv = write_anchors_csv(genes, cfg_path)
    print(f"[run]   anchors.csv → {anchors_csv}")

    # Cluster nucleotide FASTA (for blastn step)
    cluster_fna_path = None
    cluster_ref_bp   = None
    if "cluster_L" in by_role and "cluster_R" in by_role:
        cl_g = by_role["cluster_L"]
        cr_g = by_role["cluster_R"]
        cseq = extract_cluster_fna(
            tgt_cfg["fna_dir"], accession, contig,
            int(cl_g["start"]), int(cr_g["end"]),
        )
        if cseq:
            cluster_fna_path = cfg_path.with_name(cfg_path.stem + "_cluster.fna")
            with open(cluster_fna_path, "w") as f:
                f.write(f">{accession}|{contig}|{cl_g['start']}:{cr_g['end']}|cluster\n")
                for i in range(0, len(cseq), 80):
                    f.write(cseq[i:i+80] + "\n")
            cluster_ref_bp = len(cseq)
            print(f"[run]   cluster FASTA → {cluster_fna_path}  ({cluster_ref_bp:,} bp)")
        else:
            print(f"[run]   WARN: could not extract cluster FASTA for {accession} "
                  f"on {contig} {cl_g['start']}–{cr_g['end']}; check fna_dir")

    round_pairs = build_round_pairs(genes)
    func_req    = locus_cfg.get("functional_required", [])
    tag_to_gene = {g["locus_tag"]: g for g in genes}
    func_resolved = []
    for lt in func_req:
        g = tag_to_gene.get(lt)
        if g:
            func_resolved.append({"locus_tag": lt,
                                   "family":    g.get("family"),
                                   "product":   g["product"]})
        else:
            func_resolved.append({"locus_tag": lt, "family": None,
                                   "product": "NOT FOUND"})

    locus_cfg["_auto"] = {
        "built_at":    datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        "built_from":  {"db": tgt_cfg["db"], "faa": tgt_cfg["faa"]},
        "contig":      contig,
        "_anchors_faa": str(anchors_faa),
        **({"_cluster_fna":   str(cluster_fna_path),
            "cluster_ref_bp": cluster_ref_bp} if cluster_fna_path else {}),
        "expected_bp": {
            **({"L2_L1":   L2_L1_bp}  if L2_L1_bp  is not None else {}),
            **({"flank":   flank_bp}  if flank_bp   is not None else {}),
            **({"R1_R2":   R1_R2_bp}  if R1_R2_bp   is not None else {}),
            **({"cluster": cluster_bp} if cluster_bp is not None else {}),
        },
        "n_genes_ref":  len(genes),
        "anchors":      [{k: v for k, v in g.items() if k != "inner_idx"}
                         for g in genes],
        "functional_required_resolved": func_resolved,
        "round_pairs": round_pairs,
    }
    save_locus_cfg(cfg_path, locus_cfg)
    print(f"[run]   config updated → {cfg_path}")
    if g_L2: print(f"[run]   L2↔L1 = {L2_L1_bp:,} bp  ({g_L2['locus_tag']} → {flank_L})")
    if g_R2: print(f"[run]   R1↔R2 = {R1_R2_bp:,} bp  ({flank_R} → {g_R2['locus_tag']})")
    return locus_cfg


# ── Main ──────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="LocusRuler — genomic locus status surveyor",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--config",   required=True,
                    help="Path to locus JSON config (e.g. configs/example_cluster.json)")
    ap.add_argument("--settings", default="settings.toml",
                    help="Path to settings.toml (default: settings.toml in "
                         "the current directory)")
    ap.add_argument("--target",
                    help="Target name (overrides reference.target in config). "
                         "Default: the config's reference target, or the "
                         "only declared target; with several and no target "
                         "given, you are asked")
    ap.add_argument("--cpu",     type=int, default=4,
                    help="Parallel workers for BLAST jobs (default: 4)")
    ap.add_argument("--from",    dest="from_step", type=int, default=1,
                    metavar="STEP",
                    help="Resume from step N (1=build_config … 6=cassette structure)")
    ap.add_argument("--to",      dest="to_step", type=int, default=6,
                    metavar="STEP",
                    help="Stop after step N. --to 1 writes the anchor table and "
                         "stops, so it can be edited before any BLAST runs.")
    ap.add_argument("--genome",
                    help="Run only for this accession (for debugging)")
    ap.add_argument("--force",   action="store_true",
                    help="Overwrite cached BLAST results")
    ap.add_argument("--no-split-search", action="store_true",
                    help="Skip floating tblastn even if SPLIT is suspected")
    ap.add_argument("--rebuild", action="store_true",
                    help="Force rebuild of _auto block even if built_at is already set")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print steps without executing")
    ap.add_argument("--output-dir",
                    help="Override output directory (default: settings output_root)")
    # On by default: a plain run should produce everything the output guide describes.
    ap.add_argument("--no-heatmap", dest="heatmap", action="store_false",
                    help="Skip cluster_heatmap.png (render it later with "
                         "locus-ruler-heatmap)")
    ap.add_argument("--no-report", dest="report", action="store_false",
                    help="Skip locus_report.xlsx (build it later with "
                         "locus-ruler-report)")
    ap.set_defaults(heatmap=True, report=True)
    return ap.parse_args()


def main():
    args        = parse_args()
    cfg_path    = Path(args.config).resolve()
    settings_p  = Path(args.settings)
    if not settings_p.exists():
        sys.exit(f"ERROR: settings file not found: {settings_p}\n"
                  f"       Pass --settings <path>, or run from the "
                  f"directory holding settings.toml.")
    settings_p  = settings_p.resolve()
    settings    = load_settings(settings_p)

    # Resolve paths
    work_dir   = Path(settings["paths"]["work_dir"])
    output_root = Path(settings["paths"].get("output_root",
                                              str(work_dir.parent / "output")))
    if args.output_dir:
        output_root = Path(args.output_dir)

    # ── Step 1: build_config 
    if args.from_step <= 1:
        locus_cfg = step_build_config(cfg_path, settings_p, args.dry_run,
                                      force_rebuild=args.rebuild)
    else:
        locus_cfg = load_locus_cfg(cfg_path)

    locus_id   = locus_cfg["locus_id"]
    target_name = resolve_target_name(settings, args.target, locus_cfg)
    tgt_cfg    = get_target_cfg(settings, target_name)

    locus_work  = work_dir / target_name / locus_id
    locus_out   = output_root / target_name / locus_id
    locus_work.mkdir(parents=True, exist_ok=True)
    locus_out.mkdir(parents=True, exist_ok=True)

    print(f"[run] locus_id   : {locus_id}")
    print(f"[run] target     : {target_name}")
    print(f"[run] work_dir   : {locus_work}")
    print(f"[run] output_dir : {locus_out}")

    if args.dry_run:
        print("[run] --dry-run: stopping after Step 1 (no BLAST executed)")
        return

    if args.to_step <= 1:
        # Stopping here exposes the anchor table, which every later step reads.
        anchors_csv = cfg_path.with_name(cfg_path.stem + "_anchors.csv")
        print(f"\n[run] Stopped after Step 1 (--to 1). No BLAST was run.")
        if anchors_csv.exists():
            # Echo the invocation as typed (sys.argv[0]), not a bare "locus-ruler".
            print(f"[run] Anchor table: {anchors_csv}")
            print("[run] Edit family / status_role / exception / lenient, then:")
            print(f"[run]   {sys.argv[0]} --config {cfg_path} "
                  f"--settings {settings_p} --target {target_name} --from 2")
        return

    # ── Step 2: make_dbs 
    from blast import make_dbs, run_anchor_blast, run_cluster_blastn
    from ruler import run_ruler, load_ruler_results
    from content import run_content

    tools = settings.get("tools", {})
    if args.from_step <= 2:
        print("[run] Step 2: building BLAST databases …")
        db_map = make_dbs(
            tgt_cfg, work_dir, cpu=args.cpu,
            makeblastdb=tools.get("makeblastdb", "makeblastdb"),
        )
    else:
        # Reconstruct db_map from existing files, restricted to TARGET's accessions.
        import sqlite3
        tgt_db_path = tgt_cfg.get("db")
        target_accs: set[str] = set()
        if tgt_db_path and Path(tgt_db_path).exists():
            con = sqlite3.connect(tgt_db_path)
            target_accs = {r[0] for r in
                           con.execute("SELECT accession FROM genomes").fetchall()}
            con.close()

        db_base = work_dir / "blastdbs"
        db_map  = {}
        skipped = 0
        for d in db_base.iterdir():
            if not d.is_dir() or not (d / f"{d.name}.nsq").exists():
                continue
            if target_accs and d.name not in target_accs:
                skipped += 1
                continue
            db_map[d.name] = str(d / d.name)
        msg = f"[run] Step 2: loaded {len(db_map)} existing BLAST DBs for target '{target_name}'"
        if skipped:
            msg += f"  (skipped {skipped} DBs belonging to other targets)"
        print(msg)

    # Filter to a single genome if --genome specified
    if args.genome:
        if args.genome not in db_map:
            sys.exit(f"ERROR: accession '{args.genome}' has no BLAST DB in {work_dir/'blastdbs'}")
        db_map = {args.genome: db_map[args.genome]}
        print(f"[run] Filtering to single genome: {args.genome}")

    if args.to_step <= 2:
        print("[run] Stopped after Step 2 (--to 2). BLAST databases are built.")
        return

    # ── Step 3: anchor BLAST 
    if args.from_step <= 3:
        print("[run] Step 3: running anchor tblastn …")
        blast_hits = run_anchor_blast(
            locus_cfg, settings, db_map, locus_work.parent,
            cpu=args.cpu, force=args.force,
        )
    else:
        # Load cached TSVs on demand (holding every genome's is costly)
        from blast import anchor_hits_on_disk
        anchors = locus_cfg["_auto"]["anchors"]
        tags    = [a["locus_tag"] for a in anchors]
        blast_hits = anchor_hits_on_disk(locus_work.parent, locus_id,
                                         list(db_map), tags)
        print(f"[run] Step 3: using cached blast hits for {len(db_map)} genomes "
              f"(read on demand)")

    # ── Step 3b: cluster blastn 
    cluster_hits: dict = {}
    if args.from_step <= 3:
        cluster_fna = locus_cfg.get("_auto", {}).get("_cluster_fna")
        if cluster_fna:
            print("[run] Step 3b: running cluster blastn …")
            cluster_hits = run_cluster_blastn(
                locus_cfg, settings, db_map, work_dir,
                cpu=args.cpu, force=args.force,
                locus_dir=locus_work,
            )
        else:
            print("[run] Step 3b: no _cluster_fna in config; "
                  "skipping cluster blastn (run build_config first)")
    else:
        # Load cached cluster blastn TSVs on demand
        from blast import cluster_hits_on_disk
        cb_dir = locus_work / "cluster_blast"
        if cb_dir.exists():
            cluster_hits = cluster_hits_on_disk(locus_work, list(db_map))
            print(f"[run] Step 3b: using cached cluster blastn hits "
                  f"for {len(cluster_hits)} genomes (read on demand)")
        else:
            print("[run] Step 3b: no cached cluster blastn found; "
                  "ruler will use gene-level fallback")

    if args.to_step <= 3:
        print("[run] Stopped after Step 3 (--to 3). BLAST hits are cached.")
        return

    # ── Step 4: ruler 
    if args.from_step <= 4:
        print("[run] Step 4: running progressive ruler …")
        ruler_results = run_ruler(
            locus_cfg, settings, db_map, blast_hits, work_dir,
            cluster_hits=cluster_hits if cluster_hits is not None else None,
            cpu=args.cpu,
            force=args.force,
            no_split_search=args.no_split_search,
            target_name=target_name,
            locus_dir=locus_work,
        )
    else:
        print("[run] Step 4: loading cached ruler results …")
        ruler_results = load_ruler_results(locus_work)
        # Restrict to requested genomes
        if args.genome:
            ruler_results = {args.genome: ruler_results[args.genome]}

    if args.to_step <= 4:
        print("[run] Stopped after Step 4 (--to 4). Ruler results are cached.")
        return

    # ── Step 5: content 
    if args.from_step <= 5:
        print("[run] Step 5: extracting GFF content …")
        ruler_results = run_content(
            ruler_results, locus_cfg, settings, blast_hits,
            work_dir=work_dir, output_dir=locus_out,
            target_name=target_name,
            config_dir=cfg_path.parent,
            cluster_hits=cluster_hits,
        )

        # Refresh summary CSV with content columns filled in (genome_meta for species/strain)
        from ruler import write_summary, _load_ruler_genome_meta
        _summary_meta = _load_ruler_genome_meta(
            tgt_cfg.get("db", "") if tgt_cfg else "",
            list(ruler_results.keys()),
        )
        if args.genome:
            # Single-genome debug run: write to the work dir only, leaving the output dir
        # and the ruler_results.json cache untouched.
            write_summary(ruler_results, locus_cfg,
                          locus_work / f"genome_summary_{args.genome}.csv",
                          genome_meta=_summary_meta)
            print(f"[run] --genome mode: debug summary → "
                  f"{locus_work}/genome_summary_{args.genome}.csv")
            print(f"[run] --genome mode: ruler_results.json and output dir NOT updated "
                  f"(would corrupt full-run cache)")
        else:
            write_summary(ruler_results, locus_cfg, locus_work / "genome_summary.csv",
                          genome_meta=_summary_meta)
            # Also write to output dir (run_content already did this, but refresh with genome_meta)
            from writers import tables_dir

            write_summary(ruler_results, locus_cfg,
                          tables_dir(locus_out) / "genome_summary.csv",
                          genome_meta=_summary_meta)
            # Refresh JSON to include _flank_validation and other content fields
            json_path = locus_work / "ruler_results.json"
            with open(json_path, "w") as f:
                json.dump(ruler_results, f, indent=2, default=str)
    else:
        print("[run] Step 5: skipped (--from > 5)")

    if args.to_step <= 5:
        print("[run] Stopped after Step 5 (--to 5). Tables are written.")
        return

    # ── Step 6: generic cassette structure discovery
    if args.genome:
        print("[run] Step 6: skipped in --genome mode "
              "(complete canonical outputs are intentionally unchanged)")
    else:
        print("[run] Step 6: discovering cassette structures …")
        from cassette_structure import run_cassette_structure

        run_cassette_structure(
            cfg_path,
            settings_p,
            target_name=target_name,
            output_root_override=output_root,
        )

    # ── Optional presentation layer 
    if args.heatmap and not args.genome:
        print("[run] Rendering cluster_heatmap.png …")
        try:
            from heatmap import render as render_heatmap

            render_heatmap(
                ruler_results, locus_cfg, _summary_meta,
                locus_out / "cluster_heatmap.png",
            )
        except Exception as exc:
            print(f"[run] WARNING: heatmap rendering failed ({exc}); "
                  f"outputs are otherwise complete. Retry with locus-ruler-heatmap.")

    if args.report and not args.genome:
        print("[run] Building locus_report.xlsx …")
        try:
            from report import build_report

            build_report(locus_out, locus_out / "locus_report.xlsx", locus_id)
        except Exception as exc:
            print(f"[run] WARNING: report build failed ({exc}); "
                  f"outputs are otherwise complete. Retry with locus-ruler-report.")

    # ── Final status summary 
    from collections import Counter
    counts = Counter(r["status"] for r in ruler_results.values())
    print("\n[run] ── Final status counts ─────────────────────────────────")
    for status, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"        {status:20s} {n:4d}")
    print(f"        {'TOTAL':20s} {sum(counts.values()):4d}")
    print(f"\n[run] Output → {locus_out}")


if __name__ == "__main__":
    main()
