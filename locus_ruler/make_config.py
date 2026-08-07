#!/usr/bin/env python3
"""Create a minimal LocusRuler locus config from the cluster's two ends."""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
from pathlib import Path

_PKG_DIR = Path(__file__).resolve().parent
if str(_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(_PKG_DIR))

from config_utils import get_target_cfg, load_settings, resolve_target_name


def _load_table(path: Path) -> dict[str, dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"ERROR: reference table is empty: {path}")
    if "locus_tag" not in rows[0]:
        sys.exit(f"ERROR: reference table lacks a locus_tag column: {path}")
    return {row["locus_tag"]: row for row in rows}


def _normalize_flank_order(table: dict[str, dict], flank_l: str, flank_r: str) -> tuple[str, str]:
    if flank_l not in table:
        sys.exit(f"ERROR: --flank-l '{flank_l}' not found in --table")
    if flank_r not in table:
        sys.exit(f"ERROR: --flank-r '{flank_r}' not found in --table")

    left = table[flank_l]
    right = table[flank_r]
    if left.get("contig") and right.get("contig") and left["contig"] != right["contig"]:
        sys.exit("ERROR: --flank-l and --flank-r are on different contigs")

    try:
        if int(left.get("start", 0)) > int(right.get("start", 0)):
            return flank_r, flank_l
    except ValueError:
        pass
    return flank_l, flank_r


def _flanks_from_ends(db_path: str, accession: str,
                       first_tag: str, last_tag: str) -> tuple[str, str]:
    """The genes just outside the cluster, given the cluster's own ends.

    The config brackets a locus by its flanks, but naming the first and last
    gene *of* the cluster is the way people actually think about it -- and
    the way the wizard asks. Look the neighbors up rather than making the
    caller count genes outward by hand.
    """
    con = sqlite3.connect(str(db_path))
    try:
        placed = {}
        for tag in (first_tag, last_tag):
            row = con.execute(
                "SELECT contig, start, end FROM proteins "
                "WHERE genome_acc=? AND locus_tag=?", (accession, tag)).fetchone()
            if row is None:
                sys.exit(f"ERROR: locus_tag '{tag}' not found in {accession}")
            placed[tag] = row
        (contig_a, *_), (contig_b, *_) = placed[first_tag], placed[last_tag]
        if contig_a != contig_b:
            sys.exit(f"ERROR: --first and --last are on different contigs "
                     f"({contig_a} and {contig_b})")
        if first_tag == last_tag:
            sys.exit("ERROR: --first and --last are the same gene; a locus "
                     "needs at least two")

        lo = min(placed[first_tag][1], placed[last_tag][1])
        hi = max(placed[first_tag][2], placed[last_tag][2])
        left = con.execute(
            "SELECT locus_tag FROM proteins WHERE genome_acc=? AND contig=? "
            "AND end < ? ORDER BY start DESC LIMIT 1",
            (accession, contig_a, lo)).fetchone()
        right = con.execute(
            "SELECT locus_tag FROM proteins WHERE genome_acc=? AND contig=? "
            "AND start > ? ORDER BY start LIMIT 1",
            (accession, contig_a, hi)).fetchone()
    finally:
        con.close()

    if left is None or right is None:
        missing = "before --first" if left is None else "after --last"
        sys.exit(f"ERROR: no gene {missing} on {contig_a} to use as a flank; "
                 f"the cluster runs to the end of the contig")
    return left[0], right[0]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create a minimal LocusRuler JSON config for one locus."
    )
    ap.add_argument("--locus-id", required=True, help="Output locus_id")
    ap.add_argument("--accession", required=True, help="Reference genome accession")
    ap.add_argument("--settings", default="settings.toml",
                    help="Path to settings.toml (default: settings.toml in "
                         "the current directory)")
    ap.add_argument("--target", default=None,
                    help="Target name in settings.toml. Default: the only "
                         "declared target; with several and no target given, "
                         "you are asked")
    ap.add_argument("--first", dest="first_tag", default=None,
                    help="First gene OF the cluster; its flank is looked up. "
                         "Use with --last instead of --flank-l/--flank-r")
    ap.add_argument("--last", dest="last_tag", default=None,
                    help="Last gene OF the cluster; its flank is looked up")
    ap.add_argument("--flank-l", default=None,
                    help="Left flank locus_tag -- the gene just OUTSIDE the "
                         "cluster. Use with --flank-r instead of "
                         "--first/--last")
    ap.add_argument("--flank-r", default=None, help="Right flank locus_tag")
    ap.add_argument("--out", default=None,
                    help="Output JSON path (default: <locus-id>.json)")
    ap.add_argument(
        "--table",
        default=None,
        help="Optional reference table CSV from locus-ruler-genes. Used to "
             "validate and order --flank-l/--flank-r.",
    )
    args = ap.parse_args()

    ends = (args.first_tag, args.last_tag)
    flanks = (args.flank_l, args.flank_r)
    if all(ends) and any(flanks):
        sys.exit("ERROR: give either --first/--last or --flank-l/--flank-r, "
                 "not both")
    if not all(ends) and not all(flanks):
        sys.exit("ERROR: name the locus with --first and --last (the "
                 "cluster's own first and last gene), or with --flank-l and "
                 "--flank-r (the genes just outside it)")

    settings_path = Path(args.settings)
    if not settings_path.exists():
        sys.exit(f"ERROR: settings file not found: {settings_path}\n"
                 f"       Pass --settings <path>, or run from the directory "
                 f"holding settings.toml.")
    settings = load_settings(settings_path.resolve())
    target_name = resolve_target_name(settings, args.target)

    if all(ends):
        target_cfg = get_target_cfg(settings, target_name)
        flank_l, flank_r = _flanks_from_ends(
            target_cfg["db"], args.accession, args.first_tag, args.last_tag)
    else:
        flank_l, flank_r = flanks
        if args.table:
            flank_l, flank_r = _normalize_flank_order(
                _load_table(Path(args.table)), flank_l, flank_r)

    cfg = {
        "_description": (
            "Minimal LocusRuler config. The target selects the genome set; "
            "flank_L and flank_R bracket the reference locus inside accession."
        ),
        "locus_id": args.locus_id,
        "reference": {
            "accession": args.accession,
            "target": target_name,
            "flank_L": flank_l,
            "flank_R": flank_r,
        },
    }

    out_path = Path(args.out) if args.out else Path(f"{args.locus_id}.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    print(f"[make_config] wrote: {out_path}")
    if all(ends):
        print(f"[make_config] flanks from --first/--last: "
              f"{flank_l} .. {flank_r}")
    elif args.table:
        print(f"[make_config] validated flanks against: {args.table}")


if __name__ == "__main__":
    main()
