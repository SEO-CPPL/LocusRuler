#!/usr/bin/env python3
"""Create a minimal LocusRuler locus config from selected flank tags."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


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


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create a minimal LocusRuler JSON config from selected flank tags."
    )
    ap.add_argument("--locus-id", required=True, help="Output locus_id")
    ap.add_argument("--accession", required=True, help="Reference genome accession")
    ap.add_argument("--target", required=True, help="Target name in settings.toml")
    ap.add_argument("--flank-l", required=True, help="Left flank locus_tag")
    ap.add_argument("--flank-r", required=True, help="Right flank locus_tag")
    ap.add_argument("--out", required=True, help="Output JSON path")
    ap.add_argument(
        "--table",
        default=None,
        help="Optional reference table CSV from reference_table.py. Used to validate and order flanks.",
    )
    args = ap.parse_args()

    flank_l = args.flank_l
    flank_r = args.flank_r
    if args.table:
        flank_l, flank_r = _normalize_flank_order(_load_table(Path(args.table)), flank_l, flank_r)

    cfg = {
        "_description": (
            "Minimal LocusRuler config. The target selects the genome set; "
            "flank_L and flank_R bracket the reference locus inside accession."
        ),
        "locus_id": args.locus_id,
        "reference": {
            "accession": args.accession,
            "target": args.target,
            "flank_L": flank_l,
            "flank_R": flank_r,
        },
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    print(f"[make_config] wrote: {out_path}")
    if args.table:
        print(f"[make_config] validated flanks against: {args.table}")


if __name__ == "__main__":
    main()
