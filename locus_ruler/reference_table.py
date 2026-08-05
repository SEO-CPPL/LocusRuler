#!/usr/bin/env python3
"""Export a reference genome gene table to help choose LocusRuler flanks."""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
from pathlib import Path
from urllib.parse import unquote

_PKG_DIR = Path(__file__).resolve().parent
if str(_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(_PKG_DIR))

from config_utils import get_target_cfg, load_settings
from gff import _find_gff, _open_gff, _parse_attrs


def _genome_meta(db_path: str, accession: str) -> dict:
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    row = con.execute(
        "SELECT * FROM genomes WHERE accession = ?",
        (accession,),
    ).fetchone()
    con.close()
    if not row:
        sys.exit(f"ERROR: accession '{accession}' not found in {db_path}")
    return dict(row)


def _parse_reference_gff(gff_path: Path) -> list[dict]:
    rows: list[dict] = []
    with _open_gff(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig, _, ftype, start_s, end_s, _, strand, _, attr_s = parts[:9]
            if ftype not in {"CDS", "gene", "pseudogene", "rRNA", "tRNA"}:
                continue
            attrs = _parse_attrs(attr_s)
            locus_tag = attrs.get("locus_tag") or attrs.get("Name") or attrs.get("ID") or ""
            if not locus_tag:
                continue
            is_pseudo = (
                ftype == "pseudogene"
                or attrs.get("pseudo", "").lower() == "true"
                or attrs.get("pseudogene", "").lower() == "true"
            )
            if ftype == "gene" and not is_pseudo:
                continue
            product = unquote(attrs.get("product", attrs.get("Name", "")))
            rows.append({
                "contig": contig,
                "start": int(start_s),
                "end": int(end_s),
                "strand": strand,
                "locus_tag": locus_tag,
                "gene": unquote(attrs.get("gene", "")),
                "feature_type": "pseudogene" if is_pseudo else ftype,
                "product": product,
                "protein_id": attrs.get("protein_id", ""),
                "is_pseudo": "Y" if is_pseudo else "",
            })

    rows.sort(key=lambda r: (r["contig"], r["start"], r["end"], r["locus_tag"]))

    # Prefer CDS/pseudogene rows over generic gene rows when locus_tag repeats.
    priority = {"CDS": 0, "pseudogene": 1, "rRNA": 2, "tRNA": 2, "gene": 9}
    deduped: list[dict] = []
    seen: dict[tuple[str, str], int] = {}
    for row in rows:
        key = (row["contig"], row["locus_tag"])
        if key not in seen:
            seen[key] = len(deduped)
            deduped.append(row)
            continue
        old = deduped[seen[key]]
        if priority.get(row["feature_type"], 9) < priority.get(old["feature_type"], 9):
            deduped[seen[key]] = row

    for contig in sorted({r["contig"] for r in deduped}):
        ctg_rows = [r for r in deduped if r["contig"] == contig]
        for idx, row in enumerate(ctg_rows, start=1):
            row["contig_order"] = idx
            row["left_locus_tag"] = ctg_rows[idx - 2]["locus_tag"] if idx > 1 else ""
            row["right_locus_tag"] = ctg_rows[idx]["locus_tag"] if idx < len(ctg_rows) else ""
    return deduped


def _select_rows(
    rows: list[dict],
    around: str | None,
    product_regex: str | None,
    radius: int,
) -> list[dict]:
    if not around and not product_regex:
        return rows

    keep: set[int] = set()
    by_contig: dict[str, list[tuple[int, dict]]] = {}
    for i, row in enumerate(rows):
        by_contig.setdefault(row["contig"], []).append((i, row))

    def add_window(contig: str, contig_order: int) -> None:
        ctg = by_contig[contig]
        lo = max(1, contig_order - radius)
        hi = contig_order + radius
        for global_idx, row in ctg:
            if lo <= int(row["contig_order"]) <= hi:
                keep.add(global_idx)

    if around:
        matches = [r for r in rows if r["locus_tag"] == around]
        if not matches:
            sys.exit(f"ERROR: --around locus_tag '{around}' not found")
        for row in matches:
            add_window(row["contig"], int(row["contig_order"]))

    if product_regex:
        rx = re.compile(product_regex, re.IGNORECASE)
        matches = [
            r for r in rows
            if rx.search(r.get("product", "")) or rx.search(r.get("gene", ""))
        ]
        if not matches:
            sys.exit(f"ERROR: --product pattern '{product_regex}' matched no genes")
        for row in matches:
            add_window(row["contig"], int(row["contig_order"]))

    return [row for i, row in enumerate(rows) if i in keep]


def _write_csv(rows: list[dict], out_path: Path, accession: str, target: str) -> None:
    fieldnames = [
        "accession",
        "target",
        "contig",
        "contig_order",
        "locus_tag",
        "left_locus_tag",
        "right_locus_tag",
        "gene",
        "feature_type",
        "is_pseudo",
        "strand",
        "start",
        "end",
        "protein_id",
        "product",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            payload = {k: row.get(k, "") for k in fieldnames}
            payload["accession"] = accession
            payload["target"] = target
            writer.writerow(payload)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Export a reference genome gene table to help choose flank_L/flank_R."
    )
    ap.add_argument("--settings", required=True, help="Path to settings.toml")
    ap.add_argument("--target", required=True, help="Target name in settings.toml")
    ap.add_argument("--accession", required=True, help="Reference genome accession")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--around", default=None, help="Only export genes around this locus_tag")
    ap.add_argument("--product", default=None, help="Only export windows around product/gene regex matches")
    ap.add_argument("--radius", type=int, default=20, help="Genes on each side for --around/--product")
    args = ap.parse_args()

    settings = load_settings(Path(args.settings).resolve())
    target_cfg = get_target_cfg(settings, args.target)
    meta = _genome_meta(target_cfg["db"], args.accession)
    gff_path = _find_gff(target_cfg["gff_dir"], args.accession)
    if not gff_path:
        sys.exit(f"ERROR: no GFF found for {args.accession} in {target_cfg['gff_dir']}")

    rows = _parse_reference_gff(gff_path)
    selected = _select_rows(rows, args.around, args.product, args.radius)
    _write_csv(selected, Path(args.out), args.accession, args.target)
    print(f"[reference_table] accession: {args.accession}")
    print(f"[reference_table] organism : {meta.get('organism_name') or meta.get('species') or ''}")
    print(f"[reference_table] gff      : {gff_path}")
    print(f"[reference_table] wrote    : {args.out} ({len(selected)} rows)")


if __name__ == "__main__":
    main()
