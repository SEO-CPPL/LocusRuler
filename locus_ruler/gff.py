#!/usr/bin/env python3
"""GFF3 file I/O and feature parsing."""

import gzip
from pathlib import Path
from typing import Optional


def _open_gff(path: Path):
    """Open a plain or gzip-compressed GFF3 file for text reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _parse_attrs(attr_str: str) -> dict[str, str]:
    """Parse GFF3 attribute column into a dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


def _find_gff(gff_dir: str, accession: str) -> Optional[Path]:
    """Locate the GFF3 file for an accession (tries several naming conventions)."""
    d = Path(gff_dir)
    for suffix in (
        f"{accession}_genomic.gff3.gz",
        f"{accession}_genomic.gff3",
        f"{accession}_genomic.gff.gz",
        f"{accession}_genomic.gff",
        f"{accession}.gff3.gz",
        f"{accession}.gff3",
        f"{accession}.gff.gz",
        f"{accession}.gff",
    ):
        p = d / suffix
        if p.exists():
            return p
    for p in d.iterdir():
        if p.name.startswith(accession) and "gff" in p.name.lower():
            return p
    return None


def parse_gff_region(
    gff_path: Path,
    contig: str,
    region_start: int,
    region_end: int,
    feature_types: tuple[str, ...] = ("CDS", "gene", "pseudogene"),
) -> list[dict]:
    """Extract GFF features on *contig* within [region_start, region_end]."""
    genes: list[dict] = []
    try:
        with _open_gff(gff_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                seq, _, ftype, start_s, end_s, _, strand, _, attr_s = parts[:9]
                if seq != contig:
                    continue
                if ftype not in feature_types:
                    continue
                start = int(start_s)
                end   = int(end_s)
                if end < region_start or start > region_end:
                    continue
                attrs   = _parse_attrs(attr_s)
                lt      = attrs.get("locus_tag", attrs.get("Name", ""))
                product = attrs.get("product", attrs.get("Name", ""))
                product = product.replace("%2C", ",").replace("%3B", ";")
                # `pseudo=true` / `pseudogene` feature → CONFIRMED pseudogene.
                is_ps  = (
                    ftype == "pseudogene"
                    or attrs.get("pseudo", "").lower() == "true"
                )
                is_partial_ann = attrs.get("partial", "").lower() == "true"
                genes.append({
                    "contig":                seq,
                    "start":                 start,
                    "end":                   end,
                    "strand":                strand,
                    "feature_type":          ftype,
                    "locus_tag":             lt,
                    "product":               product,
                    "is_pseudo":             is_ps,
                    "is_partial_annotation": is_partial_ann,
                    "attrs":                 attr_s,
                    "family":                None,
                })
    except FileNotFoundError:
        pass
    genes.sort(key=lambda g: g["start"])
    return genes
