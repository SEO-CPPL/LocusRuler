#!/usr/bin/env python3
"""Product and gene family classification."""

import csv
import json
from pathlib import Path
from typing import Optional


def _read_rules_csv(path: Path) -> list[tuple[str, list[str]]]:
    rules: dict[str, list[str]] = {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rules.setdefault(row["family"], []).append(row["keyword"])
    return list(rules.items())


def load_rules(
    rules_file: str,
    locus_id: str = None,
    config_dir: Path = None,
) -> list[tuple[str, list[str]]]:
    """Load classification rules."""
    if locus_id and config_dir:
        local_csv = Path(config_dir) / f"{locus_id}_rules.csv"
        if local_csv.exists():
            print(f"[config] Using local classification rules: {local_csv.name}")
            return _read_rules_csv(local_csv)

    if not rules_file:
        return []
    path = Path(rules_file)
    if not path.exists():
        alt = path.with_suffix(".csv" if path.suffix == ".json" else ".json")
        if not alt.exists():
            return []
        path = alt
    if path.suffix == ".csv":
        return _read_rules_csv(path)
    data = json.loads(path.read_text())
    return [(r["family"], r["keywords"]) for r in data["family_rules"]]


def classify_product(product: str, rules: list[tuple[str, list[str]]]) -> Optional[str]:
    """Return the first matching family label, or None."""
    p = product.lower()
    for family, kws in rules:
        if any(kw.lower() in p for kw in kws):
            return family
    return None


def classify_genes(genes: list[dict], rules: list[tuple[str, list[str]]]) -> None:
    """Add 'family' field to each gene dict (in-place)."""
    for g in genes:
        g["family"] = classify_product(g["product"], rules)
