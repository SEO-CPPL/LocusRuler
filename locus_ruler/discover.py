#!/usr/bin/env python3
"""Locus discovery by browsing a reference genome."""

from __future__ import annotations

import os
import sqlite3
import sys
from pathlib import Path

# Color only when a person is watching.
_USE_COLOR = (
    sys.stdout.isatty()
    and os.environ.get("NO_COLOR") is None
    and os.environ.get("TERM", "") != "dumb"
)
_CODES = {
    "bold": "1", "dim": "2",
    "red": "31", "green": "32", "yellow": "33",
    "blue": "34", "magenta": "35", "cyan": "36",
}


def paint(text: str, *styles: str) -> str:
    """Wrap text in ANSI styles, or return it unchanged when not a terminal."""
    if not _USE_COLOR or not styles:
        return text
    codes = ";".join(_CODES[s] for s in styles if s in _CODES)
    return f"\033[{codes}m{text}\033[0m" if codes else text

# Hits further apart than this start a new candidate.
DEFAULT_GROUP_GAP_BP = 20_000
CONTEXT_GENES = 12


def search(db_path: Path, accession: str, term: str) -> list[dict]:
    """Genes matching `term` in product, gene name, or locus_tag."""
    con = sqlite3.connect(str(db_path))
    try:
        rows = con.execute(
            "SELECT locus_tag, contig, start, end, strand, product, gene_name "
            "FROM proteins WHERE genome_acc = ? AND ("
            "  lower(product) LIKE '%' || lower(?) || '%' OR "
            "  lower(COALESCE(gene_name,'')) LIKE '%' || lower(?) || '%' OR "
            "  lower(locus_tag) LIKE '%' || lower(?) || '%') "
            "ORDER BY contig, start",
            (accession, term, term, term),
        ).fetchall()
    finally:
        con.close()
    return [
        {"locus_tag": r[0], "contig": r[1], "start": r[2], "end": r[3],
         "strand": r[4], "product": r[5] or "", "gene_name": r[6] or ""}
        for r in rows
    ]


def group(hits: list[dict], max_gap: int = DEFAULT_GROUP_GAP_BP) -> list[dict]:
    """Collapse scattered hits into candidate neighborhoods."""
    groups: list[dict] = []
    for hit in hits:
        if (groups and groups[-1]["contig"] == hit["contig"]
                and hit["start"] - groups[-1]["end"] <= max_gap):
            current = groups[-1]
            current["end"] = max(current["end"], hit["end"])
            current["hits"].append(hit)
        else:
            groups.append({"contig": hit["contig"], "start": hit["start"],
                           "end": hit["end"], "hits": [hit]})
    return groups


def neighborhood(db_path: Path, accession: str, contig: str,
                  start: int, end: int, pad_genes: int = CONTEXT_GENES,
                  pad_before: int | None = None,
                  pad_after: int | None = None) -> list[dict]:
    """Genes across a span, plus context on each side.

    `pad_genes` sets both sides; `pad_before`/`pad_after` override one side
    each, so a caller can widen only the direction it needs.
    """
    pad_before = pad_genes if pad_before is None else pad_before
    pad_after = pad_genes if pad_after is None else pad_after
    con = sqlite3.connect(str(db_path))
    try:
        before = con.execute(
            "SELECT locus_tag, contig, start, end, strand, product FROM proteins "
            "WHERE genome_acc=? AND contig=? AND end < ? ORDER BY start DESC LIMIT ?",
            (accession, contig, start, pad_before)).fetchall()[::-1]
        inside = con.execute(
            "SELECT locus_tag, contig, start, end, strand, product FROM proteins "
            "WHERE genome_acc=? AND contig=? AND end >= ? AND start <= ? ORDER BY start",
            (accession, contig, start, end)).fetchall()
        after = con.execute(
            "SELECT locus_tag, contig, start, end, strand, product FROM proteins "
            "WHERE genome_acc=? AND contig=? AND start > ? ORDER BY start LIMIT ?",
            (accession, contig, end, pad_after)).fetchall()
    finally:
        con.close()

    def _rows(raw, zone):
        return [{"locus_tag": r[0], "contig": r[1], "start": r[2], "end": r[3],
                 "strand": r[4], "product": r[5] or "", "zone": zone}
                for r in raw]

    return _rows(before, "before") + _rows(inside, "hit") + _rows(after, "after")


def genome_label(db_path: Path, accession: str) -> str:
    con = sqlite3.connect(str(db_path))
    try:
        row = con.execute(
            "SELECT species, strain, assembly_level FROM genomes WHERE accession=?",
            (accession,)).fetchone()
    finally:
        con.close()
    if not row:
        return accession
    species, strain, level = (row[0] or ""), (row[1] or ""), (row[2] or "")
    return f"{species} {strain}".strip() + (f"  [{level}]" if level else "")


SHOWN_HITS = 3


def format_candidates(groups: list[dict]) -> list[str]:
    """One block per candidate, listing the genes that matched."""
    lines = []
    for index, g in enumerate(groups, start=1):
        span = g["end"] - g["start"]
        block = [
            f"  {paint(f'[{index}]', 'bold', 'cyan')} "
            f"{g['contig']}:{g['start']:,}-{g['end']:,}  "
            + paint(f"({span:,} bp, {len(g['hits'])} hit"
                    f"{'s' if len(g['hits']) != 1 else ''})", "dim")
        ]
        for hit in g["hits"][:SHOWN_HITS]:
            tag = paint(f"{hit['locus_tag']:<18s}", "cyan")
            block.append(f"      {tag} {hit['product'][:46]}")
        rest = len(g["hits"]) - SHOWN_HITS
        if rest > 0:
            block.append(paint(f"      ... and {rest} more", "dim"))
        lines.append("\n".join(block))
    return lines


def format_neighborhood(genes: list[dict]) -> list[str]:
    """One line per gene."""
    lines = []
    for index, gene in enumerate(genes, start=1):
        hit = gene["zone"] == "hit"
        mark = paint("*", "bold", "yellow") if hit else " "
        tag = paint(f"{gene['locus_tag']:<18s}", "cyan")
        coords = paint(f"{gene['start']:>9,}-{gene['end']:<9,} {gene['strand']}", "dim")
        product = gene["product"][:52]
        lines.append(
            f"  {paint(f'{index:>3}', 'bold') if hit else f'{index:>3}'}{mark} "
            f"{tag} {coords}  "
            f"{paint(product, 'yellow') if hit else product}"
        )
    return lines


def search_any_genome(db_path: Path, term: str, limit: int = 30) -> list[dict]:
    """Like `search`, but across every genome in the target, not just one."""
    con = sqlite3.connect(str(db_path))
    try:
        rows = con.execute(
            "SELECT p.locus_tag, p.genome_acc, p.contig, p.start, p.end, "
            "p.strand, p.product, g.species, g.strain "
            "FROM proteins p LEFT JOIN genomes g ON g.accession = p.genome_acc "
            "WHERE lower(p.product) LIKE '%' || lower(?) || '%' OR "
            "  lower(COALESCE(p.gene_name,'')) LIKE '%' || lower(?) || '%' OR "
            "  lower(p.locus_tag) LIKE '%' || lower(?) || '%' "
            "ORDER BY p.genome_acc, p.contig, p.start LIMIT ?",
            (term, term, term, limit + 1),
        ).fetchall()
    finally:
        con.close()
    return [
        {"locus_tag": r[0], "genome_acc": r[1], "contig": r[2], "start": r[3],
         "end": r[4], "strand": r[5], "product": r[6] or "",
         "species": r[7] or "", "strain": r[8] or ""}
        for r in rows
    ]


def find_genomes(db_path: Path, term: str) -> list[dict]:
    """Genomes whose accession, species or strain contains `term`."""
    con = sqlite3.connect(str(db_path))
    try:
        rows = con.execute(
            "SELECT accession, species, strain, assembly_level FROM genomes "
            "WHERE lower(accession) LIKE '%' || lower(?) || '%' "
            "   OR lower(COALESCE(species,'')) LIKE '%' || lower(?) || '%' "
            "   OR lower(COALESCE(strain,'')) LIKE '%' || lower(?) || '%' "
            "ORDER BY CASE assembly_level WHEN 'Complete Genome' THEN 0 ELSE 1 END,"
            "         species, strain",
            (term, term, term),
        ).fetchall()
    finally:
        con.close()
    return [
        {"accession": r[0], "species": r[1] or "", "strain": r[2] or "",
         "assembly_level": r[3] or ""}
        for r in rows
    ]
