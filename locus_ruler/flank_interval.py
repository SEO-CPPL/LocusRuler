#!/usr/bin/env python3
"""Cassette extension out to its flanks."""

from __future__ import annotations

import sqlite3
from pathlib import Path

FLANK_LEFT = "flank_left"
FLANK_RIGHT = "flank_right"


def _int(value, default=0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def piece_bounds(pieces: dict[str, dict]) -> tuple[str, int, int] | None:
    """Outer span of the accepted pieces, as (contig, lo, hi)."""
    spans = []
    for piece in pieces.values():
        contig = str(piece.get("contig") or "")
        lo, hi = _int(piece.get("s_lo")), _int(piece.get("s_hi"))
        if contig and hi > lo:
            spans.append((contig, min(lo, hi), max(lo, hi)))
    if not spans:
        return None
    contigs = {c for c, _, _ in spans}
    if len(contigs) != 1:
        return None
    contig = contigs.pop()
    return contig, min(lo for _, lo, _ in spans), max(hi for _, _, hi in spans)


def extend_to_flanks(
    rows: list[dict],
    pieces: dict[str, dict],
    reference_bp: int,
) -> tuple[str, int, int, int, int] | None:
    """Widen the piece span toward the flanks, capped at the reference length."""
    span = piece_bounds(pieces)
    if span is None or reference_bp <= 0:
        return None
    contig, lo, hi = span

    left = [r for r in rows
            if r.get("zone") == FLANK_LEFT and r.get("gene_contig") == contig]
    right = [r for r in rows
             if r.get("zone") == FLANK_RIGHT and r.get("gene_contig") == contig]

    gap_left = gap_right = -1
    if left:
        # innermost flank on the left: the one whose end sits closest below lo
        inner = max((r for r in left if _int(r.get("gene_end")) <= lo),
                    key=lambda r: _int(r.get("gene_end")), default=None)
        if inner is not None:
            gap = lo - _int(inner.get("gene_end"))
            if 0 <= gap <= reference_bp:
                lo = _int(inner.get("gene_end"))
                gap_left = gap
    if right:
        inner = min((r for r in right if _int(r.get("gene_start"), 1 << 62) >= hi),
                    key=lambda r: _int(r.get("gene_start"), 1 << 62), default=None)
        if inner is not None:
            gap = _int(inner.get("gene_start")) - hi
            if 0 <= gap <= reference_bp:
                hi = _int(inner.get("gene_start"))
                gap_right = gap
    if gap_left < 0 and gap_right < 0:
        return None
    return contig, lo, hi, gap_left, gap_right


def genes_in_span(
    db_path: Path,
    wanted: dict[str, tuple[str, int, int, int, int]],
) -> dict[str, list[dict]]:
    """Annotated genes strictly inside each extended span, from the genome DB."""
    out: dict[str, list[dict]] = {}
    db_path = Path(db_path)
    if not db_path.exists() or not wanted:
        return out
    con = sqlite3.connect(str(db_path))
    try:
        for acc, (contig, lo, hi, _gl, _gr) in wanted.items():
            # Overlap, not containment: a gene running a few bases past the flank still counts.
            rows = con.execute(
                "SELECT locus_tag, contig, start, end, strand, product "
                "FROM proteins WHERE genome_acc = ? AND contig = ? "
                "AND end > ? AND start < ? ORDER BY start",
                (acc, contig, lo, hi),
            ).fetchall()
            out[acc] = [
                {
                    "locus_tag": r[0], "gene_contig": r[1],
                    "gene_start": str(r[2]), "gene_end": str(r[3]),
                    "strand": r[4] or "", "product": r[5] or "",
                }
                for r in rows
            ]
    finally:
        con.close()
    return out


def reference_rank(anchors: list[dict]) -> dict[str, int]:
    """Family -> position in the reference locus, for orientation only."""
    rank: dict[str, int] = {}
    for index, anchor in enumerate(anchors):
        family = (anchor.get("family") or anchor.get("locus_tag") or "").strip()
        if family and family not in rank:
            rank[family] = index
    return rank


def orient(genes: list[dict], rank: dict[str, int]) -> list[dict]:
    """Put the span in reference order."""
    positions = [
        rank[(g.get("family") or "").strip()]
        for g in genes
        if (g.get("family") or "").strip() in rank
    ]
    if len(positions) < 2:
        return genes
    forward = sum(b > a for a, b in zip(positions, positions[1:]))
    backward = sum(b < a for a, b in zip(positions, positions[1:]))
    return list(reversed(genes)) if backward > forward else genes


def augment(
    ordered: list[dict],
    span_genes: list[dict],
    zone: str = "cluster",
) -> list[dict]:
    """Add what the extension reached, without dropping anything already found."""
    seen = {g.get("locus_tag") for g in ordered if g.get("locus_tag")}
    extras = []
    for gene in span_genes:
        if gene.get("locus_tag") in seen:
            continue
        row = dict(gene)
        row.update({
            "zone": zone,
            "family": "unassigned",
            "state": "",
            "label": "[unassigned]",
            "membership_source": "flank_extended_piece",
        })
        extras.append(row)
    if not extras:
        return ordered

    def _start(gene, default=0):
        return _int(gene.get("gene_start"), default)

    result = list(ordered)
    for extra in extras:
        contig = extra.get("gene_contig")
        target = len(result)
        for index, gene in enumerate(result):
            if gene.get("gene_contig") != contig:
                continue
            if _start(gene, 1 << 62) > _start(extra):
                target = index
                break
        result.insert(target, extra)
    return result


def describe(wanted: dict, total: int, reference_bp: int) -> str:
    if not wanted:
        return (f"0 of {total} genomes could reach a flank within "
                f"{reference_bp:,} bp; using accepted pieces unchanged")
    both = sum(1 for _, _, _, gl, gr in wanted.values() if gl >= 0 and gr >= 0)
    gaps = sorted(g for _, _, _, gl, gr in wanted.values() for g in (gl, gr) if g >= 0)
    return (f"{len(wanted)} of {total} genomes extended to a flank "
            f"({both} on both sides); reach {gaps[0]:,}-{gaps[-1]:,} bp "
            f"(median {gaps[len(gaps) // 2]:,}), capped at the "
            f"{reference_bp:,} bp reference locus")
