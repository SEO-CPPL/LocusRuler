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


def reach_rescued(
    span: tuple[str, int, int, int, int],
    rescued: list[dict] | None,
    reference_bp: int,
) -> tuple[str, int, int, int, int]:
    """Widen the span to a domain-rescued gene lying just outside it.

    The span comes from the accepted pieces, which by definition never held the
    rescued gene, so anything between the two would otherwise go unreported.
    The reference length caps the reach: a homolog on the far side of the
    chromosome is a separate locus, not this cassette's edge.
    """
    if not rescued:
        return span
    contig, lo, hi, gap_left, gap_right = span
    for gene in rescued:
        if gene.get("gene_contig") != contig:
            continue
        g_lo, g_hi = _int(gene.get("gene_start")), _int(gene.get("gene_end"))
        if g_hi <= g_lo:
            continue
        if 0 < lo - g_hi <= reference_bp:
            lo = g_lo
        elif 0 < g_lo - hi <= reference_bp:
            hi = g_hi
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


def contig_lengths(db_path: Path, wanted: dict[str, str]) -> dict[str, int]:
    """Last annotated base of the named contig, per genome."""
    out: dict[str, int] = {}
    db_path = Path(db_path)
    if not db_path.exists() or not wanted:
        return out
    con = sqlite3.connect(str(db_path))
    try:
        for acc, contig in wanted.items():
            row = con.execute(
                "SELECT MAX(end) FROM proteins WHERE genome_acc = ? AND contig = ?",
                (acc, contig),
            ).fetchone()
            if row and row[0]:
                out[acc] = int(row[0])
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
    membership_source: str = "flank_extended_piece",
) -> list[dict]:
    """Add what the extension reached, without dropping anything already found.

    A source row may already carry its own family/state (domain-rescued genes
    do); those are kept instead of falling back to unassigned. A gene already
    listed as context is promoted in place rather than skipped, so naming it
    is what settles which zone it belongs to.
    """
    by_tag = {g.get("locus_tag"): g for g in ordered if g.get("locus_tag")}
    extras = []
    for gene in span_genes:
        row = dict(gene)
        family = row.get("family") or "unassigned"
        row.update({
            "zone": zone,
            "family": family,
            "state": row.get("state", ""),
            "label": row.get("label") or f"[{family}]",
            "membership_source": membership_source,
        })
        already = by_tag.get(row.get("locus_tag"))
        if already is None:
            extras.append(row)
        elif family != "unassigned":
            already.update({k: row[k] for k in
                            ("zone", "family", "state", "label", "membership_source")})
    if not extras:
        return ordered

    def _start(gene, default=0):
        return _int(gene.get("gene_start"), default)

    result = list(ordered)
    for extra in extras:
        contig = extra.get("gene_contig")
        here = [g for g in result if g.get("gene_contig") == contig]
        # The list may run backwards along the contig, so follow whichever way it goes
        descending = len(here) >= 2 and _start(here[0]) > _start(here[-1])
        target = len(result)
        for index, gene in enumerate(result):
            if gene.get("gene_contig") != contig:
                continue
            position = _start(gene, 0 if descending else (1 << 62))
            if (position < _start(extra)) if descending else (position > _start(extra)):
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
