#!/usr/bin/env python3
"""Anchor coverage measured against the cohort's aligned window."""

from __future__ import annotations

from collections import Counter

# A consensus built from a handful of genomes is noise, not a cohort.
MIN_COHORT = 10
# Fraction of confident hits that must cover a query position for it to count.
DEPTH_FRACTION = 0.5

Span = tuple[int, int, int]     # (query_lo, query_hi, positions_in_consensus)


def _hit_qrange(hit: dict) -> tuple[int, int] | None:
    try:
        qs, qe = int(hit["qstart"]), int(hit["qend"])
    except (KeyError, ValueError, TypeError):
        return None
    return (qs, qe) if qs <= qe else (qe, qs)


def consensus_span(
    confident_hits: list[dict],
    depth_fraction: float = DEPTH_FRACTION,
    min_cohort: int = MIN_COHORT,
) -> Span | None:
    """Query window that at least `depth_fraction` of confident hits cover."""
    if len(confident_hits) < min_cohort:
        return None
    depth: Counter = Counter()
    for hit in confident_hits:
        span = _hit_qrange(hit)
        if span is None:
            continue
        for position in range(span[0], span[1] + 1):
            depth[position] += 1
    if not depth:
        return None
    needed = max(1, int(len(confident_hits) * depth_fraction))
    covered = sorted(position for position, d in depth.items() if d >= needed)
    if not covered:
        return None
    return covered[0], covered[-1], len(covered)


def build_spans(
    blast_hits: dict[str, dict[str, list[dict]]],
    qlens: dict[str, int],
    min_coverage: float,
    min_identity: float,
    depth_fraction: float = DEPTH_FRACTION,
    min_cohort: int = MIN_COHORT,
) -> dict[str, Span]:
    """Consensus window per anchor, from the hits that already pass the gate."""
    confident: dict[str, list[dict]] = {}
    for per_anchor in blast_hits.values():
        for tag, hits in per_anchor.items():
            qlen = qlens.get(tag)
            if not qlen:
                continue
            best = None
            best_score = -1.0
            for hit in hits or []:
                span = _hit_qrange(hit)
                if span is None:
                    continue
                try:
                    identity = float(hit.get("pident", 0))
                    score = float(hit.get("bitscore", 0))
                except (TypeError, ValueError):
                    continue
                coverage = (span[1] - span[0] + 1) / qlen
                if coverage >= min_coverage and identity >= min_identity and score > best_score:
                    best, best_score = hit, score
            if best is not None:
                confident.setdefault(tag, []).append(best)

    spans: dict[str, Span] = {}
    for tag, hits in confident.items():
        span = consensus_span(hits, depth_fraction, min_cohort)
        if span is not None:
            spans[tag] = span
    return spans


def build_spans_by_group(
    blast_hits: dict[str, dict[str, list[dict]]],
    qlens: dict[str, int],
    group_of: dict[str, str],
    min_coverage: float,
    min_identity: float,
    depth_fraction: float = DEPTH_FRACTION,
    min_cohort: int = MIN_COHORT,
) -> dict[str, dict[str, Span]]:
    """Consensus window per anchor within each group of `group_of`.

    A cohort drawn from several lineages has no single answer to how much of an
    anchor aligns. Whichever lineage brings the most genomes decides the window,
    and the others are then judged against a standard their orthologue was never
    going to meet. Grouping the consensus keeps each lineage measured against
    the way its own members align, and leaves a single-lineage run untouched.
    """
    grouped: dict[str, dict[str, dict[str, list[dict]]]] = {}
    for accession, per_anchor in blast_hits.items():
        group = group_of.get(accession) or ""
        grouped.setdefault(group, {})[accession] = per_anchor
    return {
        group: build_spans(hits, qlens, min_coverage, min_identity,
                           depth_fraction, min_cohort)
        for group, hits in grouped.items()
    }


def cohort_coverage(hit: dict, span: Span) -> float:
    """Fraction of the consensus window this hit covers (0.0 when disjoint)."""
    qrange = _hit_qrange(hit)
    if qrange is None:
        return 0.0
    lo, hi, length = span
    if length <= 0:
        return 0.0
    overlap = min(qrange[1], hi) - max(qrange[0], lo) + 1
    return max(0, overlap) / length


def describe(spans: dict[str, Span], qlens: dict[str, int]) -> list[str]:
    """One line per anchor whose consensus is shorter than the protein."""
    lines = []
    for tag, (lo, hi, length) in sorted(spans.items()):
        qlen = qlens.get(tag) or 0
        if qlen and length < qlen * 0.95:
            lines.append(
                f"{tag}: consensus {lo}-{hi} ({length} aa of {qlen} aa; "
                f"{length / qlen:.0%} of the anchor aligns across the cohort)"
            )
    return lines
