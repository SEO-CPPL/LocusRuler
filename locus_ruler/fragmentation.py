#!/usr/bin/env python3
"""Multi-contig fragmentation annotation."""

from __future__ import annotations

from typing import Any

# ── Public type constants ──────────────────────────────────────────────────────
BIOLOGICAL_SPLIT        = "BIOLOGICAL_SPLIT"
ASSEMBLY_SPLIT_CANDIDATE = "ASSEMBLY_SPLIT_CANDIDATE"
DIVERGENT_FRAGMENTED    = "DIVERGENT_FRAGMENTED"
UNCERTAIN               = "UNCERTAIN"

# NCBI assembly_level strings → tier
_DRAFT_LEVELS  = {"scaffold", "contig"}
_CLOSED_LEVELS = {"complete genome", "chromosome"}


# ── Internal helpers ──────────────────────────────────────────────────────
def _assembly_tier(assembly_level: str) -> str:
    """Map NCBI assembly_level string to 'draft', 'closed', or 'unknown'."""
    al = (assembly_level or "").strip().lower()
    if al in _DRAFT_LEVELS:
        return "draft"
    if al in _CLOSED_LEVELS:
        return "closed"
    return "unknown"


def piece_edge_info(piece: dict[str, Any], edge_margin_bp: int) -> tuple[bool | None, int | None]:
    """Return (is_edge_proximal, min_edge_distance_bp) for one piece."""
    try:
        ss   = int(piece["sstart"])
        se   = int(piece["send"])
        clen = int(piece.get("_contig_len") or 0)
    except (KeyError, TypeError, ValueError):
        return None, None

    if clen <= 0:
        return None, None

    s_lo = min(ss, se)
    s_hi = max(ss, se)

    left_dist  = s_lo - 1        # bp from contig start (0-based distance)
    right_dist = clen - s_hi     # bp from contig end

    min_dist = min(left_dist, right_dist)
    return min_dist <= edge_margin_bp, min_dist


def q_pair_max_overlap(pieces: list[dict[str, Any]]) -> float:
    """Max pairwise q-range overlap as a fraction of the shorter piece q-length."""
    if len(pieces) < 2:
        return 0.0

    intervals: list[tuple[int, int]] = []
    for p in pieces:
        try:
            qs = int(p["qstart"])
            qe = int(p["qend"])
            intervals.append((min(qs, qe), max(qs, qe)))
        except (KeyError, TypeError, ValueError):
            intervals.append((0, 0))

    max_ovlp = 0.0
    n = len(intervals)
    for i in range(n):
        for j in range(i + 1, n):
            lo_i, hi_i = intervals[i]
            lo_j, hi_j = intervals[j]
            overlap = max(0, min(hi_i, hi_j) - max(lo_i, lo_j))
            shorter = min(hi_i - lo_i, hi_j - lo_j)
            if shorter > 0:
                max_ovlp = max(max_ovlp, overlap / shorter)

    return max_ovlp


# ── Public entry point ──────────────────────────────────────────────────────
def classify_fragmentation(
    pieces: list[dict[str, Any]],
    assembly_level: str,
    coverage: float,
    settings: dict[str, Any],
) -> dict[str, Any]:
    """Classify a FRAGMENTED locus and return annotation fields."""
    q_comp_thr   = float(settings.get("q_comp_threshold",    0.70))
    q_ovlp_thr   = float(settings.get("q_overlap_threshold", 0.30))
    edge_margin  = int(  settings.get("edge_margin_bp",       500))

    contig_count = len({p.get("sseqid", "") for p in pieces})

    # Per-piece edge info
    edge_info = [piece_edge_info(p, edge_margin) for p in pieces]
    prox_flags = [ep for ep, _ in edge_info]
    dist_vals  = [d  for _, d  in edge_info]

    # same-contig FRAGMENTED → always BIOLOGICAL_SPLIT
    if contig_count <= 1:
        return _result(BIOLOGICAL_SPLIT, 0.0, prox_flags, dist_vals, assembly_level)

    q_ovlp   = q_pair_max_overlap(pieces)
    all_ep   = all(ep is True  for ep in prox_flags)
    any_ep   = any(ep is True  for ep in prox_flags)

    if coverage < q_comp_thr or q_ovlp > q_ovlp_thr:
        ftype = DIVERGENT_FRAGMENTED
    elif all_ep:
        ftype = ASSEMBLY_SPLIT_CANDIDATE
    elif not any_ep:
        ftype = BIOLOGICAL_SPLIT
    else:
        ftype = UNCERTAIN

    return _result(ftype, q_ovlp, prox_flags, dist_vals, assembly_level)


def _result(
    ftype: str,
    q_ovlp: float,
    prox_flags: list[bool | None],
    dist_vals: list[int | None],
    assembly_level: str,
) -> dict[str, Any]:
    all_ep = all(ep is True for ep in prox_flags)
    return {
        "fragmentation_type":  ftype,
        "q_pair_max_overlap":  q_ovlp,
        "all_edge_proximal":   all_ep,
        "contig_edge_support": "Y" if all_ep else ("N" if prox_flags else ""),
        "assembly_tier":       _assembly_tier(assembly_level),
        "piece_edge_proximal": prox_flags,
        "piece_edge_distance": dist_vals,
    }
