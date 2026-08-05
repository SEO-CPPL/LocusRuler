#!/usr/bin/env python3
"""Within-status hierarchical clustering for the heatmap."""

from typing import Any


def _genome_break_set(res: dict, bin_bp: int = 50) -> set[tuple[int, str]]:
    """Set of (q_pos_bin, ref_gene) pairs representing this genome's
    internal-break fingerprint.  ``bin_bp`` quantises break positions so
    nearby breaks (same biological event, slightly different alignment
    edge) collapse into a single marker."""
    out: set[tuple[int, str]] = set()
    for piece_breaks in (res.get("_internal_breaks") or []):
        for b in piece_breaks or []:
            q_lo = int(b.get("q_gap_lo") or 0)
            bin_pos = (q_lo // bin_bp) * bin_bp
            gene = b.get("ref_gene") or "(intergenic)"
            out.add((bin_pos, gene))
    return out


def informative_marker_positions(
    ruler_results: dict[str, dict],
    bin_bp: int = 50,
    min_genomes: int = 3,
) -> list[tuple[int, str]]:
    """Return marker positions that recur in ≥ ``min_genomes`` genomes."""
    counts: dict[tuple[int, str], int] = {}
    for res in ruler_results.values():
        for key in _genome_break_set(res, bin_bp):
            counts[key] = counts.get(key, 0) + 1
    return sorted([k for k, n in counts.items() if n >= min_genomes])


def cluster_within_group(
    group: list[tuple[Any, dict, dict, str]],
    markers: list[tuple[int, str]],
) -> list[tuple[Any, dict, dict, str]]:
    """Sort a ``(acc, res, meta, status)`` list by hierarchical clustering
    on the binary marker presence matrix.  See module docstring for the
    algorithm.  Returns ``group`` unchanged when clustering is not
    applicable (small group / no variation / scipy missing)."""
    if len(group) < 3 or not markers:
        return group
    try:
        import numpy as np
        from scipy.spatial.distance import pdist
        from scipy.cluster.hierarchy import (
            linkage, leaves_list, optimal_leaf_ordering,
        )
    except ImportError:
        return group

    matrix = []
    for acc, res, meta, status in group:
        s = _genome_break_set(res)
        matrix.append([1 if m in s else 0 for m in markers])
    arr = np.array(matrix, dtype=float)
    if arr.sum() == 0 or arr.std() == 0:
        return group   # no signal to cluster on
    try:
        d = pdist(arr, metric="jaccard")
        if d.sum() == 0:
            return group   # all-zero distances — leave order intact
        Z = linkage(d, method="average")
        Z = optimal_leaf_ordering(Z, d)
        order = leaves_list(Z)
        ordered = [group[i] for i in order]

        # Secondary sort within tied fingerprint blocks.
        def _info_fp(item):
            s = _genome_break_set(item[1])
            return tuple(1 if m in s else 0 for m in markers)

        def _block_key(item):
            acc_, _res, meta_, _status = item
            return (
                tuple(sorted(_genome_break_set(item[1]))),
                (meta_ or {}).get("species", ""),
                (meta_ or {}).get("strain", ""),
                acc_,
            )

        final = []
        i = 0
        while i < len(ordered):
            j = i + 1
            fp_i = _info_fp(ordered[i])
            while j < len(ordered) and _info_fp(ordered[j]) == fp_i:
                j += 1
            block = ordered[i:j]
            block.sort(key=_block_key)
            final.extend(block)
            i = j
        return final
    except Exception:
        return group
