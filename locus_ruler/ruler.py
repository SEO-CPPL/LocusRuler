#!/usr/bin/env python3
"""Locus status determination."""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional, Iterable

from duplication import (
    verify_duplications as _verify_duplications,
    _load_ruler_genome_meta,
)
from hsp_filter import cap_family_hsps as _cap_family_hsps
from weak_hsp_rescue import attach_domain_guided_weak_hsps

# ── Status constants ──────────────────────────────────────────────────────
INTACT          = "INTACT"
PSEUDOGENIZED   = "PSEUDOGENIZED"
PARTIAL_DEL     = "PARTIAL_DEL"
LARGE_DEL       = "LARGE_DEL"
FRAGMENTED      = "FRAGMENTED"
DUPLICATED      = "DUPLICATED"
ABSENT          = "ABSENT"
ANNOTATION_GAP  = "ANNOTATION_GAP"
UNKNOWN         = "UNKNOWN"

SUMMARY_COLS = [
    "genome_acc", "species", "strain", "assembly_level", "locus_id", "status",
    "synteny_class", "synteny_case",
    "L_state", "R_state",
    "L_fingerprint", "R_fingerprint",
    "L2_gene", "L1_gene", "R1_gene", "R2_gene",
    "inversion", "orientation",
    "coverage",
    "L2_L1_bp", "flank_bp", "R1_R2_bp",
    "cluster_bp",
    "expected_L2_L1_bp", "expected_flank_bp", "expected_R1_R2_bp",
    "expected_cluster_bp",
    "order_flags",
    "genes_found", "genes_expected", "pseudo_count",
    "linked_locus", "contig", "notes",
    # Multi-contig fragmentation annotation, populated for FRAGMENTED rows only.
    "fragmentation_type", "contig_edge_support", "q_pair_max_overlap",
]


# ── Generic helpers ──────────────────────────────────────────────────────
def _union_length(intervals: list[tuple[int, int]]) -> int:
    """Sum of merged lengths for a list of (lo, hi) intervals."""
    if not intervals:
        return 0
    valid = [(lo, hi) for lo, hi in intervals if hi >= lo]
    if not valid:
        return 0
    valid.sort()
    merged = [list(valid[0])]
    for s, e in valid[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return sum(e - s + 1 for s, e in merged)


# ── Hit utilities ──────────────────────────────────────────────────────
def _best_hit(rows: list[dict], min_id: float) -> Optional[dict]:
    """Return the highest-bitscore hit that meets the identity threshold."""
    candidates = [r for r in rows if float(r.get("pident", 0)) >= min_id]
    if not candidates:
        return None
    return max(candidates, key=lambda r: float(r.get("bitscore", 0)))


def _hit_span(hit: dict) -> tuple[str, int, int]:
    contig = hit["sseqid"]
    s, e = int(hit["sstart"]), int(hit["send"])
    return contig, min(s, e), max(s, e)


# ── Reference gene coordinate mapping ──────────────────────────────────────────────────────
def _reference_gene_coords(locus_cfg: dict) -> list[dict]:
    """
    Return reference cluster genes with coordinates RELATIVE to the cluster
    nucleotide chunk (cluster.fna).

    Aux-role anchors are excluded: they are auxiliary search probes that
    don't participate in cluster geometry / piece assembly.
    """
    auto    = locus_cfg.get("_auto", {})
    anchors = auto.get("anchors", [])
    cl_a    = next((a for a in anchors if a["role"] == "cluster_L"), None)
    cluster_roles = {"cluster_L", "cluster_R", "inner"}
    # Correctly identify the nucleotide chunk origin (absolute genome coord)
    coords = []
    for a in anchors:
        if a["role"] in cluster_roles:
            coords += [int(a["start"]), int(a["end"])]
    if not coords:
        return []
    origin = min(coords)

    genes: list[dict] = []
    for a in anchors:
        if a["role"] not in cluster_roles:
            continue
        s, e = int(a["start"]), int(a["end"])
        genes.append({
            "locus_tag": a["locus_tag"],
            "family":    a.get("family"),
            "role":      a["role"],
            "rel_start": min(s, e) - origin,
            "rel_end":   max(s, e) - origin,
            "exception": bool(a.get("exception", False)),
            "lenient":   bool(a.get("lenient", False)),
            "pfam":       a.get("pfam", ""),
            "pfam_split": bool(a.get("pfam_split", False)),
        })
    return sorted(genes, key=lambda g: g["rel_start"])


def _genes_overlapping(qstart: int, qend: int,
                       ref_genes: list[dict],
                       min_gene_coverage: float = 0.30) -> list[str]:
    lo, hi = min(qstart, qend), max(qstart, qend)
    out: list[str] = []
    for g in ref_genes:
        if g["rel_end"] < lo or g["rel_start"] > hi:
            continue
        ov_lo = max(g["rel_start"], lo)
        ov_hi = min(g["rel_end"],   hi)
        ov_bp = max(0, ov_hi - ov_lo + 1)
        gene_len = max(1, g["rel_end"] - g["rel_start"] + 1)
        if ov_bp / gene_len >= min_gene_coverage:
            out.append(g["locus_tag"])
    return out


def _gene_overlap_details(qstart: int, qend: int,
                          ref_genes: list[dict]) -> list[dict]:
    """Reference-gene overlaps for one HSP query range."""
    lo, hi = min(qstart, qend), max(qstart, qend)
    out: list[dict] = []
    for g in ref_genes:
        if g["rel_end"] < lo or g["rel_start"] > hi:
            continue
        ov_lo = max(g["rel_start"], lo)
        ov_hi = min(g["rel_end"], hi)
        ov_bp = max(0, ov_hi - ov_lo + 1)
        gene_len = max(1, g["rel_end"] - g["rel_start"] + 1)
        if ov_bp > 0:
            out.append({
                "locus_tag": g["locus_tag"],
                "family": g.get("family") or g["locus_tag"],
                "overlap_bp": ov_bp,
                "coverage": ov_bp / gene_len,
                "exception": bool(g.get("exception", False)),
            })
    return out


def _annotate_genes(hsps: list[dict], ref_genes: list[dict],
                    min_gene_coverage: float = 0.30) -> None:
    for h in hsps:
        qs = min(int(h["qstart"]), int(h["qend"]))
        qe = max(int(h["qstart"]), int(h["qend"]))
        h["_genes"] = _genes_overlapping(qs, qe, ref_genes, min_gene_coverage)


# ── Piece formation ──────────────────────────────────────────────────────

def _form_pieces(
    hsps: list[dict],
    max_intra_piece_gap: int,
    min_id: float,
) -> list[list[dict]]:
    """Group HSPs into pieces by contig + strand + gap."""
    filtered = [h for h in hsps if float(h.get("pident", 0)) >= min_id]
    if not filtered:
        return []

    def _strand(h): return "+" if int(h["send"]) > int(h["sstart"]) else "-"
    def _slo(h): return min(int(h["sstart"]), int(h["send"]))
    def _shi(h): return max(int(h["sstart"]), int(h["send"]))

    filtered.sort(key=lambda h: (h["sseqid"], _strand(h), _slo(h)))

    # Pass 1: linear gap-based grouping (existing behavior)
    pieces: list[list[dict]] = []
    for h in filtered:
        if pieces:
            last_group = pieces[-1]
            last       = last_group[-1]
            same_contig = h["sseqid"] == last["sseqid"]
            same_strand = _strand(h) == _strand(last)
            group_s_hi  = max(_shi(x) for x in last_group)
            gap         = _slo(h) - group_s_hi
            if same_contig and same_strand and gap <= max_intra_piece_gap:
                last_group.append(h)
                continue
        pieces.append([h])

    # Pass 2: origin-wrap merge.
    if len(pieces) < 2:
        return pieces

    buckets: dict[tuple[str, str], list[int]] = {}
    contig_lens: dict[str, int] = {}
    for i, grp in enumerate(pieces):
        c = grp[0]["sseqid"]
        s = _strand(grp[0])
        buckets.setdefault((c, s), []).append(i)
        try:
            contig_lens[c] = int(grp[0].get("slen") or 0)
        except (ValueError, TypeError):
            contig_lens[c] = 0

    wrap_merges: list[tuple[int, int]] = []   # (low_idx, high_idx)
    for (contig, strand), idxs in buckets.items():
        clen = contig_lens.get(contig, 0)
        if clen <= 0 or len(idxs) < 2:
            continue
        # Sort by piece low-coord on the subject contig.
        idxs.sort(key=lambda i: min(_slo(h) for h in pieces[i]))
        low_i, high_i = idxs[0], idxs[-1]
        first_lo = min(_slo(h) for h in pieces[low_i])
        last_hi  = max(_shi(h) for h in pieces[high_i])
        origin_gap = (clen - last_hi) + first_lo - 1
        if origin_gap <= max_intra_piece_gap:
            wrap_merges.append((low_i, high_i))

    if not wrap_merges:
        return pieces

    # Apply merges: combine pieces by index, drop the originals
    drop: set[int] = set()
    merged: list[list[dict]] = []
    for low_i, high_i in wrap_merges:
        if low_i in drop or high_i in drop:
            continue
        # Order the merged HSPs: high-coord (end-side) first, then low-coord (start-side).
        merged.append(pieces[high_i] + pieces[low_i])
        drop.add(low_i)
        drop.add(high_i)

    final: list[list[dict]] = [p for i, p in enumerate(pieces) if i not in drop]
    final.extend(merged)
    return final


def _covered_genes_union(piece_hsps: list[dict],
                         ref_genes: list[dict],
                         min_gene_coverage: float = 0.30,
                         exclude_exceptions: bool = False,
                         min_gene_coverage_lenient: Optional[float] = None) -> list[str]:
    """Reference genes covered by the HSP union of this piece."""
    covered: list[str] = []
    for g in ref_genes:
        if exclude_exceptions and g.get("exception"):
            continue
        gene_len = max(1, g["rel_end"] - g["rel_start"] + 1)
        intervals = []
        for h in piece_hsps:
            qs = min(int(h["qstart"]), int(h["qend"]))
            qe = max(int(h["qstart"]), int(h["qend"]))
            ov_lo = max(g["rel_start"], qs)
            ov_hi = min(g["rel_end"],   qe)
            if ov_hi >= ov_lo:
                intervals.append((ov_lo, ov_hi))

        # Per-gene coverage threshold (lenient genes can use a different value)
        eff_cov = min_gene_coverage
        if g.get("lenient") and min_gene_coverage_lenient is not None:
            eff_cov = min_gene_coverage_lenient
        if _union_length(intervals) / gene_len >= eff_cov:
            covered.append(g["locus_tag"])
    return covered


def _gene_union_coverages(piece_hsps: list[dict],
                          ref_genes: list[dict],
                          exclude_exceptions: bool = False) -> dict[str, float]:
    """Per-reference-gene query coverage from the union of HSP ranges."""
    covs: dict[str, float] = {}
    for g in ref_genes:
        if exclude_exceptions and g.get("exception"):
            continue
        gene_len = max(1, g["rel_end"] - g["rel_start"] + 1)
        intervals = []
        for h in piece_hsps:
            qs = min(int(h["qstart"]), int(h["qend"]))
            qe = max(int(h["qstart"]), int(h["qend"]))
            ov_lo = max(g["rel_start"], qs)
            ov_hi = min(g["rel_end"], qe)
            if ov_hi >= ov_lo:
                intervals.append((ov_lo, ov_hi))
        covs[g["locus_tag"]] = _union_length(intervals) / gene_len
    return covs


def _weighted_identity(piece_hsps: list[dict],
                       exclude_exception_regions: bool = False,
                       ref_genes: list[dict] = None) -> float:
    """Mean identity across HSPs."""
    # Simple case
    if not exclude_exception_regions or not ref_genes:
        total_w = sum(int(h.get("length", 1)) for h in piece_hsps)
        if total_w == 0: return 0.0
        weighted = sum(float(h.get("pident", 0)) * int(h.get("length", 1))
                       for h in piece_hsps)
        return weighted / total_w

    # Advanced case: filter HSP segments by exception status
    exc_intervals = [(g["rel_start"], g["rel_end"]) for g in ref_genes if g.get("exception")]

    total_len = 0
    total_pid_len = 0.0

    for h in piece_hsps:
        h_qs = min(int(h["qstart"]), int(h["qend"]))
        h_qe = max(int(h["qstart"]), int(h["qend"]))
        h_len = h_qe - h_qs + 1
        h_pid = float(h.get("pident", 0))

        # Calculate how much of this HSP overlaps with ANY exception
        overlap_exc = 0
        for es, ee in exc_intervals:
            ov_lo = max(h_qs, es)
            ov_hi = min(h_qe, ee)
            if ov_hi >= ov_lo:
                overlap_exc += (ov_hi - ov_lo + 1)

        # Inert length for this HSP
        inert_len = max(0, h_len - overlap_exc)
        if inert_len > 0:
            total_len += inert_len
            total_pid_len += (h_pid * inert_len)

    return (total_pid_len / total_len) if total_len > 0 else 0.0


def _accept_piece(
    piece_hsps: list[dict],
    ref_genes: list[dict],
    min_genes: int,
    single_gene_min_id: float,
    min_gene_coverage: float = 0.30,
    single_gene_min_coverage: float = 0.70,
    single_gene_min_id_lenient: Optional[float] = None,
    single_gene_min_coverage_lenient: Optional[float] = None,
    min_gene_coverage_lenient: Optional[float] = None,
) -> tuple[bool, list[str], Optional[str]]:
    """
    Apply 'annotated but inert' policy: exceptions don't count toward
    min_genes for Tier A or Tier B acceptance.

    Per-gene lenient override: lenient anchors (lenient=TRUE) use relaxed
    thresholds in BOTH Tier A (piece formation) AND Tier B (single-gene
    rescue). The lenient values are 40% identity + 40% coverage (typical),
    recovering highly divergent paralogs without globally relaxing intra-genus
    stringency.

    Tier A lenient: `min_gene_coverage_lenient` (e.g., 0.40) replaces
    `min_gene_coverage` (0.30) for lenient genes' covered-status check.
    Tier B lenient: `single_gene_min_id_lenient` and
    `single_gene_min_coverage_lenient` replace standard thresholds when the
    single covered gene is lenient.
    """
    # 1. Get union of covered genes, EXCLUDING exceptions (lenient-aware)
    real_covered = _covered_genes_union(
        piece_hsps, ref_genes, min_gene_coverage, exclude_exceptions=True,
        min_gene_coverage_lenient=min_gene_coverage_lenient)
    n = len(real_covered)

    # 2. Also get ALL covered genes for metadata
    all_covered = _covered_genes_union(
        piece_hsps, ref_genes, min_gene_coverage, exclude_exceptions=False,
        min_gene_coverage_lenient=min_gene_coverage_lenient)

    if n >= min_genes:
        return True, all_covered, "A"

    # For Tier B (single gene), the one gene must NOT be an exception
    if n >= 1:
        # Check weighted identity of the piece, ignoring exception regions
        pid = _weighted_identity(piece_hsps, exclude_exception_regions=True, ref_genes=ref_genes)
        covs = _gene_union_coverages(piece_hsps, ref_genes, exclude_exceptions=True)
        best_cov = max((covs.get(lt, 0.0) for lt in real_covered), default=0.0)

        # Per-gene lenient threshold: the gene chosen as the best_cov source must be lenient.
        effective_min_id = single_gene_min_id
        effective_min_cov = single_gene_min_coverage
        if real_covered and (single_gene_min_id_lenient is not None
                             or single_gene_min_coverage_lenient is not None):
            lt_lenient = {g["locus_tag"]: bool(g.get("lenient"))
                          for g in ref_genes if g.get("locus_tag")}
            if any(lt_lenient.get(lt, False) for lt in real_covered):
                if single_gene_min_id_lenient is not None:
                    effective_min_id = single_gene_min_id_lenient
                if single_gene_min_coverage_lenient is not None:
                    effective_min_cov = single_gene_min_coverage_lenient

        if pid >= effective_min_id and best_cov >= effective_min_cov:
            return True, all_covered, "B"

    return False, all_covered, None


def _consolidate_piece(piece_hsps: list[dict],
                       q_clip: Optional[tuple[int, int]] = None) -> dict:
    q_starts = [int(h["qstart"]) for h in piece_hsps]
    q_ends   = [int(h["qend"])   for h in piece_hsps]
    s_starts = [int(h["sstart"]) for h in piece_hsps]
    s_ends   = [int(h["send"])   for h in piece_hsps]

    s_lo, s_hi = min(min(s_starts), min(s_ends)), max(max(s_starts), max(s_ends))

    fwd_w = sum(int(h.get("length", 1)) for h in piece_hsps if int(h["send"]) > int(h["sstart"]))
    rev_w = sum(int(h.get("length", 1)) for h in piece_hsps if int(h["send"]) < int(h["sstart"]))
    sstart_val, send_val = (s_lo, s_hi) if fwd_w >= rev_w else (s_hi, s_lo)

    # HSP-level UNION length for query coordinates.
    intervals = []
    for h in piece_hsps:
        lo = min(int(h["qstart"]), int(h["qend"]))
        hi = max(int(h["qstart"]), int(h["qend"]))
        if q_clip is not None:
            lo = max(lo, q_clip[0])
            hi = min(hi, q_clip[1])
            if hi < lo:
                continue
        intervals.append((lo, hi))
    union_len = _union_length(intervals)

    # ── Origin-wrap detection
    contig_lens_seen: set[int] = set()
    for h in piece_hsps:
        try:
            v = int(h.get("slen") or 0)
            if v > 0:
                contig_lens_seen.add(v)
        except (ValueError, TypeError):
            pass
    clen = next(iter(contig_lens_seen), 0)

    subject_intervals: list[list[int]] = []
    if clen > 0 and len(piece_hsps) >= 2:
        # Sort HSP subject ranges by lo and find the largest internal gap.
        coords = sorted([(min(int(h["sstart"]), int(h["send"])),
                          max(int(h["sstart"]), int(h["send"])))
                         for h in piece_hsps])
        biggest_gap = 0
        gap_after = -1
        for i in range(len(coords) - 1):
            g = coords[i + 1][0] - coords[i][1]
            if g > biggest_gap:
                biggest_gap = g
                gap_after = i
        # Treat as wrap only when the origin gap is smaller than the largest internal gap.
        first_lo = coords[0][0]
        last_hi  = coords[-1][1]
        origin_gap = (clen - last_hi) + first_lo - 1
        if biggest_gap > 0 and origin_gap < biggest_gap:
            # Split at the largest gap, listing the intervals end-side first to match _form_pieces.
            end_side_lo = coords[gap_after + 1][0]
            end_side_hi = last_hi
            start_side_lo = first_lo
            start_side_hi = coords[gap_after][1]
            subject_intervals = [
                [end_side_lo,   end_side_hi],
                [start_side_lo, start_side_hi],
            ]
            # sstart/send carry the end-side interval; subject_intervals carries both.
            sstart_val, send_val = (end_side_lo, end_side_hi) if fwd_w >= rev_w \
                                    else (end_side_hi, end_side_lo)

    piece = {
        "sseqid": piece_hsps[0]["sseqid"],
        "sstart": str(sstart_val),
        "send":   str(send_val),
        "qstart": str(min(min(q_starts), min(q_ends))),
        "qend":   str(max(max(q_starts), max(q_ends))),
        "length": str(union_len),  # unioned query length
        "pident": str(round(_weighted_identity(piece_hsps), 1)),
        "_genes": [],
        "_constituent_hsp_count": len(piece_hsps),
    }
    if subject_intervals:
        piece["subject_intervals"] = subject_intervals
        piece["_wraps_origin"] = True
    if clen > 0:
        piece["_contig_len"] = clen
    return piece


def _cluster_orientation(pieces: list[dict]) -> int:
    if not pieces: return 0
    forward = reverse = 0
    for p in pieces:
        s_dir = 1 if int(p["send"]) > int(p["sstart"]) else -1
        w     = int(p.get("length", 1))
        if s_dir > 0: forward += w
        else: reverse += w
    if forward == reverse: return 0
    return 1 if forward > reverse else -1


def _flank_orientation(res: dict) -> int:
    """Cluster orientation RELATIVE TO THE FLANKING GENES.

    Returns +1 (co-oriented with the reference synteny), -1 (truly inverted),
    or 0 when it cannot be determined (zero/ambiguous cluster orientation,
    multi-contig locus, a missing flank, or flanks that do not bracket the
    cluster).

    ``_cluster_orientation`` alone only tells whether the cluster aligns
    forward or reverse to the reference ``cluster.fna`` — that is the strand
    the contig happens to be assembled/annotated on, NOT a biological event.
    A genome whose whole neighborhood (flanks + cluster) was deposited on
    the opposite strand has *both* the cluster orientation and the flank L->R
    coordinate trend flipped, so their relationship is preserved (-> +1). A
    true local inversion flips the cluster while the flanking synteny keeps
    its arrangement (-> -1).
    """
    orient = res.get("_cluster_orientation", 0) or 0
    if not orient:
        return 0
    # Only meaningful when the cluster and its flanks share one contig.
    if (res.get("contig_count") or 0) != 1:
        return 0
    left  = res.get("_L1_lo")
    if left is None:
        left = res.get("_L2_lo")
    right = res.get("_R1_lo")
    if right is None:
        right = res.get("_R2_lo")
    if left is None or right is None:
        return 0
    cmin, cmax = res.get("_cluster_pos_min"), res.get("_cluster_pos_max")
    if cmin is None or cmax is None:
        return 0
    flank_trend = 1 if right > left else -1
    # The flanks must actually bracket the cluster.
    lo_b, hi_b = (left, right) if flank_trend > 0 else (right, left)
    if not (lo_b <= cmin and cmax <= hi_b):
        return 0
    cluster_dir = 1 if orient > 0 else -1
    # In the reference, cluster.fna-forward runs from the L-flank side toward the R-flank side.
    return -1 if cluster_dir != flank_trend else 1


def _flank_relative_inversion(res: dict) -> str:
    """'Y' if the cluster is inverted relative to flanking synteny, else ''
    (co-oriented OR undeterminable). Thin wrapper over ``_flank_orientation``."""
    return "Y" if _flank_orientation(res) < 0 else ""


def _resolve_cluster_blast_cfg(settings: dict, locus_cfg: dict,
                               target_species: str = "") -> dict:
    """Merge cluster_blast config and enforce standard single-gene rescue."""
    cbcfg = {**settings.get("cluster_blast", {}),
             **locus_cfg.get("cluster_blast", {})}
    ref_genus = (locus_cfg.get("_auto", {}).get("reference_genus") or "").strip()
    target_genus = target_species.strip().split()[0] if target_species else ""
    same_genus = bool(ref_genus and target_genus and ref_genus == target_genus)

    if same_genus:
        std_id = float(settings.get("cluster_blast", {}).get(
            "single_gene_min_identity_intra", 90))
    else:
        std_id = float(settings.get("cluster_blast", {}).get(
            "single_gene_min_identity_inter", 75))
    std_cov = float(settings.get("cluster_blast", {}).get(
        "single_gene_min_coverage", 0.70))

    # Locus configs may be stricter, but should not lower the standardized single-gene rescue floor.
    cbcfg["single_gene_min_identity"] = max(
        float(cbcfg.get("single_gene_min_identity", std_id)), std_id)
    cbcfg["single_gene_min_coverage"] = max(
        float(cbcfg.get("single_gene_min_coverage", std_cov)), std_cov)
    cbcfg["_single_gene_identity_mode"] = "intra-genus" if same_genus else "inter-genus"
    return cbcfg


def _compute_internal_breaks(piece_hsps: list[dict],
                             ref_genes: list[dict]) -> list[dict]:
    if len(piece_hsps) < 2:
        return []

    def _gene_info_at(pos: int) -> tuple[Optional[str], bool]:
        g = next((g for g in ref_genes if g["rel_start"] <= pos <= g["rel_end"]), None)
        return (g["locus_tag"], g.get("exception", False)) if g else (None, False)

    srt = sorted(piece_hsps, key=lambda h: min(int(h["qstart"]), int(h["qend"])))
    breaks: list[dict] = []
    for i in range(len(srt) - 1):
        prev, nxt = srt[i], srt[i + 1]
        q_hi_prev = max(int(prev["qstart"]), int(prev["qend"]))
        q_lo_next = min(int(nxt["qstart"]),  int(nxt["qend"]))
        s_hi_prev = max(int(prev["sstart"]), int(prev["send"]))
        s_lo_next = min(int(nxt["sstart"]),  int(nxt["send"]))
        midpoint = (q_hi_prev + q_lo_next) // 2

        ref_tag, is_exc = _gene_info_at(midpoint)
        if is_exc:
            continue  # breaks inside exception genes are not shown

        s_gap = s_lo_next - s_hi_prev
        breaks.append({
            "q_gap_lo":            q_hi_prev,
            "q_gap_hi":            q_lo_next,
            "q_gap_bp":            max(0, q_lo_next - q_hi_prev),
            "s_gap_bp":            s_gap,
            "subject_insertion_bp": max(0, s_gap - max(0, q_lo_next - q_hi_prev)),
            "ref_gene":            ref_tag,
        })
    return breaks


# ── Cluster blastn HSP analysis ──────────────────────────────────────────────────────

def _analyze_cluster_hsps(
    hsps: list[dict],
    locus_cfg: dict,
    expected_bp: Optional[int],
    ref_cluster_bp: Optional[int],
    contig_hint: Optional[str] = None,
    tol_pct: float = 15.0,
    min_coverage: float = 0.15,
    min_genes: int = 2,
    min_piece_identity: float = 70.0,
    single_gene_min_identity: float = 90.0,
    min_gene_coverage: float = 0.30,
    single_gene_min_coverage: float = 0.70,
    max_intra_piece_gap_bp: int = -1,
    single_gene_min_identity_lenient: Optional[float] = None,
    single_gene_min_coverage_lenient: Optional[float] = None,
    min_piece_identity_lenient: Optional[float] = None,
    min_gene_coverage_lenient: Optional[float] = None,
    domain_rescue_enabled: bool = True,
    domain_rescue_min_identity: float = 40.0,
    domain_rescue_min_gene_coverage: float = 0.15,
    domain_rescue_max_gap_bp: Optional[int] = None,
) -> dict:
    empty_result = {
        "status":          ABSENT,
        "contig":          None,
        "pieces":          [],
        "orphans":         [],
        "coverage":        0.0,
        "span":            0,
        "pos_min":         None,
        "pos_max":         None,
        "covered_genes":   [],
        "high_confidence": False,
        "orientation":     0,
        "piece_count":     0,
        "contig_count":    0,
        "internal_breaks": [],
        "assembly_hsps": [],
        "assembly_excluded_hsps": [],
    }

    if not hsps or ref_cluster_bp is None or ref_cluster_bp == 0:
        return empty_result

    ref_genes = _reference_gene_coords(locus_cfg)
    if not ref_genes:
        out = dict(empty_result)
        out["status"]  = UNKNOWN
        out["orphans"] = list(hsps)
        return out

    # Cluster-gene query envelope, in the same coordinates as ref_genes.
    _cl_q_lo = min(g["rel_start"] for g in ref_genes)
    _cl_q_hi = max(g["rel_end"]   for g in ref_genes)

    if max_intra_piece_gap_bp is None or max_intra_piece_gap_bp < 0:
        max_intra_piece_gap_bp = ref_cluster_bp

    exception_lts = {g["locus_tag"] for g in ref_genes if g.get("exception")}
    lenient_lts   = {g["locus_tag"] for g in ref_genes if g.get("lenient")}
    assembly_hsps: list[dict] = []
    assembly_excluded: list[dict] = []
    # Effective lenient thresholds (fall back to standard if None)
    eff_lenient_id = (min_piece_identity_lenient
                     if min_piece_identity_lenient is not None else min_piece_identity)
    eff_lenient_cov = (min_gene_coverage_lenient
                     if min_gene_coverage_lenient is not None else min_gene_coverage)
    for idx, h in enumerate(hsps):
        h["_hsp_idx"] = idx
        try:
            pid = float(h.get("pident", 0))
            qs = min(int(h["qstart"]), int(h["qend"]))
            qe = max(int(h["qstart"]), int(h["qend"]))
        except (KeyError, TypeError, ValueError):
            continue
        details = _gene_overlap_details(qs, qe, ref_genes)
        h["_gene_overlap_details"] = details
        # Standard eligibility: ≥ min_gene_coverage on non-exception gene
        eligible_genes = [
            d["locus_tag"] for d in details
            if d["coverage"] >= min_gene_coverage
            and d["locus_tag"] not in exception_lts
        ]
        # Lenient eligibility: the HSP covers a lenient gene at eff_lenient_cov.
        lenient_eligible = [
            d["locus_tag"] for d in details
            if d["locus_tag"] in lenient_lts
            and d["locus_tag"] not in exception_lts
            and d["coverage"] >= eff_lenient_cov
        ]
        # Standard identity first, then the lenient threshold for divergent genes.
        if pid >= min_piece_identity and eligible_genes:
            h["_assembly_genes"] = eligible_genes
            h["_used_for_assembly"] = True
            assembly_hsps.append(h)
        elif pid >= eff_lenient_id and lenient_eligible:
            # Lenient rescue: lower identity floor for divergent lenient genes.
            h["_assembly_genes"] = sorted(set(eligible_genes) | set(lenient_eligible))
            h["_used_for_assembly"] = True
            h["_lenient_pass"] = True
            assembly_hsps.append(h)
        elif pid >= min_piece_identity:
            h["_assembly_excluded_reason"] = "intergenic_q"
            h["_used_for_assembly"] = False
            assembly_excluded.append(h)

    # Pre-filter piece formation on the lower of the standard and lenient identity.
    piece_groups = _form_pieces(assembly_hsps, max_intra_piece_gap_bp,
                                min(min_piece_identity, eff_lenient_id))

    # 1b. Cap over-represented families (paralog HSP competition)
    piece_groups = _cap_family_hsps(piece_groups, locus_cfg, ref_genes)

    # 2 + 3. Per-piece annotation and acceptance
    accepted_entries: list[dict]           = []
    orphan_hsps:     list[dict]            = []
    for grp in piece_groups:
        _annotate_genes(grp, ref_genes, min_gene_coverage)
        ok, covered, tier = _accept_piece(
            grp, ref_genes, min_genes, single_gene_min_identity,
            min_gene_coverage, single_gene_min_coverage,
            single_gene_min_id_lenient=single_gene_min_identity_lenient,
            single_gene_min_coverage_lenient=single_gene_min_coverage_lenient,
            min_gene_coverage_lenient=min_gene_coverage_lenient)
        if ok:
            accepted_entries.append({"group": grp, "covered": covered, "tier": tier})
        else:
            orphan_hsps.extend(grp)

    low_identity_orphans = [
        h for h in hsps
        if float(h.get("pident", 0)) < min_piece_identity
        and not h.get("_used_for_assembly")
    ]

    if domain_rescue_enabled and accepted_entries:
        rescue_candidates = orphan_hsps + assembly_excluded + low_identity_orphans
        accepted_groups = [e["group"] for e in accepted_entries]
        rescued_groups = attach_domain_guided_weak_hsps(
            accepted_groups,
            rescue_candidates,
            ref_genes,
            min_identity=domain_rescue_min_identity,
            min_gene_coverage=domain_rescue_min_gene_coverage,
            max_gap_bp=(max_intra_piece_gap_bp
                        if domain_rescue_max_gap_bp is None
                        else domain_rescue_max_gap_bp),
        )
        for entry, grp in zip(accepted_entries, rescued_groups):
            entry["group"] = grp

    rescued_idx = {
        int(h.get("_hsp_idx", -1))
        for e in accepted_entries
        for h in e["group"]
        if h.get("_weak_hsp_rescue")
    }
    if rescued_idx:
        existing_assembly_idx = {int(h.get("_hsp_idx", -1)) for h in assembly_hsps}
        for e in accepted_entries:
            for h in e["group"]:
                h_idx = int(h.get("_hsp_idx", -1))
                if h_idx in rescued_idx and h_idx not in existing_assembly_idx:
                    assembly_hsps.append(h)
                    existing_assembly_idx.add(h_idx)
    orphan_hsps = [
        h for h in orphan_hsps + assembly_excluded + low_identity_orphans
        if int(h.get("_hsp_idx", -1)) not in rescued_idx
    ]

    accepted_pieces: list[dict] = []
    accepted_breaks: list[list[dict]] = []
    for entry in accepted_entries:
        grp = entry["group"]
        covered = sorted(set(entry["covered"]) |
                         {lt for h in grp for lt in h.get("_assembly_genes", [])})
        consol = _consolidate_piece(grp, q_clip=(_cl_q_lo, _cl_q_hi))
        consol["_genes"] = covered
        consol["_tier"] = entry["tier"]
        consol["_hsp_indices"] = [int(h.get("_hsp_idx", -1)) for h in grp]
        consol["_weak_hsp_rescue_count"] = sum(1 for h in grp if h.get("_weak_hsp_rescue"))
        accepted_pieces.append(consol)
        accepted_breaks.append(_compute_internal_breaks(grp, ref_genes))

    if not accepted_pieces:
        out = dict(empty_result)
        out["orphans"] = orphan_hsps
        out["assembly_hsps"] = assembly_hsps
        out["assembly_excluded_hsps"] = assembly_excluded
        return out

    order = sorted(range(len(accepted_pieces)),
                   key=lambda i: min(int(accepted_pieces[i]["qstart"]), int(accepted_pieces[i]["qend"])))
    pieces = [accepted_pieces[i] for i in order]
    breaks = [accepted_breaks[i] for i in order]

    # 4. Aggregate metrics using HSP UNION across ALL accepted pieces.
    global_intervals = []
    accepted_idx = {
        int(h_idx)
        for p in pieces
        for h_idx in (p.get("_hsp_indices") or [])
        if int(h_idx) >= 0
    }
    for h in hsps:
        if int(h.get("_hsp_idx", -1)) in accepted_idx:
            qs = max(min(int(h["qstart"]), int(h["qend"])), _cl_q_lo)
            qe = min(max(int(h["qstart"]), int(h["qend"])), _cl_q_hi)
            if qe >= qs:
                global_intervals.append((qs, qe))

    # Clip to the cluster-gene envelope and clamp coverage at 1.0.
    coverage = min(_union_length(global_intervals) / ref_cluster_bp, 1.0)

    piece_count  = len(pieces)
    contigs      = [p["sseqid"] for p in pieces]
    contig_count = len(set(contigs))

    covered_genes: list[str] = []
    seen = set()
    for p in pieces:
        for g in p.get("_genes", []):
            if g not in seen:
                covered_genes.append(g)
                seen.add(g)

    if contig_count == 1:
        positions: list[int] = []
        for p in pieces:
            s, e = int(p["sstart"]), int(p["send"])
            positions += [min(s, e), max(s, e)]
        pos_min, pos_max = min(positions), max(positions)
        span, primary_contig = pos_max - pos_min, pieces[0]["sseqid"]
    else:
        pos_min, pos_max, span, primary_contig = None, None, 0, "|".join(contigs)

    if coverage < min_coverage:
        out = dict(empty_result)
        out.update({"contig": primary_contig, "pieces": pieces, "orphans": orphan_hsps,
                    "coverage": coverage, "piece_count": piece_count,
                    "contig_count": contig_count,
                    "assembly_hsps": assembly_hsps,
                    "assembly_excluded_hsps": assembly_excluded})
        return out

    if piece_count > 1:
        # Multi-piece clusters start as FRAGMENTED; _verify_duplications may upgrade
        # them to DUPLICATED.
        status = FRAGMENTED
    else:
        # Single-piece cases
        status = INTACT if coverage >= 0.85 else (PARTIAL_DEL if coverage >= 0.30 else LARGE_DEL)

    return {
        "status": status, "contig": primary_contig, "pieces": pieces, "orphans": orphan_hsps,
        "coverage": coverage, "span": span, "pos_min": pos_min, "pos_max": pos_max,
        "covered_genes": covered_genes, "high_confidence": coverage >= 0.30,
        "orientation": _cluster_orientation(pieces), "piece_count": piece_count,
        # Nested list[list[dict]] — one inner list per piece, mirroring `pieces`. Do not flatten.
        "contig_count": contig_count, "internal_breaks": breaks,
        "assembly_hsps": assembly_hsps,
        "assembly_excluded_hsps": assembly_excluded,
    }


def _make_result(accession: str, status: str = ABSENT, contig: Optional[str] = None,
                 flank_bp: Optional[int] = None, cluster_bp: Optional[int] = None, notes: str = "") -> dict:
    return {
        "accession": accession, "status": status, "contig": contig, "flank_bp": flank_bp, "cluster_bp": cluster_bp,
        "coverage": None, "_cluster_pos_min": None, "_cluster_pos_max": None, "_cluster_orientation": 0,
        "_pieces": None, "_orphans": None, "_covered_genes": None,
        "_L2_lo": None, "_L2_hi": None, "_L1_lo": None, "_L1_hi": None,
        "_R1_lo": None, "_R1_hi": None, "_R2_lo": None, "_R2_hi": None,
        "L2_L1_bp": None, "R1_R2_bp": None, "order_flags": [], "round_results": [],
        "piece_count": 0, "contig_count": 0, "_internal_breaks": [], "notes": notes,
        "genes_found": None, "genes_expected": None, "pseudo_count": None,
    }


def _flank_layer(accession: str, hits: dict[str, list[dict]], anchors: list[dict],
                 min_id: float, anchor_min_id: float) -> dict:
    def get_a(r): return next((a for a in anchors if a["role"] == r), None)
    L2_a, L1_a, CL_a, CR_a, R1_a, R2_a = get_a("flank_L2"), get_a("flank_L"), get_a("cluster_L"), get_a("cluster_R"), get_a("flank_R"), get_a("flank_R2")

    def find_h(a, use_a_min=True):
        if a is None: return None
        return _best_hit(hits.get(a["locus_tag"], []), anchor_min_id if use_a_min else min_id)

    def pos_i(h):
        if h is None: return {"found": False, "contig": None, "lo": None, "hi": None}
        c, lo, hi = _hit_span(h)
        return {"found": True, "contig": c, "lo": lo, "hi": hi}

    L2, L1, CL, CR, R1, R2 = pos_i(find_h(L2_a)), pos_i(find_h(L1_a)), pos_i(find_h(CL_a, False)), pos_i(find_h(CR_a, False)), pos_i(find_h(R1_a)), pos_i(find_h(R2_a))
    core_c = L1["contig"] if L1["found"] else (R1["contig"] if R1["found"] else None)
    L2_L1_s, L1_R1_s, R1_R2_s = L2["found"] and L1["found"] and L2["contig"]==L1["contig"], L1["found"] and R1["found"] and L1["contig"]==R1["contig"], R1["found"] and R2["found"] and R1["contig"]==R2["contig"]

    def gap(h1, h2):
        if h1["lo"] is None or h2["lo"] is None: return None
        return h2["lo"] - h1["hi"] if h1["lo"] < h2["lo"] else h1["lo"] - h2["hi"]

    CL_in = CR_in = CL_oth = CR_oth = False
    L1_CL_b = CR_R1_b = None
    if L1_R1_s:
        wl, wh = min(L1["lo"], R1["lo"]), max(L1["hi"], R1["hi"])
        if CL["found"]:
            if CL["contig"] == core_c:
                CL_in = wl <= CL["lo"] <= wh
                if CL_in: L1_CL_b = gap(L1, CL)
            else: CL_oth = True
        if CR["found"]:
            if CR["contig"] == core_c:
                CR_in = wl <= CR["lo"] <= wh
                if CR_in: CR_R1_b = gap(CR, R1)
            else: CR_oth = True

    flags = []
    if L2_L1_s and L2["lo"] > L1["lo"]: flags.append("L2>L1")
    if L1_R1_s and L1["lo"] > R1["lo"]: flags.append("L1>R1")
    if R1_R2_s and R1["lo"] > R2["lo"]: flags.append("R1>R2")

    return {
        "L2":L2,"L1":L1,"CL":CL,"CR":CR,"R1":R1,"R2":R2, "core_contig":core_c,
        "L2_L1_same":L2_L1_s,"L1_R1_same":L1_R1_s,"R1_R2_same":R1_R2_s,
        "L2_L1_bp":gap(L2,L1) if L2_L1_s else None, "L1_R1_bp":gap(L1,R1) if L1_R1_s else None, "R1_R2_bp":gap(R1,R2) if R1_R2_s else None,
        "L1_CL_bp":L1_CL_b, "CR_R1_bp":CR_R1_b, "CL_in_window":CL_in, "CR_in_window":CR_in, "CL_other_contig":CL_oth, "CR_other_contig":CR_oth, "order_flags":flags
    }


def _combine(accession, flank, blast, anchors, exp_bp):
    res = _make_result(accession)
    res.update({
        "status": blast["status"], "cluster_bp": blast.get("span"), "order_flags": flank["order_flags"],
        "contig": blast.get("contig") or flank.get("core_contig"), "coverage": blast.get("coverage"),
        "_cluster_pos_min": blast.get("pos_min"), "_cluster_pos_max": blast.get("pos_max"),
        "_cluster_orientation": blast.get("orientation", 0), "_pieces": blast.get("pieces"),
        "_orphans": blast.get("orphans"), "_covered_genes": blast.get("covered_genes"),
        "piece_count": blast.get("piece_count", 0), "contig_count": blast.get("contig_count", 0),
        "_internal_breaks": blast.get("internal_breaks") or [],
        "_assembly_hsps": blast.get("assembly_hsps") or [],
        "_assembly_excluded_hsps": blast.get("assembly_excluded_hsps") or [],
    })
    for k, info in (("_L2", flank["L2"]), ("_L1", flank["L1"]), ("_R1", flank["R1"]), ("_R2", flank["R2"])):
        res[k + "_lo"], res[k + "_hi"] = info.get("lo"), info.get("hi")

    cov, span_s, exp_s = blast.get("coverage") or 0.0, f"{blast['span']:,}" if blast.get("span") else "?", f"{exp_bp.get('cluster', 0):,}"
    parts = [f"coverage={cov:.0%} cluster_span={span_s} bp (expected {exp_s})"]
    if not flank["L1"]["found"] and not flank["R1"]["found"]: parts.append("no flanks (blast-only)")
    res["notes"] = "; ".join(parts + (flank["order_flags"] or []))
    return res


def _assess_genome(accession, hits, locus_cfg, settings, cluster_hsps=None,
                   genome_meta=None):
    auto, anchors, exp_bp = locus_cfg.get("_auto", {}), locus_cfg.get("_auto", {}).get("anchors", []), locus_cfg.get("_auto", {}).get("expected_bp", {})
    bcfg, rcfg = dict(settings.get("blast", {})), dict(settings.get("ruler", {}))
    # merge with locus overrides
    bcfg.update(locus_cfg.get("blast", {})); rcfg.update(locus_cfg.get("ruler", {}))
    meta = genome_meta or {}
    cbcfg = _resolve_cluster_blast_cfg(settings, locus_cfg, meta.get("species", ""))

    if cluster_hsps is None: return _make_result(accession, UNKNOWN, notes="cluster_hsps missing")

    blast = _analyze_cluster_hsps(cluster_hsps, locus_cfg, exp_bp.get("cluster"), auto.get("cluster_ref_bp"),
                                 tol_pct=float(rcfg.get("bp_tolerance_pct", 15)),
                                 min_genes=int(cbcfg.get("min_genes", 2)),
                                 min_piece_identity=float(cbcfg.get("min_piece_identity", 70)),
                                 single_gene_min_identity=float(cbcfg.get("single_gene_min_identity", 90)),
                                 min_gene_coverage=float(cbcfg.get("min_gene_coverage", 0.30)),
                                 single_gene_min_coverage=float(cbcfg.get("single_gene_min_coverage", 0.70)),
                                 single_gene_min_identity_lenient=float(cbcfg.get("single_gene_min_identity_lenient", 40)),
                                 single_gene_min_coverage_lenient=float(cbcfg.get("single_gene_min_coverage_lenient", 0.40)),
                                 min_piece_identity_lenient=float(cbcfg.get("min_piece_identity_lenient", 40)),
                                 min_gene_coverage_lenient=float(cbcfg.get("min_gene_coverage_lenient", 0.40)),
                                 domain_rescue_enabled=bool(cbcfg.get("domain_rescue_enabled", True)),
                                 domain_rescue_min_identity=float(cbcfg.get("domain_rescue_min_identity", 40)),
                                 domain_rescue_min_gene_coverage=float(cbcfg.get("domain_rescue_min_gene_coverage", 0.15)),
                                 domain_rescue_max_gap_bp=(
                                     None if cbcfg.get("domain_rescue_max_gap_bp") in (None, "", -1)
                                     else int(cbcfg.get("domain_rescue_max_gap_bp"))
                                 ),
                                 min_coverage=float(cbcfg.get("min_coverage", 0.15)))

    has_L1 = any(a["role"] == "flank_L" for a in anchors)
    has_R1 = any(a["role"] == "flank_R" for a in anchors)

    if not has_L1 and not has_R1:
        res = _make_result(accession, blast["status"], blast["contig"], cluster_bp=blast["span"])
        res.update({"coverage": blast["coverage"], "_pieces": blast["pieces"], "_orphans": blast["orphans"],
                    "_covered_genes": blast["covered_genes"], "piece_count": blast["piece_count"],
                    "contig_count": blast["contig_count"], "_internal_breaks": blast["internal_breaks"],
                    "_assembly_hsps": blast.get("assembly_hsps") or [],
                    "_assembly_excluded_hsps": blast.get("assembly_excluded_hsps") or []})
        return res

    flank = _flank_layer(accession, hits, anchors, float(bcfg.get("min_identity", 50)), float(rcfg.get("anchor_min_identity", 40)))
    return _combine(accession, flank, blast, anchors, exp_bp)


def run_ruler(locus_cfg, settings, db_map, blast_hits, work_dir, cluster_hits=None, cpu=4, force=False, no_split_search=False, target_name: str = "", locus_dir: "Path | None" = None):
    locus_id = locus_cfg["locus_id"]
    tgt_name   = target_name or locus_cfg.get("reference", {}).get("target", "")
    tgt_cfg_ws = next((t for t in settings.get("targets", []) if t["name"] == tgt_name), None)
    genome_meta_ws = _load_ruler_genome_meta(
        tgt_cfg_ws.get("db", "") if tgt_cfg_ws else "",
        list(db_map.keys()),
    )

    results = {
        acc: _assess_genome(
            acc, blast_hits.get(acc, {}), locus_cfg, settings,
            cluster_hits.get(acc, []) if cluster_hits is not None else None,
            genome_meta=genome_meta_ws.get(acc, {}),
        )
        for acc in db_map
    }
    _verify_duplications(results, locus_cfg, settings, target_name=target_name)

    _locus_dir = locus_dir if locus_dir is not None else work_dir / locus_id
    _locus_dir.mkdir(parents=True, exist_ok=True)
    with open(_locus_dir / "ruler_results.json", "w") as f: json.dump(results, f, indent=2, default=str)
    write_summary(results, locus_cfg, _locus_dir / "genome_summary.csv",
                  genome_meta=genome_meta_ws)
    return results

def write_summary(results, locus_cfg, out_path, genome_meta=None):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    auto   = locus_cfg.get("_auto", {})
    exp_bp = auto.get("expected_bp", {})
    _meta  = genome_meta or {}
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_COLS)
        w.writeheader()
        for acc in sorted(results):
            res  = results[acc]
            fgi  = res.get("_flank_gene_info") or {}
            fs   = res.get("_flank_states") or {}
            meta = _meta.get(acc, {})

            flags = res.get("order_flags") or []

            w.writerow({
                "genome_acc":          acc,
                "species":             meta.get("species", ""),
                "strain":              meta.get("strain",  ""),
                "assembly_level":      meta.get("assembly_level", ""),
                "locus_id":            locus_cfg["locus_id"],
                "status":              res.get("status", ""),
                "synteny_class":       res.get("_synteny_class", ""),
                "synteny_case":        res.get("_synteny_case", ""),
                "L_state":             fs.get("left", ""),
                "R_state":             fs.get("right", ""),
                "L_fingerprint":       res.get("_L_fingerprint", ""),
                "R_fingerprint":       res.get("_R_fingerprint", ""),
                "L2_gene":             fgi.get("L2", ""),
                "L1_gene":             fgi.get("L1", ""),
                "R1_gene":             fgi.get("R1", ""),
                "R2_gene":             fgi.get("R2", ""),
                "inversion":           _flank_relative_inversion(res),
                "orientation":         res.get("_cluster_orientation", ""),
                "coverage":            f"{res.get('coverage') or 0:.1%}",
                "L2_L1_bp":            res.get("L2_L1_bp", ""),
                "flank_bp":            res.get("flank_bp", ""),
                "R1_R2_bp":            res.get("R1_R2_bp", ""),
                "cluster_bp":          res.get("cluster_bp", ""),
                "expected_L2_L1_bp":   exp_bp.get("L2_L1", ""),
                "expected_flank_bp":   exp_bp.get("flank", ""),
                "expected_R1_R2_bp":   exp_bp.get("R1_R2", ""),
                "expected_cluster_bp": exp_bp.get("cluster", ""),
                "order_flags":         "; ".join(flags) if flags else "",
                "genes_found":         res.get("genes_found", ""),
                "genes_expected":      res.get("genes_expected", ""),
                "pseudo_count":        res.get("pseudo_count", ""),
                "linked_locus":        res.get("linked_locus", ""),
                "contig":              res.get("contig", ""),
                "notes":               res.get("notes", ""),
                # v1.0 fragmentation annotation
                "fragmentation_type":  res.get("_fragmentation_type", ""),
                "contig_edge_support": res.get("_contig_edge_support", ""),
                "q_pair_max_overlap":  (
                    f"{res['_q_pair_max_overlap']:.2f}"
                    if "_q_pair_max_overlap" in res else ""
                ),
            })

def load_ruler_results(locus_dir: "Path"):
    """Load ruler_results.json from *locus_dir*."""
    from pathlib import Path as _Path
    with open(_Path(locus_dir) / "ruler_results.json") as f:
        return json.load(f)
