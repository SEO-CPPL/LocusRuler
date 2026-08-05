#!/usr/bin/env python3
"""Per-gene state classification."""

from typing import Optional


def inner_anchor_tblastn_metrics(
    anchor_hits: list[dict],
    piece_layouts: list[dict],
    min_identity_floor: float = 30.0,
) -> dict:
    """Best tblastn hit metrics for one reference inner anchor at the cluster locus."""
    blank = {"cov": 0.0, "pid": 0.0, "hit": None}
    if not anchor_hits or not piece_layouts:
        return blank

    best = None
    best_score = 0.0

    def _pl_intervals(pl: dict) -> list[tuple[int, int]]:
        """Return the list of (lo, hi) subject intervals for a piece layout."""
        ivs = pl.get("subject_intervals")
        if ivs:
            return [(int(lo), int(hi)) for lo, hi in ivs]
        return [(int(pl["p_lo"]), int(pl["p_hi"]))]

    for h in anchor_hits:
        try:
            pid    = float(h.get("pident", 0))
            length = int(h.get("length", 0))
            qlen   = int(h.get("qlen", 0))
        except (TypeError, ValueError):
            continue
        if pid < min_identity_floor:
            continue
        h_contig = h.get("sseqid")
        try:
            s = min(int(h["sstart"]), int(h["send"]))
            e = max(int(h["sstart"]), int(h["send"]))
        except (KeyError, TypeError, ValueError):
            continue
        for pl in piece_layouts:
            if pl["contig"] != h_contig:
                continue
            # Hit qualifies if it overlaps any of the piece's subject intervals.
            if not any(e >= lo and s <= hi for lo, hi in _pl_intervals(pl)):
                continue
            score = length * pid
            if score > best_score:
                best_score = score
                # Preserve raw query coords for downstream split-CDS detection.
                try:
                    qs = int(h.get("qstart") or 0)
                    qe = int(h.get("qend") or 0)
                except (TypeError, ValueError):
                    qs = qe = 0
                best = {"pident": pid, "length": length, "qlen": qlen,
                        "lo": s, "hi": e, "contig": h_contig,
                        "qstart": qs, "qend": qe}
            break

    if best is None:
        return blank
    cov = best["length"] / max(1, best["qlen"])
    return {"cov": min(1.0, cov), "pid": best["pident"], "hit": best}


def cluster_dna_present(
    cluster_hsps_for_acc: list[dict],
    ref_anchor_q_range: tuple[int, int],
    piece_layouts: Optional[list[dict]] = None,
    threshold: float = 0.30,
) -> bool:
    return cluster_dna_metrics(
        cluster_hsps_for_acc, ref_anchor_q_range,
        piece_layouts=piece_layouts, threshold=threshold,
    )["present"]


def cluster_dna_metrics(
    cluster_hsps_for_acc: list[dict],
    ref_anchor_q_range: tuple[int, int],
    piece_layouts: Optional[list[dict]] = None,
    threshold: float = 0.30,
) -> dict:
    """Return cluster blastn DNA coverage for one anchor q-range."""
    qs, qe = ref_anchor_q_range
    qs, qe = min(qs, qe), max(qs, qe)
    gene_len = max(1, qe - qs + 1)

    def _in_a_piece(h_contig: str, s_lo: int, s_hi: int) -> bool:
        if piece_layouts is None:
            return True
        for pl in piece_layouts:
            if pl["contig"] != h_contig:
                continue
            ivs = pl.get("subject_intervals")
            intervals = ([(int(lo), int(hi)) for lo, hi in ivs] if ivs
                         else [(int(pl["p_lo"]), int(pl["p_hi"]))])
            for lo, hi in intervals:
                if s_hi >= lo and s_lo <= hi:
                    return True
        return False

    intervals: list[tuple[int, int]] = []
    for h in cluster_hsps_for_acc or []:
        try:
            hqs = min(int(h["qstart"]), int(h["qend"]))
            hqe = max(int(h["qstart"]), int(h["qend"]))
            hs_lo = min(int(h["sstart"]), int(h["send"]))
            hs_hi = max(int(h["sstart"]), int(h["send"]))
        except (KeyError, TypeError, ValueError):
            continue
        if hqe < qs or hqs > qe:
            continue
        if not _in_a_piece(h.get("sseqid", ""), hs_lo, hs_hi):
            continue
        intervals.append((max(qs, hqs), min(qe, hqe)))
    if not intervals:
        return {"present": False, "cov": 0.0, "intervals": []}
    intervals.sort()
    merged: list[list[int]] = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    total = sum(e - s + 1 for s, e in merged)
    cov = total / gene_len
    return {
        "present": cov >= threshold,
        "cov": cov,
        "intervals": [(int(s), int(e)) for s, e in merged],
    }


def classify_gene_state(
    tblastn_cov: float,
    tblastn_pid: float,
    cluster_dna_present: bool,
    gff_pseudo: bool,
    intact_coverage: float = 0.70,
    intact_identity: float = 50.0,
    min_coverage: float    = 0.30,
    min_identity: float    = 30.0,
    lenient: bool          = False,
) -> str:
    """tblastn-primary gene state classifier.

    Decision rules (first match wins):
      1. gff_pseudo                                  -> PSEUDOGENE
      2. cov >= intact_cov AND pid >= intact_pid      -> INTACT
      3. cov >= intact_cov AND pid >= min_pid         -> DIVERGENT
      4. min_cov <= cov < intact_cov AND pid >= min   -> DIVERGENT
      5. cov < min_cov AND pid >= intact_pid AND dna  -> DIVERGENT
      6. no meaningful protein/domain signal          -> ABSENT
      7. else                                         -> ABSENT

    ``intact_identity`` is dataset-aware: pass the intra- or inter-genus value
    as appropriate (see content.run_content for the selection logic).
    """
    if gff_pseudo:
        return "PSEUDOGENE"

    if tblastn_cov >= intact_coverage and tblastn_pid >= intact_identity:
        return "INTACT"

    if tblastn_cov >= intact_coverage and tblastn_pid >= min_identity:
        return "DIVERGENT"

    if min_coverage <= tblastn_cov < intact_coverage and tblastn_pid >= min_identity:
        return "DIVERGENT"

    if tblastn_cov < min_coverage and tblastn_pid >= intact_identity and cluster_dna_present:
        return "DIVERGENT"

    return "ABSENT"
