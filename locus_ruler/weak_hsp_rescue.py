#!/usr/bin/env python3
"""Domain-configured weak HSP attachment."""

from typing import Optional


def _strand(h: dict) -> str:
    return "+" if int(h["send"]) > int(h["sstart"]) else "-"


def _srange(h: dict) -> tuple[int, int]:
    return min(int(h["sstart"]), int(h["send"])), max(int(h["sstart"]), int(h["send"]))


def _span_gap(a_lo: int, a_hi: int, b_lo: int, b_hi: int) -> int:
    if a_hi < b_lo:
        return b_lo - a_hi
    if b_hi < a_lo:
        return a_lo - b_hi
    return 0


def _group_span(group: list[dict]) -> tuple[str, str, int, int]:
    contig = group[0]["sseqid"]
    strand = _strand(group[0])
    lows: list[int] = []
    highs: list[int] = []
    for h in group:
        lo, hi = _srange(h)
        lows.append(lo)
        highs.append(hi)
    return contig, strand, min(lows), max(highs)


def _domain_enabled_lts(ref_genes: list[dict]) -> set[str]:
    out: set[str] = set()
    for g in ref_genes:
        if g.get("exception"):
            continue
        if g.get("pfam"):
            out.add(g["locus_tag"])
    return out


def _hsp_domain_gene_hits(
    h: dict,
    ref_genes: list[dict],
    enabled_lts: set[str],
    min_gene_coverage: float,
) -> list[dict]:
    qs = min(int(h["qstart"]), int(h["qend"]))
    qe = max(int(h["qstart"]), int(h["qend"]))
    hits: list[dict] = []
    for g in ref_genes:
        lt = g["locus_tag"]
        if lt not in enabled_lts:
            continue
        if g["rel_end"] < qs or g["rel_start"] > qe:
            continue
        ov = max(0, min(g["rel_end"], qe) - max(g["rel_start"], qs) + 1)
        gene_len = max(1, g["rel_end"] - g["rel_start"] + 1)
        cov = ov / gene_len
        if cov >= min_gene_coverage:
            hits.append({
                "locus_tag": lt,
                "family": g.get("family") or lt,
                "coverage": cov,
            })
    return hits


def attach_domain_guided_weak_hsps(
    accepted_groups: list[list[dict]],
    candidate_hsps: list[dict],
    ref_genes: list[dict],
    *,
    min_identity: float = 40.0,
    min_gene_coverage: float = 0.15,
    max_gap_bp: Optional[int] = None,
) -> list[list[dict]]:
    """Attach weak HSPs to already accepted groups."""
    if not accepted_groups or not candidate_hsps:
        return accepted_groups

    enabled_lts = _domain_enabled_lts(ref_genes)
    if not enabled_lts:
        return accepted_groups

    group_spans = [_group_span(g) for g in accepted_groups]
    rescued_by_group: dict[int, list[dict]] = {}
    seen_candidates: set[int] = set()

    for h in candidate_hsps:
        try:
            h_idx = int(h.get("_hsp_idx", -1))
            pid = float(h.get("pident", 0))
        except (TypeError, ValueError):
            continue
        if h_idx in seen_candidates or pid < min_identity:
            continue

        domain_hits = _hsp_domain_gene_hits(
            h, ref_genes, enabled_lts, min_gene_coverage)
        if not domain_hits:
            continue

        h_contig = h.get("sseqid", "")
        h_strand = _strand(h)
        h_lo, h_hi = _srange(h)
        best_group: Optional[int] = None
        best_gap: Optional[int] = None
        for i, (g_contig, g_strand, g_lo, g_hi) in enumerate(group_spans):
            if h_contig != g_contig or h_strand != g_strand:
                continue
            gap = _span_gap(h_lo, h_hi, g_lo, g_hi)
            if max_gap_bp is not None and gap > max_gap_bp:
                continue
            if best_gap is None or gap < best_gap:
                best_group = i
                best_gap = gap

        if best_group is None:
            continue

        h["_used_for_assembly"] = True
        h["_weak_hsp_rescue"] = True
        h["_assembly_excluded_reason"] = ""
        h["_rescue_reason"] = "domain_configured_weak_hsp"
        h["_rescue_gap_bp"] = best_gap
        h["_rescue_domain_hits"] = domain_hits
        h["_assembly_genes"] = sorted({d["locus_tag"] for d in domain_hits})
        rescued_by_group.setdefault(best_group, []).append(h)
        seen_candidates.add(h_idx)

    if not rescued_by_group:
        return accepted_groups

    out: list[list[dict]] = []
    for i, grp in enumerate(accepted_groups):
        additions = rescued_by_group.get(i, [])
        if additions:
            existing = {int(h.get("_hsp_idx", -1)) for h in grp}
            merged = list(grp)
            for h in additions:
                if int(h.get("_hsp_idx", -1)) not in existing:
                    merged.append(h)
            out.append(sorted(merged, key=lambda h: int(h.get("_hsp_idx", 10**9))))
        else:
            out.append(grp)
    return out
