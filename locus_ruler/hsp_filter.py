#!/usr/bin/env python3
"""HSP-level family count cap."""

from typing import Optional


def cap_family_hsps(
    piece_groups: list[list[dict]],
    locus_cfg: dict,
    ref_genes: list[dict],
) -> list[list[dict]]:
    """Cap each anchor family's HSP count within each piece to the number
    of times that family appears in the reference cluster.

    Family membership is determined by qstart/qend overlap with reference
    gene coordinates (same logic as _covered_genes_union in ruler.py).
    Excess HSPs (lowest bitscore) are discarded.  Exception-flagged anchors
    are never capped.

    Returns a new list of piece groups; empty groups after filtering are
    removed.
    """
    anchors = locus_cfg.get("_auto", {}).get("anchors", [])

    cluster_roles = {"cluster_L", "cluster_R", "inner"}
    ref_family_count: dict[str, int] = {}
    exception_lts: set[str] = set()

    # Map locus_tag → family; build ref_family_count per family
    lt_to_family: dict[str, str] = {}
    for a in anchors:
        if a.get("role") not in cluster_roles:
            continue
        lt  = a.get("locus_tag", "")
        fam = a.get("family") or lt
        if not fam or not lt:
            continue
        lt_to_family[lt] = fam
        if a.get("exception"):
            exception_lts.add(lt)
        else:
            ref_family_count[fam] = ref_family_count.get(fam, 0) + 1

    if not ref_family_count or not ref_genes:
        return piece_groups

    # Build per-gene family lookup (locus_tag → family) and exception set
    gene_family: dict[str, str] = {}
    exception_fams: set[str] = set()
    for g in ref_genes:
        lt  = g["locus_tag"]
        fam = lt_to_family.get(lt, lt)
        gene_family[lt] = fam
        if g.get("exception"):
            exception_fams.add(fam)

    result: list[list[dict]] = []
    for grp in piece_groups:
        filtered = _cap_group(grp, ref_genes, gene_family,
                              ref_family_count, exception_fams)
        if filtered:
            result.append(filtered)
    return result


def _hsp_covered_families(
    h: dict,
    ref_genes: list[dict],
    gene_family: dict[str, str],
    min_gene_coverage: float = 0.30,
) -> list[str]:
    """Return the list of families whose reference gene this HSP covers
    (≥ min_gene_coverage of that gene's length)."""
    qs = min(int(h.get("qstart", 0)), int(h.get("qend", 0)))
    qe = max(int(h.get("qstart", 0)), int(h.get("qend", 0)))
    families: list[str] = []
    for g in ref_genes:
        if g["rel_end"] < qs or g["rel_start"] > qe:
            continue
        ov = max(0, min(g["rel_end"], qe) - max(g["rel_start"], qs) + 1)
        gene_len = max(1, g["rel_end"] - g["rel_start"] + 1)
        if ov / gene_len >= min_gene_coverage:
            fam = gene_family.get(g["locus_tag"])
            if fam:
                families.append(fam)
    return families


def _cap_group(
    hsps: list[dict],
    ref_genes: list[dict],
    gene_family: dict[str, str],
    ref_family_count: dict[str, int],
    exception_fams: set[str],
) -> list[dict]:
    """Filter one piece group so each family appears at most ref_count times."""
    # Assign each HSP to its covered families
    hsp_families: list[list[str]] = []
    for h in hsps:
        fams = _hsp_covered_families(h, ref_genes, gene_family)
        # Only consider families that are capped and not exceptions
        capped_fams = [f for f in fams
                       if f in ref_family_count and f not in exception_fams]
        hsp_families.append(capped_fams)

    # Group HSP indices by their primary capped family.
    family_indices: dict[str, list[int]] = {}
    uncapped: list[int] = []

    for i, fams in enumerate(hsp_families):
        if fams:
            for fam in fams:
                family_indices.setdefault(fam, []).append(i)
        else:
            uncapped.append(i)

    keep: set[int] = set(uncapped)

    for fam, idxs in family_indices.items():
        cap = ref_family_count[fam]
        # Deduplicate (same HSP may appear under multiple families)
        unique_idxs = list(dict.fromkeys(idxs))
        if len(unique_idxs) <= cap:
            keep.update(unique_idxs)
            continue
        sorted_idxs = sorted(
            unique_idxs,
            key=lambda i: float(hsps[i].get("bitscore", 0)),
            reverse=True,
        )
        keep.update(sorted_idxs[:cap])
        n_dropped = len(unique_idxs) - cap
        print(f"[hsp_filter] capped family '{fam}': "
              f"kept {cap}, dropped {n_dropped} paralog HSP(s)")

    return [hsps[i] for i in sorted(keep)]
