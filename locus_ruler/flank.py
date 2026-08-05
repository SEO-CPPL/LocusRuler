#!/usr/bin/env python3
"""Flank validation utilities."""

import subprocess
from pathlib import Path
from typing import Optional


def get_contig_lengths(
    blastdb_path: str,
    blastdbcmd: str = "blastdbcmd",
) -> dict[str, int]:
    """Return ``{contig_id: length}`` for every sequence in a BLAST DB."""
    try:
        result = subprocess.run(
            [blastdbcmd, "-db", blastdb_path, "-entry", "all",
             "-outfmt", "%t\t%l"],
            capture_output=True, text=True, check=True, timeout=60,
        )
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return {}
    out: dict[str, int] = {}
    for line in result.stdout.splitlines():
        if "\t" not in line:
            continue
        title, length_s = line.rsplit("\t", 1)
        first_word = title.strip().split()[0] if title.strip() else ""
        if not first_word:
            continue
        try:
            out[first_word] = int(length_s)
        except ValueError:
            continue
    return out


def gff_gene_at(gff_genes: list[dict], lo: int, hi: int) -> Optional[dict]:
    """Return the GFF gene whose interval overlaps [lo, hi]."""
    if not gff_genes or lo is None or hi is None:
        return None
    best_overlap = 0
    best_gene = None
    for g in gff_genes:
        if g["end"] < lo or g["start"] > hi:
            continue
        ov = min(g["end"], hi) - max(g["start"], lo)
        if ov > best_overlap:
            best_overlap = ov
            best_gene = g
    return best_gene


# ── Split CDS detection helpers ──────────────────────────────────────────────────────
def _adjacent_same_strand_cds_pair(
    gff_genes: list[dict],
    best_gene: dict,
) -> tuple[Optional[dict], Optional[dict]]:
    """Return the immediately-adjacent same-strand CDS pair (prev, next)."""
    if not gff_genes or not best_gene:
        return (None, None)
    contig = best_gene.get("contig")
    strand = best_gene.get("strand")
    same_contig = [g for g in gff_genes if g.get("contig") == contig]
    same_contig.sort(key=lambda g: (int(g["start"]), int(g["end"])))

    # Locate best_gene by exact (start, end, strand)
    bs, be = int(best_gene["start"]), int(best_gene["end"])
    try:
        i = next(
            idx for idx, g in enumerate(same_contig)
            if int(g["start"]) == bs and int(g["end"]) == be
            and g.get("strand") == strand
        )
    except StopIteration:
        return (None, None)

    def _is_cds(g):
        ft = g.get("feature_type")
        return ft in ("CDS", "pseudogene")

    prev_cds = None
    for j in range(i - 1, -1, -1):
        g = same_contig[j]
        if not _is_cds(g):
            continue
        if g.get("strand") == strand:
            prev_cds = g
        break  # first CDS we encounter (regardless of strand match) is "next door"

    next_cds = None
    for j in range(i + 1, len(same_contig)):
        g = same_contig[j]
        if not _is_cds(g):
            continue
        if g.get("strand") == strand:
            next_cds = g
        break

    return (prev_cds, next_cds)


def adjacent_pseudo_partner(
    gff_genes: list[dict],
    best_gene: Optional[dict],
) -> Optional[dict]:
    """Option B: return adjacent same-strand CDS marked GFF pseudo, if any."""
    if not best_gene:
        return None
    prev_cds, next_cds = _adjacent_same_strand_cds_pair(gff_genes, best_gene)
    for partner in (prev_cds, next_cds):
        if partner and partner.get("is_pseudo"):
            return partner
    return None


def adjacent_split_cds_pair(
    anchor_hits: list[dict],
    gff_genes: list[dict],
    best_hit: Optional[dict],
    best_gene: Optional[dict],
    min_partner_pid: float = 25.0,
    min_partner_aa: int = 80,
    max_ref_overlap_frac: float = 0.30,
) -> Optional[dict]:
    """Option A: return adjacent same-strand CDS that hits the SAME reference
    with reference-coverage range complementary (N/C) to the best hit.

    Distinguishes split CDS (two HSPs span complementary parts of the
    reference protein) from tandem paralog (two HSPs over the same reference
    region).  Returns the partner CDS dict on detection, ``None`` otherwise.

    Parameters
    ----------
    anchor_hits : list[dict]
        All tblastn HSP rows for one reference anchor against one genome,
        pre-filtered to pid ≥ 25%.  Each row carries `qstart, qend, sstart,
        send, sseqid, pident, length`.
    gff_genes : list[dict]
        Same target genome's GFF gene records (per-contig list).
    best_hit : dict
        The selected best HSP (must carry `lo, hi, contig` plus the original
        `qstart, qend` from blast.py — see caller-provided wrapper).
    best_gene : dict
        The GFF gene that overlaps the best hit (`gff_gene_at(...)` result).
    min_partner_pid : float
        Minimum tblastn pid for a partner HSP to count.
    min_partner_aa : int
        Minimum protein length of the candidate partner CDS — filters out
        tiny ORF artifacts.  Protein length is read from the GFF gene's
        `(end - start + 1) / 3` if `protein_len` is unavailable.
    max_ref_overlap_frac : float
        If (overlap_ref / union_ref) < this, the two HSPs are considered
        complementary (= split).  Above this threshold they overlap in the
        reference and we treat them as tandem paralog, NOT split.
    """
    if not anchor_hits or not best_hit or not best_gene:
        return None
    best_qs = best_hit.get("qstart") or best_hit.get("_qstart")
    best_qe = best_hit.get("qend")   or best_hit.get("_qend")
    if best_qs is None or best_qe is None:
        return None
    try:
        best_qs = int(best_qs)
        best_qe = int(best_qe)
    except (TypeError, ValueError):
        return None
    best_ref_lo = min(best_qs, best_qe)
    best_ref_hi = max(best_qs, best_qe)

    prev_cds, next_cds = _adjacent_same_strand_cds_pair(gff_genes, best_gene)
    partners = [p for p in (prev_cds, next_cds) if p]
    if not partners:
        return None

    def _hsp_on_cds(h: dict, cds: dict) -> bool:
        try:
            hs_lo = min(int(h["sstart"]), int(h["send"]))
            hs_hi = max(int(h["sstart"]), int(h["send"]))
        except (KeyError, TypeError, ValueError):
            return False
        if h.get("sseqid") != cds.get("contig"):
            return False
        c_lo = int(cds["start"])
        c_hi = int(cds["end"])
        return hs_hi >= c_lo and hs_lo <= c_hi

    for partner in partners:
        # Length filter (CDS bp → approx protein aa)
        partner_aa = max(0, (int(partner["end"]) - int(partner["start"]) + 1) // 3)
        if partner_aa < min_partner_aa:
            continue
        partner_hsps = [
            h for h in anchor_hits
            if float(h.get("pident") or 0) >= min_partner_pid
            and _hsp_on_cds(h, partner)
        ]
        if not partner_hsps:
            continue
        ph = max(partner_hsps, key=lambda h: int(h.get("length") or 0))
        try:
            p_qs = int(ph.get("qstart") or 0)
            p_qe = int(ph.get("qend")   or 0)
        except (TypeError, ValueError):
            continue
        if not (p_qs and p_qe):
            continue
        p_lo = min(p_qs, p_qe)
        p_hi = max(p_qs, p_qe)
        ov = max(0, min(best_ref_hi, p_hi) - max(best_ref_lo, p_lo))
        un = max(best_ref_hi, p_hi) - min(best_ref_lo, p_lo)
        if un <= 0:
            continue
        if (ov / un) < max_ref_overlap_frac:
            return partner  # complementary N/C ranges → split CDS pair
    return None
