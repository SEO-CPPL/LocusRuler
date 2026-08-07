#!/usr/bin/env python3
"""Cluster heatmap rendering."""

_DESCRIPTION = """\
Strain by reference-gene matrix of per-piece coverage, internal breaks, and
multi-piece overlap, read from ruler_results.json.

Cell fill follows the synteny palette:
    both flanks canonical    navy        #1f77b4
    one side shared group    sky-blue    per-group
    both sides singleton     gray        #dddddd
    gene not covered         light gray  #eeeeee

Overlays: diagonal hatch marks a gene covered by two or more pieces; a gold
asterisk marks an internal break, positioned where the alignment split.

Rows group by status, then cluster within each group on the internal-break
marker matrix, so visually identical rows sit together.

Only --config is required. --settings defaults to settings.toml in the
current directory, --target defaults to the config's reference target (or
the only declared target; with several and no target given, you are asked),
and --work-dir/--out default to the standard settings
paths.work_dir/paths.output_root layout for that target and locus.
Set --target explicitly whenever the scan target differs from the
reference's own dataset -- otherwise the reference's database is used,
which gives the wrong row labels.
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path

_PKG_DIR = Path(__file__).resolve().parent
if str(_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(_PKG_DIR))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from heatmap_cluster import (
    informative_marker_positions as _informative_marker_positions,
    cluster_within_group as _cluster_within_group,
)
from config_utils import load_locus_cfg


# ── Constants ──────────────────────────────────────────────────────
SYNTENY_COLORS: dict[str, list[str]] = {
    "SYNTENIC":     ["#1f77b4", "#4eb3d3", "#7bccc4", "#a8ddb5"],   # blue-green ramp
    "INVERTED":     ["#6a3d9a", "#9467bd", "#b39ddb", "#d1c4e9"],   # purple ramp
    "REORDERED":    ["#b07a00", "#e0a800", "#f0c869", "#f5dca0"],   # amber ramp
    "TRANSLOCATED": ["#c43d3d", "#ff7f0e", "#ffa94d", "#ffd8a8"],   # red-orange ramp
}
# Alias so PIECE_COLORS still resolves to the syntenic-piece color family.
PIECE_COLORS     = SYNTENY_COLORS["SYNTENIC"]
COL_EMPTY        = "#eeeeee"
COL_BREAK        = "#e6b800"   # gold - internal break marker
COL_PSEUDO       = "#d62728"   # red border - GFF-confirmed pseudogene
COL_DIVERGENT    = "#f1c40f"   # yellow border - divergent gene signal
COL_PARTIAL      = "#ff7f0e"   # orange
DUP_HATCH        = "////"      # diagonal hatch for "multi-piece covers this gene"
EXC_HATCH        = "...."      # dotted hatch for 'annotated but inert' exception gene
                                # - kept colorless because color = piece identity
STAR_SIZE        = 90          # fixed star size; position inside cell encodes location
STATE_BORDER_WIDTH = 2.0       # thick border so pseudo / candidate read at a glance

STATUS_TINT = {
    "COMPLETE":    "#e2efda",
    "CONDITIONAL": "#ddebf7",
    "DIVERGENT":   "#eadcf8",
    "DECAYED":     "#fff2cc",
    "ABSENT":      "#d9d9d9",
    "UNKNOWN":     "#f3f3f3",
}
STATUS_ORDER = ["COMPLETE", "CONDITIONAL", "DIVERGENT", "DECAYED", "ABSENT", "UNKNOWN"]

# Flank pair state color map.
FLANK_STATE_TINT = {
    "FULL":    "#4daf4a",
    "HALF":    "#ffd400",
    "EDGE":    "#9aa3a8",
    "MISSING": "#e41a1c",
    "":        "#ffffff",
}


# ── Reference gene layout ──────────────────────────────────────────────────────
def aux_anchors(locus_cfg: dict) -> list[dict]:
    """Return aux locus anchors (role='aux') for separate column rendering."""
    anchors = locus_cfg.get("_auto", {}).get("anchors", [])
    out = []
    for a in anchors:
        if a.get("role") != "aux":
            continue
        fam = a.get("family") or a.get("locus_tag", "aux")
        out.append({
            "locus_tag": a.get("locus_tag", ""),
            "family":    fam,
            "label":     fam,
            "product":   a.get("product", ""),
            "source_genome": a.get("genome_acc", ""),
        })
    return out


def reference_genes(locus_cfg: dict) -> list[dict]:
    """
    Ordered list of reference cluster genes (cluster_L, inner, cluster_R)
    with coordinates RELATIVE to cluster_L.start, matching the query-space
    used by ruler.py's internal-break records.
    """
    anchors = locus_cfg.get("_auto", {}).get("anchors", [])
    cluster_roles = {"cluster_L", "cluster_R", "inner"}
    cl = next((a for a in anchors if a["role"] == "cluster_L"), None)
    origin = int(cl["start"]) if cl else 0
    rows = [a for a in anchors if a["role"] in cluster_roles]
    rows.sort(key=lambda a: int(a["start"]))
    out = []
    for a in rows:
        fam = a.get("family") or ""
        rs  = int(a["start"]) - origin
        re_ = int(a["end"])   - origin
        # Column label: family, else the first significant product word, else locus_tag.
        if fam:
            label = fam
        else:
            prod = (a.get("product") or "").strip()
            label = _short_product_label(prod) or a["locus_tag"]
        out.append({
            "locus_tag": a["locus_tag"],
            "family":    fam,
            "label":     label,
            "rel_start": rs,
            "rel_end":   re_,
            "exception": bool(a.get("exception")),
        })
    return out


_STOPWORDS_FOR_LABEL = {
    "family", "domain-containing", "domain", "containing", "protein",
    "putative", "type",
}

def _short_product_label(product: str) -> str:
    """First informative word of a GFF product description, truncated."""
    if not product:
        return ""
    p = product.replace("%2C", ",").replace("%3B", ";")
    for tok in p.split():
        clean = tok.strip(",;()[]").lower()
        if clean and clean not in _STOPWORDS_FOR_LABEL:
            return tok.strip(",;()[]")[:12]
    return ""


# ── Per-genome cell value derivation ──────────────────────────────────────────────────────
def piece_index_for_gene(pieces: list[dict], gene_tag: str) -> list[int]:
    """Indices of pieces that cover the given reference gene."""
    hits = []
    for i, p in enumerate(pieces or []):
        if gene_tag in (p.get("_genes") or []):
            hits.append(i)
    return hits


def _owning_piece_index(pieces: list[dict], gene: dict, s_range: dict | None = None) -> int | None:
    """Piece whose reference-q range overlaps the gene the most."""
    if not pieces:
        return None

    # Pass 1: Subject overlap (target coordinates)
    if s_range:
        best_i_s, best_ov_s = None, 0
        s_contig = s_range.get("contig")
        s_lo, s_hi = s_range.get("s_lo", 0), s_range.get("s_hi", 0)
        s_lo, s_hi = min(s_lo, s_hi), max(s_lo, s_hi)
        for i, p in enumerate(pieces):
            if p.get("sseqid") != s_contig:
                continue
            p_lo = int(p.get("sstart", 0))
            p_hi = int(p.get("send", 0))
            p_lo, p_hi = min(p_lo, p_hi), max(p_lo, p_hi)

            ov = max(0, min(s_hi, p_hi) - max(s_lo, p_lo))
            # Also accept proximity within 10 kb when there is no strict overlap.
            if ov > 0:
                if ov > best_ov_s:
                    best_i_s, best_ov_s = i, ov
            else:
                dist = max(0, max(s_lo, p_lo) - min(s_hi, p_hi))
                # Score nearby non-overlapping genes by proximity so they latch onto the nearest piece.
                if dist < 15000:
                    prox_score = 1.0 / (dist + 1.0)  # small positive value
                    if prox_score > best_ov_s:
                        best_i_s, best_ov_s = i, prox_score

        if best_i_s is not None:
            return best_i_s

    # Pass 2 & 3: Query coordinates
    rs, re_ = gene["rel_start"], gene["rel_end"]
    g_mid = (rs + re_) / 2.0
    best_i, best_ov = None, 0
    near_i, near_dist = None, float("inf")
    for i, p in enumerate(pieces):
        # Pieces store BLAST query coordinates in `qstart` / `qend`, in either order by strand.
        qs = int(p.get("qstart") or 0)
        qe = int(p.get("qend") or 0)
        q_lo, q_hi = (qs, qe) if qs <= qe else (qe, qs)
        if q_hi <= q_lo:
            continue
        ov = max(0, min(re_, q_hi) - max(rs, q_lo))
        if ov > best_ov:
            best_i, best_ov = i, ov
        # Distance from gene mid-point to closest end of this piece.
        p_mid = (q_lo + q_hi) / 2.0
        dist = abs(g_mid - p_mid)
        if dist < near_dist:
            near_i, near_dist = i, dist
    # Prefer overlap winner; fall back to nearest piece by midpoint.
    return best_i if best_i is not None else near_i


def break_fractions_in_gene(
    internal_breaks_per_piece: list[list[dict]],
    gene: dict,
) -> list[float]:
    """
    Return the fractional positions (0-1) of internal breaks that fall
    inside this reference gene.  Fraction = (break_midpoint - gene.rel_start)
    / gene_length.  This lets the renderer place the star at the actual
    position of the BLAST discontinuity, not just in the cell center.
    """
    rs, re_ = gene["rel_start"], gene["rel_end"]
    gene_len = max(1, re_ - rs)
    out: list[float] = []
    for piece_breaks in (internal_breaks_per_piece or []):
        for b in piece_breaks or []:
            if b.get("ref_gene") != gene["locus_tag"]:
                continue
            q_lo = int(b.get("q_gap_lo") or 0)
            q_hi = int(b.get("q_gap_hi") or 0)
            mid  = (q_lo + q_hi) / 2.0
            frac = (mid - rs) / gene_len
            frac = max(0.05, min(0.95, frac))   # keep inside the cell
            out.append(frac)
    return out


# ── Genome metadata loader (light) ──────────────────────────────────────────────────────
def load_genome_meta(db_path: str, accessions: list[str]) -> dict[str, dict]:
    if not db_path or not Path(db_path).exists() or not accessions:
        return {}
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    placeholders = ",".join("?" * len(accessions))
    rows = con.execute(
        f"SELECT accession, species, strain FROM genomes "
        f"WHERE accession IN ({placeholders})",
        accessions,
    ).fetchall()
    con.close()
    return {
        r["accession"]: {
            "species": r["species"] or "",
            "strain":  r["strain"]  or "",
        } for r in rows
    }


def row_label(acc: str, meta: dict) -> str:
    sp = meta.get("species", "") or acc
    parts = sp.split()
    abbr = f"{parts[0][0]}. {parts[1]}" if len(parts) > 1 else sp
    strain = meta.get("strain", "")[:18]
    return f"{abbr} {strain}".strip()


def _resolve_label_orientations(items) -> dict[str, int]:
    """Per-genome heatmap gene-order direction (+1 ascending / -1 descending),
    with synteny propagation for undeterminable genomes.

    Base call = the rigorous flank-relative orientation
    (``ruler._flank_orientation``: +1 co-oriented, -1 inverted, 0 = cannot be
    determined directly, e.g. multi-contig / displaced / rearranged flanks).

    For a base-0 genome we do NOT blindly assume reference order (that can mis-
    number a genome that is actually inverted). Instead we propagate from other
    QUERY genomes:

      1. If it shares a *distinctive* flank-synteny group (a REF_LT / sky-blue
         group - genomes whose flank carries the same non-canonical neighbor
         gene), inherit the signed consensus of the determinable members of
         that group. This is the "if it's sky-blue and those strains read
         right=n, read right=n too" rule.
      2. Otherwise fall back to the GLOBAL consensus of all determinable query
         genomes at this locus (still query-derived, not a reference guess).
      3. Only a genuinely isolated, no-signal genome defaults to +1.

    The conservative ``inversion`` flag in genome_summary.csv is unaffected --    it reports only directly-confirmed inversions (``_flank_relative_inversion``).
    """
    from ruler import _flank_orientation
    from synteny import compute_ctx_groups, _ctx_as_frozenset

    accs = [acc for acc, _res, _m, _s in items]
    per_genome_ctx = [(_res.get("_piece_side_contexts") or [])
                      for _acc, _res, _m, _s in items]
    base = {acc: _flank_orientation(_res) for acc, _res, _m, _s in items}

    raw_to_group, _ = compute_ctx_groups(per_genome_ctx)

    group_votes: dict = {}          # REF_LT group id -> signed vote sum
    genome_gids: dict = {}          # acc -> set of REF_LT group ids it touches
    global_vote = 0
    for acc, grp in zip(accs, per_genome_ctx):
        gids = set()
        for l, r in grp:
            for gid in (raw_to_group.get((_ctx_as_frozenset(l), "left")),
                        raw_to_group.get((_ctx_as_frozenset(r), "right"))):
                if gid and gid[0] == "REF_LT":   # only groups that get a distinct color
                    gids.add(gid)
        genome_gids[acc] = gids
        b = base[acc]
        if b != 0:
            v = 1 if b > 0 else -1
            global_vote += v
            for gid in gids:
                group_votes[gid] = group_votes.get(gid, 0) + v

    resolved: dict[str, int] = {}
    for acc in accs:
        b = base[acc]
        if b != 0:
            resolved[acc] = b
            continue
        s = sum(group_votes.get(gid, 0) for gid in genome_gids[acc])
        if s != 0:
            resolved[acc] = 1 if s > 0 else -1
        elif global_vote != 0:
            resolved[acc] = 1 if global_vote > 0 else -1
        else:
            resolved[acc] = 1
    return resolved


# ── Main render ──────────────────────────────────────────────────────
def render(ruler_results: dict[str, dict],
           locus_cfg: dict,
           genome_meta: dict[str, dict],
           out_path: Path,
           include_unknown: bool = False,
           max_pieces: int | None = None) -> None:
    """Render the strain x ref-gene heatmap."""

    ref_genes = reference_genes(locus_cfg)
    if not ref_genes:
        sys.exit("[heatmap] reference cluster has no cluster_L/inner/cluster_R anchors")
    aux_list = aux_anchors(locus_cfg)
    n_aux = len(aux_list)
    AUX_GAP = 0.6  # visual gap (in column-widths) between cluster cols and aux cols

    # Collect genomes
    items = []
    dropped_n_pieces = 0
    for acc, res in ruler_results.items():
        status = res.get("status", "UNKNOWN")
        if status == "UNKNOWN" and not include_unknown:
            continue
        n_pieces = len(res.get("_pieces") or [])
        if max_pieces is not None and n_pieces > max_pieces:
            dropped_n_pieces += 1
            continue
        meta = genome_meta.get(acc, {})
        items.append((acc, res, meta, status))
    if dropped_n_pieces:
        print(f"[heatmap] dropped {dropped_n_pieces} genome(s) with > {max_pieces} pieces")

    # Group by status, then cluster within each group by internal-break pattern.
    markers = _informative_marker_positions(ruler_results)
    by_status: dict[str, list] = {}
    for it in items:
        by_status.setdefault(it[3], []).append(it)

    sorted_items: list = []
    for status in STATUS_ORDER:
        grp = by_status.get(status)
        if not grp:
            continue
        # Pre-sort alphabetically as deterministic tiebreaker, then cluster.
        grp.sort(key=lambda t: (t[2].get("species", ""),
                                 t[2].get("strain", ""), t[0]))
        sorted_items.extend(_cluster_within_group(grp, markers))
    # Statuses outside STATUS_ORDER are appended after the ordered ones.
    for status, grp in by_status.items():
        if status not in STATUS_ORDER:
            sorted_items.extend(grp)
    items = sorted_items

    # Per-genome gene-order direction, propagated where it cannot be called directly.
    label_orient = _resolve_label_orientations(items)

    n_rows = len(items)
    n_cols = len(ref_genes)
    if n_rows == 0:
        sys.exit("[heatmap] no genomes to plot")

    # Figure size scales with data - include aux columns in width budget
    _total_col_w = n_cols + (AUX_GAP + n_aux if n_aux else 0)
    fig, ax = plt.subplots(figsize=(max(10, 0.45 * _total_col_w + 6),
                                    max(6,  0.22 * n_rows + 2)))

    # Three strips left of the body: status, left synteny, right synteny.
    X_STATUS = -4.00
    X_LFP    = -3.00
    X_RFP    = -2.00
    CHIP_W       = 0.85
    CHIP_HW      = CHIP_W / 2.0
    FP_CHIP_W    = 0.85          # total width of a segmented FP chip
    FP_CHIP_HW   = FP_CHIP_W / 2.0

    # Build per-piece palettes for L and R contexts.
    from synteny import (
        shared_side_palette, group_counts_lookup, render_context_label,
        is_canonical_side,
    )
    all_piece_ctx = [res.get("_piece_side_contexts") or []
                     for acc, res, _, _ in items]
    canonical_blue = SYNTENY_COLORS["SYNTENIC"][0]
    # Contexts sharing a ref_lt token are unioned into one group and one color.
    shared_pal = shared_side_palette(
        all_piece_ctx, canonical_color=canonical_blue)
    _group_count = group_counts_lookup(all_piece_ctx)

    def _strip_edge(ctx):
        return frozenset(ctx) - {"EDGE"}

    def _seg_color(ctx, side):
        return shared_pal.get((frozenset(ctx), side), "#dddddd")

    def _group_n(ctx, side):
        """Group-total count for the given raw ctx, used as the
        cell-fill frequency winner key."""
        return _group_count.get((frozenset(ctx), side), 0)

    # ── Piece cell-fill rule
    _GREY      = "#dddddd"   # singleton / unknown
    _INNER_KEY = frozenset({"INNER"})

    def _piece_fill(l_ctx, r_ctx):
        l_key = _strip_edge(l_ctx)
        r_key = _strip_edge(r_ctx)
        # Group-level counts (union-find merged) - see synteny.compute_ctx_groups.
        l_n   = _group_n(l_ctx, "left")
        r_n   = _group_n(r_ctx, "right")
        l_col = _seg_color(l_ctx, "left")
        r_col = _seg_color(r_ctx, "right")
        # Rule 1: any canonical edge anchors the piece to a canonical locus.
        if (is_canonical_side(l_ctx, "left")
                or is_canonical_side(r_ctx, "right")):
            return canonical_blue
        # Rule 2: INNER deference - yields to the more specific colored side.
        if l_key == _INNER_KEY and r_col != _GREY:
            return r_col
        if r_key == _INNER_KEY and l_col != _GREY:
            return l_col
        # Rule 2.5: non-gray beats gray regardless of frequency.
        if l_col == _GREY and r_col != _GREY:
            return r_col
        if r_col == _GREY and l_col != _GREY:
            return l_col
        # Rule 3: frequency wins; L wins ties.
        return l_col if l_n >= r_n else r_col

    def _status_display(res, status: str) -> tuple[str, str]:
        """Visible row-strip label/tint without changing the status count."""
        return status[:4], STATUS_TINT.get(status, "#ffffff")

    # ── Status strip + segmented fingerprint chips
    for i, (acc, res, meta, status) in enumerate(items):
        status_label, tint = _status_display(res, status)
        ax.add_patch(plt.Rectangle((X_STATUS - CHIP_HW, i - CHIP_HW),
                                   CHIP_W, CHIP_W,
                                   facecolor=tint, edgecolor="#888",
                                   linewidth=0.4))
        ax.text(X_STATUS, i, status_label,
                ha="center", va="center", fontsize=5.5, color="#333")

        per_piece_ctx = res.get("_piece_side_contexts") or []
        n_pieces = len(per_piece_ctx)
        for cx, side_key, side_name in ((X_LFP, 0, "left"),
                                          (X_RFP, 1, "right")):
            if n_pieces == 0:
                # Nothing to color: draw an empty box so the column still reads as a chip slot.
                ax.add_patch(plt.Rectangle(
                    (cx - FP_CHIP_HW, i - CHIP_HW),
                    FP_CHIP_W, CHIP_W,
                    facecolor="#ffffff", edgecolor="#888",
                    linewidth=0.3))
                continue
            seg_w = FP_CHIP_W / n_pieces
            for k, ctx_pair in enumerate(per_piece_ctx):
                col = _seg_color(ctx_pair[side_key], side_name)
                ax.add_patch(plt.Rectangle(
                    (cx - FP_CHIP_HW + k * seg_w, i - CHIP_HW),
                    seg_w, CHIP_W,
                    facecolor=col, edgecolor="#888",
                    linewidth=0.15))
            # Outer border around the full chip so multi-segment chips stay visually grouped.
            ax.add_patch(plt.Rectangle(
                (cx - FP_CHIP_HW, i - CHIP_HW),
                FP_CHIP_W, CHIP_W,
                facecolor="none", edgecolor="#444",
                linewidth=0.35))

    # ── Gene cells

    for i, (acc, res, meta, status) in enumerate(items):
        # Final ABSENT may still carry pieces or pseudogene cells; only UNKNOWN suppresses content.
        is_unavailable = status == "UNKNOWN"
        pieces = [] if is_unavailable else (res.get("_pieces") or [])
        breaks = res.get("_internal_breaks") or []
        gene_states = res.get("_gene_states") or {}
        gene_states_by_family = res.get("_gene_states_by_family") or {}
        gene_s_ranges_by_family = res.get("_gene_s_ranges_by_family") or {}
        per_piece_ctx = [] if is_unavailable else (res.get("_piece_side_contexts") or [])
        cell_piece_owner: dict[int, int] = {}
        for jj, gg in enumerate(ref_genes):
            cov = piece_index_for_gene(pieces, gg["locus_tag"])
            if cov:
                cell_piece_owner[jj] = sorted(cov)[0]
                continue
            ffam = gg.get("family") or ""
            gstate = (
                gene_states_by_family.get(ffam)
                if ffam and ffam in gene_states_by_family
                else gene_states.get(gg["locus_tag"], "ABSENT")
            )
            if gstate not in ("PSEUDOGENE", "DIVERGENT"):
                continue
            rel = (
                (res.get("_gene_relations_by_family") or {}).get(ffam, {})
                if ffam else {}
            )
            explicit_piece_idx = rel.get("piece_idx")
            s_range = (
                gene_s_ranges_by_family.get(ffam)
                if ffam and ffam in gene_s_ranges_by_family
                else res.get("_gene_s_ranges", {}).get(gg["locus_tag"])
            )
            try:
                owner = int(explicit_piece_idx)
            except (TypeError, ValueError):
                owner = _owning_piece_index(pieces, gg, s_range=s_range)
            if owner is not None:
                cell_piece_owner[jj] = owner

        # Gene-order numbering follows the flank-relative orientation, falling back to
        # reference order when that orientation is undetermined.
        reverse_all = label_orient.get(acc, 1) < 0
        piece_order_labels: dict[int, int] = {}
        for p_idx in sorted(set(cell_piece_owner.values())):
            members = [jj for jj, owner in cell_piece_owner.items() if owner == p_idx]
            if not members:
                continue
            labels = (range(len(members), 0, -1) if reverse_all
                      else range(1, len(members) + 1))
            for jj, lab in zip(sorted(members), labels):
                piece_order_labels[jj] = lab
        for j, g in enumerate(ref_genes):
            covering   = piece_index_for_gene(pieces, g["locus_tag"])
            fam = g.get("family") or ""
            gene_state = (
                gene_states_by_family.get(fam)
                if fam and fam in gene_states_by_family
                else gene_states.get(g["locus_tag"], "ABSENT")
            )
            hatch = ""
            if len(covering) >= 2:
                hatch = DUP_HATCH
            elif g.get("exception"):
                hatch = EXC_HATCH

            # Fill: synteny-palette color of the covering / owning piece.
            def _fill_for(piece_idx: int) -> str:
                if piece_idx < len(per_piece_ctx):
                    l_ctx, r_ctx = per_piece_ctx[piece_idx]
                    return _piece_fill(l_ctx, r_ctx)
                return canonical_blue

            if covering:
                # DUPLICATED loci: covering = [0, 1, ...].
                if len(covering) == 1:
                    fill = _fill_for(covering[0])
                else:
                    _fills = [_fill_for(p) for p in sorted(covering)]
                    if canonical_blue in _fills:
                        fill = canonical_blue
                    else:
                        fill = next(
                            (f for f in _fills if f != _GREY and f != COL_EMPTY),
                            _fills[0],
                        )
            elif gene_state in ("PSEUDOGENE", "DIVERGENT"):
                # _owning_piece_index falls back to the nearest piece, so None means no pieces.
                rel = (
                    (res.get("_gene_relations_by_family") or {}).get(fam, {})
                    if fam else {}
                )
                explicit_piece_idx = rel.get("piece_idx")
                s_range = (
                    gene_s_ranges_by_family.get(fam)
                    if fam and fam in gene_s_ranges_by_family
                    else res.get("_gene_s_ranges", {}).get(g["locus_tag"])
                )
                try:
                    owner = int(explicit_piece_idx)
                except (TypeError, ValueError):
                    owner = _owning_piece_index(pieces, g, s_range=s_range)
                fill = _fill_for(owner) if owner is not None else COL_EMPTY
            else:
                fill = COL_EMPTY

            # Border: state takes priority.
            if gene_state == "PSEUDOGENE":
                edge_col, edge_lw = COL_PSEUDO, STATE_BORDER_WIDTH
            elif gene_state == "DIVERGENT":
                edge_col, edge_lw = COL_DIVERGENT, STATE_BORDER_WIDTH
            elif not covering:
                edge_col, edge_lw = "#bbb", 0.3
            else:
                edge_col, edge_lw = "#333", 0.4

            # Two passes: fill at zorder 1, then the border on top at zorder 2.
            ax.add_patch(plt.Rectangle(
                (j - 0.45, i - 0.45), 0.9, 0.9,
                facecolor=fill, edgecolor="none",
                hatch=hatch, zorder=1))
            ax.add_patch(plt.Rectangle(
                (j - 0.45, i - 0.45), 0.9, 0.9,
                facecolor="none", edgecolor=edge_col,
                linewidth=edge_lw, zorder=2))

            # Star markers: one per internal break, at its fractional position within the gene.
            fracs = break_fractions_in_gene(breaks, g)
            for f in fracs:
                ax.scatter(
                    [j - 0.4 + f * 0.8], [i],
                    s=STAR_SIZE,
                    marker="*", color=COL_BREAK,
                    edgecolors="#7a5b00", linewidths=0.6, zorder=5,
                )
            if j in piece_order_labels and fill != COL_EMPTY:
                ax.text(
                    j, i, str(piece_order_labels[j]),
                    ha="center", va="center",
                    fontsize=5.5, fontweight="bold",
                    color="#111111", zorder=6,
                )

    # ── Aux locus cells (right of cluster cols, after AUX_GAP)
    if n_aux:
        aux_x0 = n_cols + AUX_GAP   # x of first aux column
        # Vertical separator line between cluster cols and aux cols
        ax.axvline(x=n_cols - 0.5 + AUX_GAP / 2, color="#666",
                   linewidth=0.8, linestyle=":", zorder=0)
        for i, (acc, res, meta, status) in enumerate(items):
            aux_hits = res.get("_aux_hits", {}) or {}
            for k, a in enumerate(aux_list):
                fam = a["family"]
                hit = aux_hits.get(fam)
                jx = aux_x0 + k
                if hit and hit.get("state") == "INTACT":
                    fill = canonical_blue
                    edge_col, edge_lw = "#333", 0.4
                elif hit and hit.get("state") == "DIVERGENT":
                    fill = "#aec7e8"     # light blue (degraded ortholog)
                    edge_col, edge_lw = COL_DIVERGENT, STATE_BORDER_WIDTH
                else:
                    fill = COL_EMPTY
                    edge_col, edge_lw = "#bbb", 0.3
                ax.add_patch(plt.Rectangle(
                    (jx - 0.45, i - 0.45), 0.9, 0.9,
                    facecolor=fill, edgecolor="none", zorder=1))
                ax.add_patch(plt.Rectangle(
                    (jx - 0.45, i - 0.45), 0.9, 0.9,
                    facecolor="none", edgecolor=edge_col,
                    linewidth=edge_lw, zorder=2))

    # ── Axis labels & limits
    _xmax = (n_cols - 0.3) if not n_aux else (n_cols + AUX_GAP + n_aux - 0.3)
    ax.set_xlim(-4.55, _xmax)
    ax.set_ylim(-0.6, n_rows - 0.2)

    # Strip column headers and gene names share one rotated xtick set so they align.
    strip_positions = [X_STATUS, X_LFP, X_RFP]
    strip_labels    = ["Status", "L_synteny", "R_synteny"]
    aux_positions = ([n_cols + AUX_GAP + k for k in range(n_aux)]
                     if n_aux else [])
    aux_labels    = [f"[aux] {a['label']}" for a in aux_list] if n_aux else []
    ax.set_xticks(strip_positions + list(range(n_cols)) + aux_positions)
    ax.set_xticklabels(strip_labels + [g["label"] for g in ref_genes] + aux_labels,
                       rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(
        [row_label(acc, meta) for acc, _, meta, _ in items],
        fontsize=6,
    )
    ax.invert_yaxis()
    _ref_acc = (locus_cfg.get("reference", {}).get("accession")
                or locus_cfg.get("_auto", {}).get("built_from", {}).get("accession")
                or "reference")
    ax.set_xlabel(f"Reference cluster gene  ({_ref_acc} family label)", fontsize=9)
    ax.set_title(f"{locus_cfg.get('locus_id','')}  — per-gene coverage "
                 f"({n_rows} genomes × {n_cols} reference genes)", fontsize=11)

    # ── Legend (categorised)
    def _section(title: str):
        return mpatches.Patch(facecolor="none", edgecolor="none",
                              label=f"── {title} ──")
    def _blank():
        return mpatches.Patch(facecolor="none", edgecolor="none", label="")

    legend = [
        _section("Cell fill: piece flank synteny"),
        mpatches.Patch(facecolor=canonical_blue, edgecolor="#333",
                       label="Both flanks canonical (= reference)"),
        mpatches.Patch(facecolor="#aec7e8", edgecolor="#333",
                       label="One side shared with ≥1 other strain\n"
                             "(sky-blue = group shared with chip)"),
        mpatches.Patch(facecolor="#dddddd", edgecolor="#333",
                       label="Both sides singleton / no anchor"),
        mpatches.Patch(facecolor=COL_EMPTY, edgecolor="#bbb",
                       label="Gene not covered"),
        _blank(),

        # 2. Cell border - gene state
        _section("Cell border: gene state"),
        mpatches.Patch(facecolor=SYNTENY_COLORS["SYNTENIC"][0], edgecolor=COL_PSEUDO,
                       linewidth=STATE_BORDER_WIDTH,
                       label="Confirmed pseudogene (red border)"),
        mpatches.Patch(facecolor=SYNTENY_COLORS["SYNTENIC"][0], edgecolor=COL_DIVERGENT,
                       linewidth=STATE_BORDER_WIDTH,
                       label="Decayed gene signal (yellow border)"),
        _blank(),

        # 3. Cell overlay
        _section("Cell overlay"),
        mpatches.Patch(facecolor=SYNTENY_COLORS["SYNTENIC"][0], edgecolor="#333",
                       hatch=DUP_HATCH,
                       label="≥2 pieces cover this gene (hatch)"),
        plt.scatter([], [], s=STAR_SIZE, marker="*",
                    color=COL_BREAK, edgecolors="#7a5b00",
                    label="Internal break  — (position inside cell\n"
                          "= location of break in the gene)"),
        _blank(),

        # 4. Status strip
        _section("Status strip (leftmost chip)"),
        mpatches.Patch(facecolor=STATUS_TINT["COMPLETE"],   edgecolor="#888", label="COMPLETE"),
        mpatches.Patch(facecolor=STATUS_TINT["CONDITIONAL"], edgecolor="#888", label="CONDITIONAL"),
        mpatches.Patch(facecolor=STATUS_TINT["DIVERGENT"],  edgecolor="#888", label="DIVERGENT"),
        mpatches.Patch(facecolor=STATUS_TINT["DECAYED"],    edgecolor="#888", label="DECAYED"),
        mpatches.Patch(facecolor=STATUS_TINT["ABSENT"],     edgecolor="#888", label="ABSENT"),
        mpatches.Patch(facecolor=STATUS_TINT["UNKNOWN"],    edgecolor="#888", label="UNKNOWN"),
        _blank(),

        # 5. L_synteny / R_synteny chip strips (left of heatmap)
        _section("L_synteny / R_synteny chip strips"),
        mpatches.Patch(facecolor=SYNTENY_COLORS["SYNTENIC"][0], edgecolor="#666",
                       label="Canonical = reference flank\n"
                             "(any of L1/L2/L1+L2 on left,\n"
                             " R1/R2/R1+R2 on right)"),
        mpatches.Patch(facecolor="#aec7e8", edgecolor="#666",
                       label="Shared non-canonical group\n"
                             "(≥2 strains share this flank;\n"
                             " sky-blue = palette priority 1)"),
        mpatches.Patch(facecolor="#dddddd", edgecolor="#666",
                       label="Singleton / no anchor / no GFF"),
        mpatches.Patch(facecolor="none", edgecolor="#444",
                       label="Chip is segmented N-way for N pieces\n"
                             "(left→right within chip = piece order)"),
    ]
    ax.legend(handles=legend, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), fontsize=7, frameon=False,
              handlelength=1.6, handletextpad=0.6, labelspacing=0.5)

    plt.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[heatmap] wrote {out_path}")


# ── CLI ──────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=_DESCRIPTION,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config",   required=True, help="Locus JSON config")
    ap.add_argument("--settings", default="settings.toml",
                    help="Path to settings.toml (default: settings.toml in "
                         "the current directory). Needed to resolve --target, "
                         "the metadata DB, and the default --work-dir/--out.")
    ap.add_argument("--target",   default=None,
                    help="Target name in settings.toml. Default: the config's "
                         "reference target, or the only declared target; "
                         "with several and no target given, you are asked. "
                         "Set this explicitly when scanning a target "
                         "different from the reference's own dataset.")
    ap.add_argument("--db",       default=None,
                    help="(optional) SQLite DB for species/strain row labels, "
                         "overriding the one resolved from --target")
    ap.add_argument("--work-dir", default=None,
                    help="LocusRuler work dir for this locus, containing "
                         "ruler_results.json (default: settings paths.work_dir"
                         "/<target>/<locus_id>)")
    ap.add_argument("--out",      default=None,
                    help="Output PNG path (default: settings paths.output_root"
                         "/<target>/<locus_id>/cluster_heatmap.png)")
    ap.add_argument("--output-dir", default=None,
                    help="Override output_root from settings; only affects "
                         "the default --out (ignored if --out is also given)")
    ap.add_argument("--include-unknown", action="store_true",
                    help="Plot UNKNOWN-status genomes too")
    ap.add_argument("--max-pieces", type=int, default=4,
                    help="Drop genomes whose accepted piece count exceeds N "
                         "(default: 4).  Cross-genus runs can fragment a "
                         "single cluster into many short HSPs; capping at 4 "
                         "removes those noisy rows.  Pass 0 (or any large "
                         "number) to disable.")
    args = ap.parse_args()

    cfg_path = Path(args.config).resolve()
    locus_cfg = load_locus_cfg(cfg_path)
    locus_id = locus_cfg.get("locus_id")

    settings = None
    settings_path = Path(args.settings)
    needs_settings = not (args.work_dir and args.out and args.db)
    if needs_settings or args.target:
        if not settings_path.exists():
            sys.exit(f"ERROR: settings file not found: {settings_path}\n"
                      f"       Pass --settings <path>, or supply --work-dir, "
                      f"--out and --db to skip it entirely.")
        from config_utils import load_settings, get_target_cfg, resolve_target_name
        settings = load_settings(settings_path)

    target_name = args.target
    if target_name is None and settings is not None:
        target_name = resolve_target_name(settings, None, locus_cfg)

    if args.work_dir:
        work_dir = Path(args.work_dir).resolve()
    else:
        work_dir = Path(settings["paths"]["work_dir"]) / target_name / locus_id

    if args.out:
        out_path = Path(args.out).resolve()
    else:
        output_root = (
            Path(args.output_dir) if args.output_dir
            else Path(settings["paths"]["output_root"])
        )
        out_path = output_root / target_name / locus_id / "cluster_heatmap.png"

    rrjson = work_dir / "ruler_results.json"
    if not rrjson.exists():
        sys.exit(f"ruler_results.json not found: {rrjson}\n"
                  f"       Run locus-ruler for this config/target first.")

    with open(rrjson) as f:
        ruler_results = json.load(f)

    # Row-label DB, in order: --db, then --target (resolved via settings), then the config.
    db_path = args.db
    if not db_path and target_name and settings is not None:
        db_path = get_target_cfg(settings, target_name)["db"]
    if not db_path:
        db_path = (locus_cfg.get("_auto", {})
                            .get("built_from", {})
                            .get("db"))
    accessions = sorted(ruler_results.keys())
    genome_meta = load_genome_meta(db_path, accessions) if db_path else {}

    # Quality check: warn when most genomes resolve to an empty species, leaving accession labels.
    if genome_meta:
        missing = sum(1 for a in accessions
                      if not genome_meta.get(a, {}).get("species"))
        if missing > len(accessions) * 0.5:
            sys.stderr.write(
                f"[heatmap] WARN: {missing}/{len(accessions)} genomes have "
                f"no species/strain in the DB ({db_path}).  Row labels will "
                f"fall back to accessions.  Did you mean --target NAME "
                f"--settings settings.toml?\n")

    # The filter is n_pieces > max_pieces, so treat 0 or less as "no cap".
    max_pieces = args.max_pieces if args.max_pieces and args.max_pieces > 0 else None
    render(ruler_results, locus_cfg, genome_meta, out_path,
           include_unknown=args.include_unknown,
           max_pieces=max_pieces)


if __name__ == "__main__":
    main()
