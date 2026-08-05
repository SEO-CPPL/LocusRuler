#!/usr/bin/env python3
"""Per-piece per-side context analysis."""

from collections import Counter
from typing import Optional, Iterable, Any

# ── internal label mapping ──────────────────────────────────────────────────────
# Order for sorted join: L2+L1|...|R1+R2
_LABEL_ORDER = ["L2", "L1", "INNER", "R1", "R2", "EDGE"]


def _piece_span(piece: dict) -> tuple[int, int, bool]:
    """(min, max, is_plus_strand) from a piece/HSP dict."""
    s, e = int(piece["sstart"]), int(piece["send"])
    return min(s, e), max(s, e), e > s


def row_fingerprints(
    contexts: list[tuple[frozenset[str], frozenset[str]]],
) -> tuple[tuple[frozenset[str], ...], tuple[frozenset[str], ...]]:
    L = tuple(c[0] for c in contexts); R = tuple(c[1] for c in contexts)
    return L, R


# ── rendering ──────────────────────────────────────────────────────
def render_context_label(ctx: frozenset[str]) -> str:
    if not ctx: return "?"
    return "+".join(sorted(ctx, key=lambda x: _LABEL_ORDER.index(x) if x in _LABEL_ORDER else 99))


def render_fingerprint(fp: tuple[frozenset[str], ...]) -> str:
    return "|".join(render_context_label(c) for c in fp)


# ── Palette assignment ──────────────────────────────────────────────────────
_GROUP_PALETTE: tuple[str, ...] = (
    "#aec7e8", "#2ca02c", "#98df8a", "#9467bd", "#c5b0d5",
    "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#17becf", "#9edae5",
)
_OTHER_COLOR     = "#dddddd"
_CANONICAL_COLOR = "#eaeaea"

_CANONICAL_L_SETS = frozenset({frozenset({"L1"}), frozenset({"L2"}), frozenset({"L1", "L2"})})
_CANONICAL_R_SETS = frozenset({frozenset({"R1"}), frozenset({"R2"}), frozenset({"R1", "R2"})})

def is_canonical_side(ctx: frozenset[str], side: str) -> bool:
    stripped = frozenset(ctx) - {"EDGE"}
    return stripped in (_CANONICAL_L_SETS if side == "left" else _CANONICAL_R_SETS)

_TAB20_PAIRS = {
    "#1f77b4": "#aec7e8", "#ff7f0e": "#ffbb78", "#2ca02c": "#98df8a", "#d62728": "#ff9896",
    "#9467bd": "#c5b0d5", "#8c564b": "#c49c94", "#e377c2": "#f7b6d2", "#7f7f7f": "#c7c7c7",
    "#bcbd22": "#dbdb8d", "#17becf": "#9edae5",
}

def _ctx_as_frozenset(ctx) -> frozenset[str]:
    """Accept list / set / frozenset, return frozenset."""
    return ctx if isinstance(ctx, frozenset) else frozenset(ctx)


def _initial_group_key(raw: frozenset, side: str):
    """First-pass group key for a raw (EDGE-stripped) ctx."""
    s = (raw - {"EDGE"}) if isinstance(raw, frozenset) else _ctx_as_frozenset(raw) - {"EDGE"}
    if not s or s == frozenset({"?"}):
        return ("OTHER", side)
    canonical = {"L1", "L2"} if side == "left" else {"R1", "R2"}
    if s & canonical:
        return ("CANONICAL", side)
    if s == frozenset({"INNER"}):
        return ("INNER", side)
    return ("REF_LT", side, s)


def compute_ctx_groups(per_piece_contexts) -> tuple[dict, dict]:
    """Cluster all (raw_ctx, side) pairs into groups for palette assignment."""
    materialised = [list(grp) for grp in per_piece_contexts]

    unique: set[tuple[frozenset, str]] = set()
    for grp in materialised:
        for l, r in grp:
            unique.add((_ctx_as_frozenset(l), "left"))
            unique.add((_ctx_as_frozenset(r), "right"))

    # Bucket ctxs by initial key class so REF_LT ones can be union-find'd
    ref_lt_pairs: list[tuple[frozenset, str]] = []
    fixed_keys: dict[tuple[frozenset, str], tuple] = {}
    for raw, side in unique:
        k = _initial_group_key(raw, side)
        if k[0] == "REF_LT":
            ref_lt_pairs.append((raw, side))
        else:
            # Drop the side-disambiguator from CANONICAL / INNER / OTHER, keeping side itself.
            fixed_keys[(raw, side)] = k

    # Union-find on REF_LT ctxs: edge iff same side AND nonempty token overlap
    parent: dict = {p: p for p in ref_lt_pairs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, p1 in enumerate(ref_lt_pairs):
        s1 = p1[0] - {"EDGE"}
        for p2 in ref_lt_pairs[i + 1:]:
            if p1[1] != p2[1]:
                continue
            s2 = p2[0] - {"EDGE"}
            if s1 & s2:
                union(p1, p2)

    # Resolve final group id for every unique (raw_ctx, side)
    raw_to_group: dict[tuple[frozenset, str], tuple] = {}
    for raw, side in unique:
        if (raw, side) in fixed_keys:
            raw_to_group[(raw, side)] = fixed_keys[(raw, side)]
        else:
            root_raw, root_side = find((raw, side))
            # Group id = ("REF_LT", side, sorted root tokens) — hashable & stable
            raw_to_group[(raw, side)] = (
                "REF_LT", root_side,
                tuple(sorted(root_raw - {"EDGE"})),
            )

    # Group occurrence counts
    counts: Counter = Counter()
    for grp in materialised:
        for l, r in grp:
            counts[raw_to_group[(_ctx_as_frozenset(l), "left")]] += 1
            counts[raw_to_group[(_ctx_as_frozenset(r), "right")]] += 1

    return raw_to_group, dict(counts)


def shared_side_palette(
    per_piece_contexts,
    canonical_color: str = "#1f77b4",
    top_n: int | None = None,
) -> dict:
    """Single color map across L AND R contexts combined."""
    materialised = [list(grp) for grp in per_piece_contexts]
    raw_to_group, group_counts = compute_ctx_groups(materialised)

    # Stable ordering: highest count first, then side, then a member's label.
    def _example_label(group_id):
        for (raw, side), gid in raw_to_group.items():
            if gid == group_id:
                return render_context_label(_ctx_as_frozenset(raw) - {"EDGE"})
        return ""

    ordered_groups = sorted(
        group_counts.items(),
        key=lambda kv: (-kv[1], kv[0][1] if len(kv[0]) > 1 else "",
                        _example_label(kv[0])),
    )
    cap = top_n if top_n is not None else len(_GROUP_PALETTE)

    color_for_group: dict = {}
    palette_idx = 0
    for group_id, n in ordered_groups:
        kind = group_id[0]
        if kind == "CANONICAL":
            color_for_group[group_id] = canonical_color
            continue
        if kind == "OTHER":
            color_for_group[group_id] = _OTHER_COLOR
            continue
        if n < 2 or palette_idx >= cap:
            color_for_group[group_id] = _OTHER_COLOR
            continue
        color_for_group[group_id] = _GROUP_PALETTE[palette_idx]
        palette_idx += 1

    # Build raw → color map (every raw ctx inherits its group's color)
    final: dict = {}
    for (raw, side), gid in raw_to_group.items():
        final[(raw, side)] = color_for_group.get(gid, _OTHER_COLOR)
    return final


def group_counts_lookup(per_piece_contexts) -> dict[tuple[frozenset, str], int]:
    """Return ``{(raw_ctx, side): group_total_count}`` — every raw ctx maps
    to the total occurrences of its UNION-FIND GROUP (not its raw frozenset).
    Used by ``heatmap._piece_fill`` so the frequency-winner tiebreaker
    operates on the merged group's count, not the singleton raw count."""
    materialised = [list(grp) for grp in per_piece_contexts]
    raw_to_group, group_counts = compute_ctx_groups(materialised)
    return {
        (raw, side): group_counts.get(gid, 0)
        for (raw, side), gid in raw_to_group.items()
    }

