#!/usr/bin/env python3
"""Cell-state aggregation for the final locus status."""

VALID_STATUS_ROLES = {"CORE", "ASSOCIATED", "SHARED", "CONTEXT", "IGNORE"}
SCORED_STATUS_ROLES = {"CORE", "ASSOCIATED", "SHARED"}
VALID_GENE_STATES = {"INTACT", "DIVERGENT", "PSEUDOGENE", "ABSENT"}


def default_status_role(role: str, exception: bool = False) -> str:
    """Return the default status_role for an anchor row that omits one."""
    if exception:
        return "IGNORE"
    if role in {"cluster_L", "cluster_R", "inner"}:
        return "CORE"
    return "CONTEXT"


def normalize_status_role(value: str, role: str, exception: bool = False) -> str:
    normalized = (value or "").strip().upper()
    if not normalized:
        return default_status_role(role, exception)
    if normalized not in VALID_STATUS_ROLES:
        allowed = ", ".join(sorted(VALID_STATUS_ROLES))
        raise ValueError(
            f"invalid status_role '{value}' for role={role!r}; expected one of {allowed}"
        )
    return normalized


def classify_locus_status(
    anchor_states: dict[str, str],
    anchors: list[dict],
) -> tuple[str, dict]:
    """Aggregate curated anchor cell states into the final locus status."""
    scored: list[tuple[str, str, str]] = []
    core: list[tuple[str, str, str]] = []
    shared: list[tuple[str, str, str]] = []

    for anchor in anchors:
        role = normalize_status_role(
            anchor.get("status_role", ""),
            anchor.get("role", ""),
            bool(anchor.get("exception")),
        )
        if role not in SCORED_STATUS_ROLES:
            continue
        locus_tag = anchor.get("locus_tag", "")
        state = (anchor_states.get(locus_tag) or "ABSENT").upper()
        if state not in VALID_GENE_STATES:
            raise ValueError(
                f"invalid gene state '{state}' for status anchor '{locus_tag}'"
            )
        item = (locus_tag, role, state)
        scored.append(item)
        if role == "CORE":
            core.append(item)
        elif role == "SHARED":
            shared.append(item)

    if not core:
        raise ValueError("anchor CSV must contain at least one status_role=CORE row")

    core_lost = {"PSEUDOGENE", "ABSENT"}
    if all(state in core_lost for _, _, state in core):
        status = "ABSENT"
    elif any(state in core_lost for _, _, state in core):
        status = "DECAYED"
    elif any(state == "DIVERGENT" for _, _, state in core):
        status = "DIVERGENT"
    elif any(state != "INTACT" for _, _, state in shared):
        status = "CONDITIONAL"
    else:
        status = "COMPLETE"

    detail = {
        "scored": [
            {"locus_tag": locus_tag, "status_role": role, "state": state}
            for locus_tag, role, state in scored
        ],
        "counts": {
            state: sum(1 for _, _, observed in scored if observed == state)
            for state in sorted(VALID_GENE_STATES)
        },
    }
    return status, detail


def apply_locus_coverage_floor(
    status: str,
    coverage: float | None,
    min_decayed_coverage: float = 0.30,
) -> str:
    """Demote tiny locus fragments that only reached DECAYED by gene signal."""
    if status != "DECAYED":
        return status
    if coverage is None:
        return status
    if coverage < min_decayed_coverage:
        return "ABSENT"
    return status


def summarize_status_detail(detail: dict) -> dict[str, int]:
    """Return output counts for the same curated cells tracked by status."""
    counts = detail.get("counts") or {}
    expected = sum(int(counts.get(state, 0)) for state in VALID_GENE_STATES)
    absent = int(counts.get("ABSENT", 0))
    return {
        "genes_found": expected - absent,
        "genes_expected": expected,
        "pseudo_count": int(counts.get("PSEUDOGENE", 0)),
    }
