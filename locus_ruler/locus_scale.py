#!/usr/bin/env python3
"""Distances expressed relative to the reference locus."""

from __future__ import annotations

DERIVE = -1
DEFAULT_REFERENCE_BP = 30_000


def reference_length(locus_cfg: dict, default: int = DEFAULT_REFERENCE_BP) -> int:
    """Reference cluster length in bp, as recorded by build_config."""
    auto = (locus_cfg or {}).get("_auto") or {}
    try:
        value = int(auto.get("cluster_ref_bp") or 0)
    except (TypeError, ValueError):
        value = 0
    return value or default


def resolve(configured, reference_bp: int, fraction: float, minimum: int = 0) -> int:
    """A configured window in bp; a negative value derives it from the reference."""
    try:
        value = int(configured)
    except (TypeError, ValueError):
        value = DERIVE
    if value >= 0:
        return value
    return max(minimum, int(round(reference_bp * fraction)))


def describe(name: str, configured, resolved: int, reference_bp: int) -> str:
    """One run-log line giving the effective window."""
    if _is_derived(configured):
        share = resolved / reference_bp if reference_bp else 0.0
        return (f"{name} = {resolved:,} bp "
                f"({share:.0%} of the {reference_bp:,} bp reference locus)")
    return f"{name} = {resolved:,} bp (fixed by settings)"


def _is_derived(configured) -> bool:
    try:
        return int(configured) < 0
    except (TypeError, ValueError):
        return True
