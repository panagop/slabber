"""Signed color scale for interactive maps and PyNite renders."""

FIELD_COLORS = ("#335E9F", "#FFFFFF", "#B15324")


def field_extent(low: float, high: float) -> float:
    return max(abs(low), abs(high), 1e-12)
