"""Censoring-aware comparison of assay-matched measurements."""

from __future__ import annotations

import math


def _interval(measurement: dict, scale: str) -> tuple[float, float]:
    relation = measurement["relation"]
    value = float(measurement["value"])
    if scale == "log10":
        if value <= 0:
            raise ValueError("log10 outcome comparison requires positive values")
        value = math.log10(value)
    elif scale != "linear":
        raise ValueError(f"Unsupported comparison scale: {scale}")
    if relation == "=":
        return value, value
    if relation == "<":
        return -math.inf, value
    if relation == ">":
        return value, math.inf
    raise ValueError(f"Unsupported relation: {relation}")


def _display_bound(value: float) -> float | None:
    return None if math.isinf(value) else round(value, 6)


def compare_measurements(
    parent: dict,
    child: dict,
    *,
    higher_is_better: bool,
    equivalence_margin: float = 0.0,
    comparison_scale: str = "linear",
) -> dict:
    """Compare like-for-like values without turning bounds into point values."""
    comparability_fields = ("endpoint", "state", "units", "assay_context")
    mismatches = [
        field
        for field in comparability_fields
        if parent.get(field) != child.get(field)
    ]
    if mismatches:
        return {
            "assay_comparability": "not_comparable",
            "mismatched_fields": mismatches,
            "classification": "indeterminate",
            "delta_lower": None,
            "delta_upper": None,
        }

    parent_low, parent_high = _interval(parent, comparison_scale)
    child_low, child_high = _interval(child, comparison_scale)
    delta_low = child_low - parent_high
    delta_high = child_high - parent_low
    if not higher_is_better:
        delta_low, delta_high = -delta_high, -delta_low

    if delta_low > equivalence_margin:
        classification = "improved"
    elif delta_high < -equivalence_margin:
        classification = "worsened"
    elif (
        math.isfinite(delta_low)
        and math.isfinite(delta_high)
        and delta_low >= -equivalence_margin
        and delta_high <= equivalence_margin
    ):
        classification = "comparable"
    else:
        classification = "indeterminate"

    return {
        "assay_comparability": "exact_context_match",
        "mismatched_fields": [],
        "classification": classification,
        "delta_lower": _display_bound(delta_low),
        "delta_upper": _display_bound(delta_high),
        "delta_orientation": "child_minus_parent; positive_is_better",
        "comparison_scale": comparison_scale,
        "censoring_preserved": (
            parent["relation"] != "=" or child["relation"] != "="
        ),
    }
