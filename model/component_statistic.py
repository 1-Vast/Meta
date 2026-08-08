"""Gauge-separated component statistic retained from the F6I research branch.

This module does not construct the biological surface and is not exported from
``model``.  It only preserves the verified algebraic separation

    prediction(P, L; S) = biological_surface(P, L) + location(S)

where the support location is one-dimensional, bounded, permutation invariant,
and independent of the protein representation.  External affinity admission is
still required before this statistic may be connected to the production state.
"""
from __future__ import annotations

import numpy as np


def support_location(
    residual: np.ndarray,
    *,
    ridge: float,
    bound: float = 0.5,
) -> float:
    """Return the bounded scalar location identified by a support multiset."""
    values = np.asarray(residual, dtype=np.float64)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("residual must be a non-empty finite vector")
    if not np.isfinite(ridge) or ridge < 0.0:
        raise ValueError("ridge must be finite and non-negative")
    if not np.isfinite(bound) or bound <= 0.0:
        raise ValueError("bound must be finite and positive")
    location = float(values.sum()) / (values.size + float(ridge))
    return float(np.clip(location, -float(bound), float(bound)))


def component_prediction(
    biological_surface: np.ndarray,
    support_residual: np.ndarray,
    *,
    ridge: float,
    bound: float = 0.5,
) -> tuple[np.ndarray, float]:
    """Add the support location to a finite protein-dependent query surface."""
    surface = np.asarray(biological_surface, dtype=np.float64)
    if surface.ndim != 1 or surface.size == 0 or not np.all(np.isfinite(surface)):
        raise ValueError("biological_surface must be a non-empty finite vector")
    location = support_location(support_residual, ridge=ridge, bound=bound)
    return surface + location, location
