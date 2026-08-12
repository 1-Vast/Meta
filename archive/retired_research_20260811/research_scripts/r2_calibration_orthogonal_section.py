"""Small calibration-orthogonal ridge primitives for the R2 diagnosis.

These functions are research-only.  They make the per-target intercept explicit
and force the remaining ridge section to act on centered support residuals and
centered coordinates.  They do not authorize a production model change.
"""
from __future__ import annotations

import torch
from torch import Tensor


def ridge_support_weights(coordinates_s: Tensor, coordinates_q: Tensor,
                          ridge: float) -> Tensor:
    """Return query-by-support weights for a strictly positive dual ridge."""
    if ridge <= 0:
        raise ValueError("ridge must be strictly positive")
    if coordinates_s.ndim != 2 or coordinates_q.ndim != 2:
        raise ValueError("support and query coordinates must be matrices")
    if coordinates_s.shape[1] != coordinates_q.shape[1]:
        raise ValueError("support and query coordinate dimensions differ")
    gram = coordinates_s @ coordinates_s.T
    identity = torch.eye(len(gram), dtype=gram.dtype, device=gram.device)
    return torch.linalg.solve(
        gram + ridge * identity, coordinates_s @ coordinates_q.T
    ).T


def calibration_orthogonal_prediction(
    mu_s: Tensor,
    coordinates_s: Tensor,
    y_s: Tensor,
    mu_q: Tensor,
    coordinates_q: Tensor,
    ridge: float,
) -> tuple[Tensor, Tensor, Tensor]:
    """Return prediction, explicit calibration, and ligand-specific correction.

    Centering makes a constant support residual flow only through the explicit
    intercept.  The task-specific ridge state remains at most rank(M) <= k.
    """
    if len(y_s) == 0:
        raise ValueError("support cannot be empty")
    if len(mu_s) != len(y_s) or len(coordinates_s) != len(y_s):
        raise ValueError("support tensors have inconsistent lengths")
    if len(mu_q) != len(coordinates_q):
        raise ValueError("query tensors have inconsistent lengths")
    residual = y_s - mu_s
    calibration = residual.mean()
    centered_s = coordinates_s - coordinates_s.mean(dim=0, keepdim=True)
    centered_q = coordinates_q - coordinates_s.mean(dim=0, keepdim=True)
    weights = ridge_support_weights(centered_s, centered_q, ridge)
    specific = weights @ (residual - calibration)
    return mu_q + calibration + specific, calibration, specific
