"""The valid CDF-band polytope B, and convex assembly.

Contract rows 2.1 (B), 2.2 (B(z) assembly).

Band vector layout (ENGINEERING_CHOICE):
    beta = concat(lower[0..M], upper[0..M])  in R^{2(M+1)}
    ||beta||_B = max_r |beta_r|              (band sup norm)

Declared polytope (ENGINEERING_CHOICE; required: linear, bounded, implies
K(beta) nonempty):

    0 <= l_j <= u_j <= 1
    l_j <= l_{j+1},  u_j <= u_{j+1}
    u_M = 1

THEOREM_EXACT: B is compact convex, and p in Delta_m  =>  B(z)p in B.
"""

from __future__ import annotations

import numpy as np
import torch

TOL = 1e-9


# ----------------------------------------------------------------------
# splitting / joining
# ----------------------------------------------------------------------
def split(beta):
    """beta (..., 2G) -> (lower (..., G), upper (..., G))."""
    G = beta.shape[-1] // 2
    return beta[..., :G], beta[..., G:]


def join(lower, upper):
    if isinstance(lower, torch.Tensor):
        return torch.cat([lower, upper], dim=-1)
    return np.concatenate([lower, upper], axis=-1)


# ----------------------------------------------------------------------
# validity
# ----------------------------------------------------------------------
def validity_report(beta, tol: float = 1e-7) -> dict:
    """Exact check of every defining inequality of B. Works for np or torch."""
    if isinstance(beta, torch.Tensor):
        b = torch.atleast_2d(beta.detach())
        lo, up = split(b)
        d_lo = torch.diff(lo, dim=-1)
        d_up = torch.diff(up, dim=-1)
        rep = {
            "lower_nonneg": float(lo.amin().item()),
            "lower_le_upper": float((up - lo).amin().item()),
            "upper_le_one": float(up.amax().item()),
            "lower_monotone": float(d_lo.amin().item()) if d_lo.numel() else 0.0,
            "upper_monotone": float(d_up.amin().item()) if d_up.numel() else 0.0,
            "upper_terminal_dev": float((up[..., -1] - 1.0).abs().amax().item()),
        }
    else:
        b = np.atleast_2d(np.asarray(beta))
        lo, up = split(b)
        d_lo = np.diff(lo, axis=-1)
        d_up = np.diff(up, axis=-1)
        rep = {
            "lower_nonneg": float(np.min(lo)),
            "lower_le_upper": float(np.min(up - lo)),
            "upper_le_one": float(np.max(up)),
            "lower_monotone": float(np.min(d_lo)) if d_lo.size else 0.0,
            "upper_monotone": float(np.min(d_up)) if d_up.size else 0.0,
            "upper_terminal_dev": float(np.max(np.abs(up[..., -1] - 1.0))),
        }
    rep["valid"] = bool(
        rep["lower_nonneg"] >= -tol
        and rep["lower_le_upper"] >= -tol
        and rep["upper_le_one"] <= 1.0 + tol
        and rep["lower_monotone"] >= -tol
        and rep["upper_monotone"] >= -tol
        and rep["upper_terminal_dev"] <= tol
    )
    return rep


def is_valid(beta, tol: float = 1e-7) -> bool:
    return validity_report(beta, tol)["valid"]


def assert_valid(beta, name: str = "band", tol: float = 1e-7) -> None:
    rep = validity_report(beta, tol)
    if not rep["valid"]:
        raise AssertionError(f"{name} left the valid polytope B: {rep}")


# ----------------------------------------------------------------------
# construction recipes (ENGINEERING_CHOICE; validity is structural)
# ----------------------------------------------------------------------
def band_from_shape(u: np.ndarray, width: float) -> np.ndarray:
    """Given a nondecreasing shape u with u[0]>=0, u[-1]=1, return (max(0,u-w), u).

    Validity is automatic:
      lower = max(0, u - w) is nondecreasing, >= 0, <= u;
      u nondecreasing, <= 1, u[-1] = 1.
    """
    u = np.asarray(u, dtype=np.float64)
    u = np.maximum.accumulate(np.clip(u, 0.0, 1.0))
    u[-1] = 1.0
    lo = np.maximum(0.0, u - float(width))
    lo = np.maximum.accumulate(lo)
    return join(lo, u)


def logistic_shape(grid: np.ndarray, center: float, scale: float) -> np.ndarray:
    phi = 1.0 / (1.0 + np.exp(-(grid - center) / scale))
    denom = phi[-1] - phi[0]
    if denom <= 1e-12:
        return np.linspace(0.0, 1.0, len(grid))
    return (phi - phi[0]) / denom


def band_from_ecdf(ecdf: np.ndarray, eps: float) -> np.ndarray:
    """Monotone DKW-style envelope around an empirical CDF -> element of B."""
    ecdf = np.asarray(ecdf, dtype=np.float64)
    up = np.maximum.accumulate(np.clip(ecdf + eps, 0.0, 1.0))
    up[-1] = 1.0
    lo = np.maximum.accumulate(np.clip(ecdf - eps, 0.0, 1.0))
    lo = np.minimum(lo, up)
    lo = np.maximum.accumulate(lo)
    return join(lo, up)


# ----------------------------------------------------------------------
# convex assembly  beta = B p    (THEOREM_EXACT: stays in B)
# ----------------------------------------------------------------------
def assemble(B: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
    """B: (..., 2G, m+1) or (m+1, 2G); p: (..., m+1)  ->  beta: (..., 2G).

    No clipping. If the result is invalid, either a column is invalid or p is
    not in Delta_m -- both are bugs, not things to repair here.
    """
    if B.dim() == 2 and B.shape[0] == p.shape[-1]:
        # (m+1, 2G)
        return p @ B
    # (..., 2G, m+1) batched
    return torch.einsum("...rk,...k->...r", B, p)


def band_sup_norm(beta1, beta2):
    if isinstance(beta1, torch.Tensor):
        return (beta1 - beta2).abs().amax(dim=-1)
    return np.max(np.abs(np.asarray(beta1) - np.asarray(beta2)), axis=-1)
