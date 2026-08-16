"""Coherent target-level Bayesian linear residual posterior (few-shot target adaptation).

One posterior per target support set (not a different local solve per query). For support
features Phi_s, residuals r_s = y_s - base(l_s), reliability weights R and observation noise
sigma_n^2, the Gaussian-linear posterior over the residual weights theta is

    Sigma = (Lambda + Phi_s^T R Phi_s / sigma_n^2)^{-1}
    m     = Sigma Phi_s^T R r_s / sigma_n^2

with Lambda a diagonal prior precision using SEPARATE positive precisions for the intercept,
linear-chemistry and nonlinear-chemistry feature blocks. Per query,

    predictive_mean     = base(l_q) + phi_q^T m
    epistemic_variance  = phi_q^T Sigma phi_q
    total_variance      = aleatoric_variance + epistemic_variance   (aleatoric from the head)

All solves are Cholesky-based; matrices are never explicitly inverted. Chemical similarity is
carried by the features (Option A), not by query-dependent support weights — so this is one
coherent target-level posterior, not a per-query weighted MAP.

k=0 (no support) gives m=0 and the predictive mean falls back exactly to the base predictor; the
epistemic variance is the prior predictive phi_q^T Lambda^{-1} phi_q.
"""
from __future__ import annotations

import torch
from torch import nn

JITTER = 1e-6  # symmetric PD jitter for the Cholesky factor
MAX_RELIABILITY = 1.0


def _validated_reliability(R: torch.Tensor | None, n: int, *, device, dtype) -> torch.Tensor:
    """Normalized confidence score in [0, 1]; zero means exact exclusion."""
    if R is None:
        return torch.ones(n, device=device, dtype=dtype)
    rel = R.to(device=device, dtype=dtype).reshape(-1)
    if rel.numel() != n:
        raise ValueError("support reliability length must match support rows")
    if not bool(torch.isfinite(rel).all()):
        raise ValueError("support reliability must be finite")
    if not bool((rel >= 0).all()):
        raise ValueError("support reliability must be non-negative")
    if not bool((rel <= MAX_RELIABILITY).all()):
        raise ValueError("support reliability is a normalized confidence score in [0, 1]")
    return rel


class BayesianResidualPosterior(nn.Module):
    def __init__(self, rep: int = 128, pdim: int = 64, quadratic: bool = False,
                 linear: bool = True,
                 prec_init: tuple[float, float, float] = (0.1, 1.0, 5.0),
                 noise_init: float = 1.0):
        super().__init__()
        self.rep = rep
        self.pdim = pdim
        self.linear = linear
        self.quadratic = quadratic
        # feature blocks: intercept (always) + linear chemistry + nonlinear chemistry (ablatable)
        self.feat_dim = 1 + (pdim if linear else 0) + (pdim if quadratic else 0)
        self.proj = nn.Linear(rep, pdim, bias=False)
        # separate positive prior precisions: [intercept, linear chemistry, nonlinear chemistry]
        self.log_prec = nn.Parameter(torch.log(torch.tensor(prec_init, dtype=torch.float32)))
        self.log_noise = nn.Parameter(torch.log(torch.tensor(float(noise_init))))

    # ---------------------------------------------------------------- features / prior
    def features(self, l_rep: torch.Tensor) -> torch.Tensor:
        """l_rep: (..., rep) -> (..., feat_dim). Design = [1, z, z^2] (quadratic) or [1, z]."""
        z = self.proj(l_rep)
        parts = [torch.ones_like(z[..., :1])]
        if self.linear:
            parts.append(z)
        if self.quadratic:
            parts.append(z * z)
        return torch.cat(parts, dim=-1)

    def prior_precision(self, device, dtype) -> torch.Tensor:
        p = torch.exp(self.log_prec.clamp(-8.0, 8.0))
        blocks = [p[0:1]]
        if self.linear:
            blocks.append(p[1].expand(self.pdim))
        if self.quadratic:
            blocks.append(p[2].expand(self.pdim))
        return torch.cat(blocks).to(device=device, dtype=dtype)

    def noise_var(self) -> torch.Tensor:
        return torch.exp(self.log_noise.clamp(-6.0, 6.0)) ** 2

    # ---------------------------------------------------------------- posterior solve
    def solve(self, Phi_s: torch.Tensor, r_s: torch.Tensor, R: torch.Tensor | None = None):
        """Return posterior mean m (feat_dim,) and Cholesky factor L of the precision A.

        Phi_s: (k, d) support features (k may be 0). r_s: (k,) residual targets.
        R: (k,) non-negative reliability weights (default ones).
        """
        d = self.feat_dim
        lam = self.prior_precision(Phi_s.device, Phi_s.dtype)          # (d,)
        sigma2 = self.noise_var().to(Phi_s.dtype)
        R = _validated_reliability(
            R, Phi_s.shape[0], device=Phi_s.device, dtype=Phi_s.dtype
        )
        RPhi = R.unsqueeze(-1) * Phi_s                                 # (k, d)
        A = torch.diag(lam) + (Phi_s.transpose(-1, -2) @ RPhi) / sigma2
        A = A + JITTER * torch.eye(d, device=Phi_s.device, dtype=Phi_s.dtype)
        b = (Phi_s.transpose(-1, -2) @ (R * r_s)) / sigma2            # (d,)
        L = torch.linalg.cholesky(A)                                  # A = L L^T, lower
        m = torch.cholesky_solve(b.unsqueeze(-1), L).squeeze(-1)      # (d,)
        return m, L

    def predict(self, Phi_q: torch.Tensor, m: torch.Tensor, L: torch.Tensor):
        """Predictive correction and epistemic variance for queries.

        mean_correction = Phi_q m ; epistemic = diag(Phi_q Sigma Phi_q^T) = || L^{-1} Phi_q^T ||^2.
        """
        mean_correction = Phi_q @ m                                   # (Q,)
        W = torch.linalg.solve_triangular(L, Phi_q.transpose(-1, -2), upper=False)  # (d, Q)
        epistemic = (W * W).sum(dim=0)                                # (Q,)
        return mean_correction, epistemic

    # ---------------------------------------------------------------- covariance (for tests/audit)
    def covariance(self, L: torch.Tensor) -> torch.Tensor:
        """Full posterior covariance Sigma = A^{-1} from its Cholesky factor (audit only)."""
        d = L.shape[-1]
        eye = torch.eye(d, device=L.device, dtype=L.dtype)
        return torch.cholesky_solve(eye, L)

    # ---------------------------------------------------------------- episode forward
    def forward(self, l_rep_q: torch.Tensor, l_rep_s: torch.Tensor | None,
                r_s: torch.Tensor | None, base_q: torch.Tensor,
                R_s: torch.Tensor | None = None) -> dict:
        Phi_q = self.features(l_rep_q)
        effective_k = 0
        if l_rep_s is not None:
            if r_s is None or r_s.reshape(-1).numel() != l_rep_s.shape[0]:
                raise ValueError("support residuals must match support representations")
            rel = _validated_reliability(
                R_s, l_rep_s.shape[0], device=l_rep_s.device, dtype=l_rep_s.dtype
            )
            active = rel > 0
            l_rep_s = l_rep_s[active]
            r_s = r_s.reshape(-1)[active]
            rel = rel[active]
            effective_k = int(active.sum().item())
        else:
            rel = None
        if effective_k == 0:                                        # strict k = 0 path
            zero = torch.zeros(0, self.feat_dim, device=l_rep_q.device, dtype=l_rep_q.dtype)
            m, L = self.solve(zero, torch.zeros(0, device=l_rep_q.device, dtype=l_rep_q.dtype))
        else:
            m, L = self.solve(self.features(l_rep_s), r_s, rel)
        mean_correction, epistemic = self.predict(Phi_q, m, L)
        return {"pred": base_q + mean_correction, "correction": mean_correction,
                "epistemic": epistemic, "m": m, "L": L,
                "effective_k": effective_k}


class JointPosterior(nn.Module):
    """Condition ligand and protein feature blocks in one Gaussian update."""

    @staticmethod
    def _validate(
        query: torch.Tensor,
        support: torch.Tensor,
        residual: torch.Tensor,
        precision: torch.Tensor,
        noise: torch.Tensor,
    ) -> torch.Tensor:
        if query.ndim != 2 or support.ndim != 2:
            raise ValueError("query and support designs must be matrices")
        if query.shape[1] != support.shape[1]:
            raise ValueError("query and support designs must have equal width")
        residual = residual.float().reshape(-1)
        if residual.numel() != support.shape[0]:
            raise ValueError("support residuals must match support rows")
        if precision.shape != (support.shape[1], support.shape[1]):
            raise ValueError("prior precision must match the design width")
        if noise.ndim != 0 or not bool(torch.isfinite(noise)) or float(noise) <= 0.0:
            raise ValueError("observation noise must be a positive scalar")
        return residual

    def condition(
        self,
        query: torch.Tensor,
        support: torch.Tensor,
        residual: torch.Tensor,
        precision: torch.Tensor,
        noise: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        """Return the exact finite-feature posterior mean and latent variance."""

        residual = self._validate(query, support, residual, precision, noise)
        eye = torch.eye(
            precision.shape[0], device=precision.device, dtype=precision.dtype
        )
        updated = precision + support.T @ support / noise
        chol = torch.linalg.cholesky(updated + JITTER * eye)
        rhs = support.T @ residual / noise
        weight = torch.cholesky_solve(rhs[:, None], chol).squeeze(1)
        mean = query @ weight
        solved = torch.linalg.solve_triangular(chol, query.T, upper=False)
        variance = solved.square().sum(dim=0)
        return {
            "mean": mean,
            "variance": variance,
            "weight": weight,
            "precisionchol": chol,
        }

    @staticmethod
    def logevidence(
        support: torch.Tensor,
        residual: torch.Tensor,
        covariance: torch.Tensor,
        noise: torch.Tensor,
    ) -> torch.Tensor:
        """Evaluate the support marginal likelihood up to a shared constant."""

        residual = residual.float().reshape(-1)
        if support.ndim != 2 or residual.numel() != support.shape[0]:
            raise ValueError("support design and residuals must align")
        if covariance.shape != (support.shape[1], support.shape[1]):
            raise ValueError("prior covariance must match the design width")
        rows = support.shape[0]
        eye = torch.eye(rows, device=support.device, dtype=support.dtype)
        marginal = support @ covariance @ support.T + noise * eye
        chol = torch.linalg.cholesky(marginal + JITTER * eye)
        solved = torch.cholesky_solve(residual[:, None], chol).squeeze(1)
        return -0.5 * (
            residual @ solved + 2.0 * torch.log(torch.diagonal(chol)).sum()
        )


__all__ = ["BayesianResidualPosterior", "JointPosterior"]


