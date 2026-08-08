"""Mathematical support utilities for the frozen model implementation.

This file consolidates law-class, loss, target, reference, and readout helpers.
It does not modify or replace the authoritative read-only theory.
"""
from __future__ import annotations

# Source section: lawclass
"""K(beta): the set-valued output of the frozen operator.

Contract row 2.4.

THEOREM_EXACT (asserted by the frozen theory, not re-derived here):
    K(beta) is a nonempty, compact, convex, W1-closed subset of Delta(V), and
    d_H^{W1}(K(beta), K(beta')) <= D_V^val * ||beta - beta'||_B + 2h.

ENGINEERING_CHOICE (this module):
    * the computational representation of K(beta) = (lower CDF, upper CDF, mesh h)
    * the quadrature   W1(P,P') = h * sum_{j<M} |F_j - F'_j|
      (exact for laws supported on the fixed grid, which is the class the band
      constrains)
    * the exact dynamic program for the directed Hausdorff distance.

PROHIBITED (contract row 2.4): replacing this set by one Gaussian, one variance,
or one interval width.
"""


from dataclasses import dataclass

import numpy as np

from . import bands


@dataclass
class LawClass:
    """Computational representation of K(beta) on the fixed mesh."""

    lower: np.ndarray      # (G,) lower CDF band
    upper: np.ndarray      # (G,) upper CDF band
    grid: np.ndarray       # (G,) fixed mesh points
    h: float               # fixed mesh -- never refined

    # -- required exposure ------------------------------------------------
    @property
    def lower_cdf(self) -> np.ndarray:
        return self.lower

    @property
    def upper_cdf(self) -> np.ndarray:
        return self.upper

    @property
    def valid(self) -> bool:
        return bands.is_valid(bands.join(self.lower, self.upper))

    @property
    def mesh(self) -> float:
        return self.h

    # -- extremal members (both are genuine CDFs on V) --------------------
    def upper_extremal_cdf(self) -> np.ndarray:
        """F = upper. Nondecreasing, in [0,1], F(a_max)=1, inside the band."""
        return self.upper.copy()

    def lower_extremal_cdf(self) -> np.ndarray:
        """F = lower, forced to reach 1 at the last grid point."""
        f = self.lower.copy()
        f[-1] = 1.0
        return np.maximum.accumulate(f)

    @staticmethod
    def cdf_to_pmf(F: np.ndarray) -> np.ndarray:
        return np.diff(np.concatenate([[0.0], F]))

    def contains(self, F: np.ndarray, tol: float = 1e-9) -> bool:
        F = np.asarray(F, dtype=np.float64)
        mono = np.all(np.diff(F) >= -tol)
        return bool(mono and F[-1] >= 1.0 - tol
                    and np.all(F >= self.lower - tol)
                    and np.all(F <= self.upper + tol))

    def width(self) -> float:
        return float(np.mean(self.upper - self.lower))

    def __repr__(self) -> str:
        return (f"K(beta): G={len(self.grid)}, h={self.h:.5f}, "
                f"valid={self.valid}, mean_width={self.width():.4f}")


def from_band(beta, grid: np.ndarray, h: float) -> LawClass:
    beta = np.asarray(beta, dtype=np.float64).reshape(-1)
    lo, up = bands.split(beta)
    return LawClass(lower=lo.copy(), upper=up.copy(), grid=np.asarray(grid), h=float(h))


# ----------------------------------------------------------------------
# W1 and Hausdorff-W1 between band classes
# ----------------------------------------------------------------------
def w1_between_cdfs(F: np.ndarray, Fp: np.ndarray, h: float) -> float:
    """W1(P,P') = int |F - F'| dv, exact for grid-supported laws."""
    return float(h * np.abs(F[:-1] - Fp[:-1]).sum())


def _candidate_levels(lo, up, lo2, up2) -> np.ndarray:
    lv = np.concatenate([lo, up, lo2, up2, [0.0, 1.0]])
    return np.unique(np.clip(lv, 0.0, 1.0))


def directed_hausdorff_w1(beta, beta2, h: float) -> float:
    """sup_{F in K(beta)} inf_{F' in K(beta2)} W1(F, F').

    Inner infimum is attained by the pointwise clip  F' = median(lo2, F, up2),
    which is feasible: it is monotone (pointwise median of monotone functions),
    lies in [lo2, up2], and equals 1 at the last grid point (up2[-1]=1 and
    F[-1]=1).  Hence

        inf_{F'} W1 = h * sum_j d_j(F_j),
        d_j(t) = max(lo2_j - t, t - up2_j, 0)   (convex, piecewise linear).

    The outer supremum is a separable maximisation over monotone F with
    lo <= F <= up.  Because each d_j is convex piecewise linear with breakpoints
    in {lo2_j, up2_j} and bounds in {lo_j, up_j}, an optimal solution exists on
    the finite candidate-level lattice, and the DP below is exact on it.
    """
    lo, up = bands.split(np.asarray(beta, dtype=np.float64).reshape(-1))
    lo2, up2 = bands.split(np.asarray(beta2, dtype=np.float64).reshape(-1))
    G = len(lo)
    levels = _candidate_levels(lo, up, lo2, up2)
    L = len(levels)

    NEG = -1e18
    dp = np.full(L, NEG)
    # j = 0
    feas = (levels >= lo[0] - 1e-12) & (levels <= up[0] + 1e-12)
    d0 = np.maximum.reduce([lo2[0] - levels, levels - up2[0], np.zeros(L)])
    dp[feas] = d0[feas]

    for j in range(1, G - 1):          # terms j = 0..G-2 enter the W1 quadrature
        run = np.maximum.accumulate(dp)          # best over levels <= current
        dj = np.maximum.reduce([lo2[j] - levels, levels - up2[j], np.zeros(L)])
        feas = (levels >= lo[j] - 1e-12) & (levels <= up[j] + 1e-12)
        nxt = np.where(feas, run + dj, NEG)
        nxt = np.where(run <= NEG / 2, NEG, nxt)
        dp = nxt

    best = float(np.max(dp))
    if best <= NEG / 2:
        raise RuntimeError("no feasible monotone CDF inside the band -- band invalid")
    return h * best


def hausdorff_w1(beta, beta2, h: float) -> float:
    """Exact (on the candidate lattice) Hausdorff-W1 distance between K(beta), K(beta')."""
    return max(directed_hausdorff_w1(beta, beta2, h),
               directed_hausdorff_w1(beta2, beta, h))


def theory_stability_bound(beta, beta2, h: float, D_V_val: float) -> float:
    """THEOREM_EXACT bound: d_H^{W1} <= D_V^val ||beta-beta'||_B + 2h."""
    return float(D_V_val * bands.band_sup_norm(beta, beta2) + 2.0 * h)


def sample_member_cdf(beta, rng: np.random.Generator) -> np.ndarray:
    """Draw a random monotone CDF inside the band (used to stress-test the DP)."""
    lo, up = bands.split(np.asarray(beta, dtype=np.float64).reshape(-1))
    G = len(lo)
    F = np.empty(G)
    prev = 0.0
    for j in range(G):
        a = max(lo[j], prev)
        b = up[j]
        if b < a:
            b = a
        F[j] = a + (b - a) * rng.random()
        prev = F[j]
    F[-1] = 1.0
    return np.maximum.accumulate(F)


# Source section: loss
"""The band loss L : B x V -> [0, inf), and the positive-ridge objective.

Contract rows 3.4 (L), 3.5 (L_0), 3.6 (J_mu), 4.3 (empirical risk).

Instantiation (ENGINEERING_CHOICE, verified in 04_training/loss_derivation.md):

    H_j(Y) = 1{v_j >= Y}                       (step CDF of the point target)
    c_j(beta,Y) = max(l_j - H_j, H_j - u_j, 0) (band-miss at grid point j)

    L(beta,Y) = (1/G) sum_j c_j  +  lambda_w * (1/G) sum_j (u_j - l_j)

Verified properties (all three are required by Foundations sec.4):
    convex in beta       -- c_j is a max of three affine functions; width is affine
    bounded              -- L <= 1 + lambda_w =: L_bar
    Lipschitz            -- |L(beta,Y)-L(beta',Y)| <= (1 + 2 lambda_w) ||beta-beta'||_B
                            uniformly in Y

Objective actually optimised (THEOREM_EXACT form, Target sec.3 / Meta-Learning sec.3):

    L(B(z) p, Y) + (mu/2) ||p||^2 ,     mu > 0.

NEGATIVE RIDGE IS FORBIDDEN (FAILURE_HISTORY sec.1).
"""


import numpy as np
import torch

from . import bands


# ----------------------------------------------------------------------
# step CDF of the point target
# ----------------------------------------------------------------------
def step_cdf(Y, grid):
    """H_j = 1{v_j >= Y}. Y: (...,) ; grid: (G,) -> (..., G)."""
    if isinstance(Y, torch.Tensor):
        g = grid if isinstance(grid, torch.Tensor) else torch.as_tensor(
            grid, dtype=Y.dtype, device=Y.device)
        return (g >= Y.unsqueeze(-1)).to(Y.dtype)
    g = np.asarray(grid, dtype=np.float64)
    Ya = np.asarray(Y, dtype=np.float64)
    return (g >= Ya[..., None]).astype(np.float64)


# ----------------------------------------------------------------------
# L
# ----------------------------------------------------------------------
def band_loss(beta, Y, grid, lambda_w: float):
    """L(beta, Y). beta: (..., 2G); Y: (...,). Returns (...,)."""
    lo, up = bands.split(beta)
    H = step_cdf(Y, grid)
    if isinstance(beta, torch.Tensor):
        zero = torch.zeros((), dtype=beta.dtype, device=beta.device)
        c = torch.maximum(torch.maximum(lo - H, H - up), zero)
        cover = c.mean(dim=-1)
        width = (up - lo).mean(dim=-1)
        return cover + lambda_w * width
    c = np.maximum(np.maximum(lo - H, H - up), 0.0)
    return c.mean(axis=-1) + lambda_w * (up - lo).mean(axis=-1)


def ridge(p, mu: float):
    """(mu/2)||p||^2 -- POSITIVE ridge on the coefficient vector."""
    if isinstance(p, torch.Tensor):
        return 0.5 * mu * (p * p).sum(dim=-1)
    return 0.5 * mu * (np.asarray(p) ** 2).sum(axis=-1)


def regularized_objective(beta, p, Y, grid, lambda_w: float, mu: float):
    """L(beta,Y) + (mu/2)||p||^2 -- the exact frozen per-task objective."""
    assert mu > 0.0, "positive ridge required"
    return band_loss(beta, Y, grid, lambda_w) + ridge(p, mu)


# ----------------------------------------------------------------------
# L_0 and the exact coefficient-linear form
# ----------------------------------------------------------------------
def L0_from_conditional_cdf(beta, Gcond, lambda_w: float):
    """L_0(z, beta) = E[L(beta,Y) | zeta = z], given G_j = P(Y <= v_j | z).

    Derivation (exact, see loss_derivation.md sec.5):
        E[c_j] = G_j (1 - u_j) + (1 - G_j) l_j
    hence L_0 is AFFINE in beta.
    """
    lo, up = bands.split(np.asarray(beta, dtype=np.float64))
    Gc = np.asarray(Gcond, dtype=np.float64)
    cover = (Gc * (1.0 - up) + (1.0 - Gc) * lo).mean(axis=-1)
    width = (up - lo).mean(axis=-1)
    return cover + lambda_w * width


def coefficient_linear_form(Bmat: np.ndarray, Gcond: np.ndarray, lambda_w: float):
    """Return (a, const) with L_0(z, B p) = a . p + const.

    Bmat: (2G, m+1) band-assembly matrix (columns are bands).
    Gcond: (G,) conditional CDF of Y given z on the fixed grid.
    """
    Bmat = np.asarray(Bmat, dtype=np.float64)
    G2, K = Bmat.shape
    G = G2 // 2
    lo = Bmat[:G, :]          # (G, K)
    up = Bmat[G:, :]          # (G, K)
    Gc = np.asarray(Gcond, dtype=np.float64).reshape(G, 1)
    a = ((up * (lambda_w - Gc)) + (lo * (1.0 - Gc - lambda_w))).mean(axis=0)
    const = float(Gc.mean())
    return a, const


# ----------------------------------------------------------------------
# simplex projection (used for the exact target and for reference fitting)
# ----------------------------------------------------------------------
def project_simplex(v: np.ndarray) -> np.ndarray:
    """Euclidean projection of v onto {p >= 0, sum p = 1} (Duchi et al.)."""
    v = np.asarray(v, dtype=np.float64)
    if v.ndim == 1:
        return _proj1(v)
    return np.stack([_proj1(row) for row in v.reshape(-1, v.shape[-1])]
                    ).reshape(v.shape)


def _proj1(v: np.ndarray) -> np.ndarray:
    n = v.shape[0]
    u = np.sort(v)[::-1]
    css = np.cumsum(u) - 1.0
    ind = np.arange(1, n + 1)
    cond = u - css / ind > 0
    rho = ind[cond][-1]
    theta = css[rho - 1] / rho
    return np.maximum(v - theta, 0.0)


# ----------------------------------------------------------------------
# empirical constants (numerical verification of L_bar and L_lip)
# ----------------------------------------------------------------------
def empirical_lipschitz(sample_bands, Ys, grid, lambda_w: float) -> float:
    """max_{i<j} |L(b_i,Y)-L(b_j,Y)| / ||b_i-b_j||_B over supplied probes."""
    worst = 0.0
    n = len(sample_bands)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.max(np.abs(sample_bands[i] - sample_bands[j])))
            if d < 1e-9:
                continue
            for Y in Ys:
                li = band_loss(sample_bands[i], np.array(Y), grid, lambda_w)
                lj = band_loss(sample_bands[j], np.array(Y), grid, lambda_w)
                worst = max(worst, abs(float(li) - float(lj)) / d)
    return worst


# Source section: target
"""The sole target g_mu^*(z) and the risks R_mu, E_mu.

Contract rows 3.6, 4.3, 4.8.

    J_mu(z,p) = L_0(z, B(z) p) + (mu/2)||p||^2 ,   g_mu^*(z) = argmin_{p in Delta_m} J_mu

For the instantiated loss, L_0 is AFFINE in beta (loss.py), hence

    L_0(z, B(z)p) = a(z).p + const   =>   J_mu = (mu/2)|| p + a/mu ||^2 + const
    =>  g_mu^*(z) = Proj_{Delta_m}( -a(z)/mu )      -- CLOSED FORM, exact.

This is a property of the chosen ENGINEERING_CHOICE loss, not an extra theory
claim: strong convexity, uniqueness and measurability are exactly what the
frozen theory proves, and the closed form is simply how they realise here.

The conditional CDF  G_j = P(Y <= v_j | zeta = z)  is required. It is
    * exact by enumeration on the analytic synthetic family (finite Z),
    * estimated by kernel-conditioning on continuous families (clearly labelled
      as a numerical estimate of the target, never as a theorem).
"""


import numpy as np
import torch



def target_from_conditional_cdf(Bmat: np.ndarray, Gcond: np.ndarray,
                                lambda_w: float, mu: float) -> np.ndarray:
    """Exact g_mu^*(z). Bmat: (2G, m+1); Gcond: (G,)."""
    assert mu > 0.0
    a, _ = coefficient_linear_form(Bmat, Gcond, lambda_w)
    return project_simplex(-a / mu)


def J_mu(Bmat: np.ndarray, Gcond: np.ndarray, p: np.ndarray,
         lambda_w: float, mu: float) -> float:
    a, const = coefficient_linear_form(Bmat, Gcond, lambda_w)
    return float(a @ p + const + 0.5 * mu * float(p @ p))


class ExactConditional:
    """Exact conditional CDF table for a finite-Z family.

    z values are hashed to a key; P(Y <= v | z) is computed by enumerating the
    generative model. Used by the analytic family in synthetic.py.
    """

    def __init__(self, grid: np.ndarray):
        self.grid = np.asarray(grid, dtype=np.float64)
        self.table = {}

    @staticmethod
    def key(z: np.ndarray, decimals: int = 9) -> bytes:
        return np.round(np.asarray(z, dtype=np.float64), decimals).tobytes()

    def set(self, z, Gcond):
        self.table[self.key(z)] = np.asarray(Gcond, dtype=np.float64)

    def get(self, z):
        return self.table[self.key(z)]

    def has(self, z) -> bool:
        return self.key(z) in self.table


class KernelConditional:
    """Nadaraya-Watson estimate of P(Y <= v | z) from a large reference pool.

    ENGINEERING_CHOICE / numerical estimate. Used only to *measure* the target
    on continuous families; never used inside the training objective, and never
    presented as theorem-backed.
    """

    def __init__(self, z_pool: np.ndarray, Y_pool: np.ndarray, grid: np.ndarray,
                 bandwidth: float = 0.06):
        self.Z = np.asarray(z_pool, dtype=np.float64)
        self.grid = np.asarray(grid, dtype=np.float64)
        self.H = (np.asarray(Y_pool, dtype=np.float64)[:, None] <= self.grid[None, :])
        self.H = self.H.astype(np.float64)
        self.bw = float(bandwidth)

    def __call__(self, z: np.ndarray) -> np.ndarray:
        z = np.asarray(z, dtype=np.float64)
        d2 = ((self.Z - z[None, :]) ** 2).sum(axis=1)
        w = np.exp(-d2 / (2.0 * self.bw ** 2))
        s = w.sum()
        if s < 1e-12:
            return self.H.mean(axis=0)
        Gc = (w[:, None] * self.H).sum(axis=0) / s
        return np.maximum.accumulate(np.clip(Gc, 0.0, 1.0))

    def batch(self, Zq: np.ndarray) -> np.ndarray:
        return np.stack([self(z) for z in np.asarray(Zq, dtype=np.float64)], axis=0)


# ----------------------------------------------------------------------
# risks
# ----------------------------------------------------------------------
def population_risk(Bmats: np.ndarray, Gconds: np.ndarray, P: np.ndarray,
                    lambda_w: float, mu: float, weights=None) -> float:
    """R_mu(F) = E_zeta[ J_mu(zeta, F(zeta)) ], evaluated on a sample of z."""
    n = len(P)
    vals = np.empty(n)
    for i in range(n):
        vals[i] = J_mu(Bmats[i], Gconds[i], P[i], lambda_w, mu)
    if weights is None:
        return float(vals.mean())
    w = np.asarray(weights, dtype=np.float64)
    return float((vals * w).sum() / w.sum())


def excess_risk(Bmats, Gconds, P, lambda_w, mu, weights=None) -> float:
    """E_mu(F) = R_mu(F) - R_mu(g_mu^*) >= 0."""
    n = len(P)
    star = np.stack([target_from_conditional_cdf(Bmats[i], Gconds[i], lambda_w, mu)
                     for i in range(n)], axis=0)
    return (population_risk(Bmats, Gconds, P, lambda_w, mu, weights)
            - population_risk(Bmats, Gconds, star, lambda_w, mu, weights))


def coefficient_l2_error(P: np.ndarray, Pstar: np.ndarray) -> float:
    """|| F - g_mu^* ||_{L^2(mu_zeta)} estimated on the sample."""
    d = np.linalg.norm(np.asarray(P) - np.asarray(Pstar), axis=-1)
    return float(np.sqrt((d ** 2).mean()))


# Source section: reference
"""The EXACT frozen hypothesis class H_N -- reference model.

Contract row 4.4. This is the class for which eps_approx(N) and Gamma_N are
proved. It is implemented literally:

    Omega_N = (Delta_m)^{N_N},  D_N = (m+1) nu_N
    phi_nu  = tensor-product piecewise-multilinear hats, phi_nu >= 0, sum = 1
    F_omega(z) = sum_nu phi_nu(z) omega_nu

CLAIM LEVEL: THEOREM_EXACT (this *is* H_N when ``coords`` spans all of Z).
When ``coords`` is a strict subset, the object is the exact H_N for the
*reduced* statistic z restricted to those coordinates -- still an exact frozen
sieve, but on a declared sub-cube. That is stated wherever it is used.

Parameterisation is ``projected``: node parameters are raw vectors projected
onto Delta_m, so the realised class is exactly Omega_N (boundary included), not
a dense interior subset.
"""


import numpy as np
import torch
import torch.nn as nn

from .config import MetaSieveConfig as DeploymentConfig
from .meta_operator import MultilinearSieve, assert_exact_hn_node_budget
from .runtime import DEFAULT_DTYPE, require_cuda


class FrozenHN(nn.Module):
    """F_omega(z) = sum_nu phi_nu(z) omega_nu, omega in Omega_N."""

    def __init__(self, cfg: DeploymentConfig, res_N: int, coords=None, seed: int = 0,
                 node_param: str = "projected", init_scale: float = 0.0,
                 device=None, dtype=DEFAULT_DTYPE):
        super().__init__()
        device = require_cuda(device)
        self.cfg = cfg
        self.coords = list(range(cfg.d_z)) if coords is None else list(coords)
        self.res_N = int(res_N)
        assert_exact_hn_node_budget(cfg, self.res_N, self.coords)
        self.sieve = MultilinearSieve(self.coords, res_N, cfg.n_coef,
                                      init_scale, seed, node_param=node_param,
                                      device=device, dtype=dtype)
        self.full_Z = (len(self.coords) == cfg.d_z
                       and set(self.coords) == set(range(cfg.d_z)))

    # -- frozen quantities ----------------------------------------------
    @property
    def nu_N(self) -> int:
        return self.sieve.n_nodes

    @property
    def D_N(self) -> int:
        return self.cfg.n_coef * self.nu_N

    def mesh_diameter(self) -> float:
        """max d_Z-diameter of a mesh cell, with d_Z = Euclidean on the cube."""
        d = len(self.coords)
        return float(np.sqrt(d) / self.res_N)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.sieve(z)

    def omega(self) -> torch.Tensor:
        return self.sieve.omega()

    # -- interpolation witness (Approximation sec.2) ---------------------
    def node_coordinates(self) -> np.ndarray:
        """(nu_N, d_coords) grid node positions in [0,1]^{d_coords}."""
        side = self.res_N + 1
        axes = [np.linspace(0.0, 1.0, side) for _ in self.coords]
        mesh = np.meshgrid(*axes, indexing="ij")
        return np.stack([mm.reshape(-1) for mm in mesh], axis=-1)

    def set_witness(self, target_at_nodes: np.ndarray) -> None:
        """omega_nu <- g_mu^*(nu). The exact witness of Approximation sec.2.

        Requires target_at_nodes in Delta_m; then omega in Omega_N, so the
        witness is a genuine member of H_N.
        """
        if self.sieve.node_param != "projected":
            raise ValueError("set_witness requires projected-node parameterisation")
        t = np.asarray(target_at_nodes, dtype=np.float64)
        assert t.shape == (self.nu_N, self.cfg.n_coef)
        assert np.all(t >= -1e-9) and np.allclose(t.sum(1), 1.0, atol=1e-8), \
            "witness values must lie in Delta_m"
        with torch.no_grad():
            self.sieve.node_logits.copy_(torch.as_tensor(
                t, dtype=self.sieve.node_logits.dtype,
                device=self.sieve.node_logits.device))

    def summary(self) -> dict:
        return {
            "class": "H_N (exact frozen sieve)",
            "coords": self.coords,
            "full_Z": self.full_Z,
            "res_N": self.res_N,
            "nu_N": self.nu_N,
            "D_N": self.D_N,
            "mesh_diameter": self.mesh_diameter(),
            "claim_level": "THEOREM_EXACT" if self.full_Z
                           else "THEOREM_EXACT (for the reduced statistic z|coords)",
        }


# Source section: readout
"""Scalar readouts from K(beta).

CLAIM LEVEL: ENGINEERING_CHOICE, and OUTSIDE the frozen operator.

The frozen theory defines NO unique scalar selector from K(beta) (contract row
4.1 / sec.6). Nothing here modifies A(F,z), nothing here enters the training
objective, and no theorem is claimed for any of it. Several readouts are
provided precisely so that sensitivity to the choice can be reported (T15).
"""


import numpy as np

from . import bands


def _quantile_from_cdf(F: np.ndarray, grid: np.ndarray, q: float = 0.5) -> float:
    idx = int(np.searchsorted(F, q, side="left"))
    idx = min(idx, len(grid) - 1)
    return float(grid[idx])


def _mean_from_cdf(F: np.ndarray, grid: np.ndarray) -> float:
    pmf = np.diff(np.concatenate([[0.0], F]))
    pmf = np.clip(pmf, 0.0, None)
    s = pmf.sum()
    return float((grid * pmf).sum() / s) if s > 0 else float(grid.mean())


def _assert_midband_member(lo: np.ndarray, up: np.ndarray, mid: np.ndarray) -> None:
    """Fail closed unless ``mid`` is a CDF member of the supplied band."""
    tol = 1e-9
    if not (np.all(np.isfinite(mid))
            and np.all(mid >= lo - tol)
            and np.all(mid <= up + tol)
            and np.all(np.diff(mid, axis=-1) >= -tol)
            and np.all(np.abs(mid[..., -1] - 1.0) <= tol)):
        raise AssertionError("midband CDF is not a member of K(beta)")


def midband_cdf(beta) -> np.ndarray:
    """Return the canonical midpoint CDF member of ``K(beta)``.

    The midpoint is constructed pointwise from the two band boundaries.  Only
    the terminal CDF value is set to one; global normalization would move the
    other entries outside the declared band.
    """
    beta = np.asarray(beta, dtype=np.float64)
    bands.assert_valid(beta, name="midband readout input")
    lo, up = bands.split(beta)
    mid = 0.5 * (lo + up)
    mid[..., -1] = 1.0
    _assert_midband_member(lo, up, mid)
    return mid


def _assert_midband_member_torch(lo, up, mid) -> None:
    """Torch equivalent of the NumPy member check used in batch evaluation."""
    import torch

    tol = 1e-9
    valid = (torch.isfinite(mid).all()
             and (mid >= lo - tol).all()
             and (mid <= up + tol).all()
             and (torch.diff(mid, dim=-1) >= -tol).all()
             and (mid[..., -1] - 1.0).abs().le(tol).all())
    if not bool(valid.item()):
        raise AssertionError("midband CDF is not a member of K(beta)")


def midband_cdf_torch(beta):
    """Torch batch form of :func:`midband_cdf` for CUDA evaluation."""
    import torch

    if not isinstance(beta, torch.Tensor):
        raise TypeError("beta must be a torch.Tensor")
    bands.assert_valid(beta, name="midband readout input")
    lo, up = bands.split(beta)
    mid = 0.5 * (lo + up)
    mid = mid.clone()
    mid[..., -1] = 1.0
    _assert_midband_member_torch(lo, up, mid)
    return mid


def readout_midband_median(beta, grid):
    mid = midband_cdf(beta)
    return _quantile_from_cdf(mid, grid, 0.5)


def readout_midband_mean(beta, grid):
    mid = midband_cdf(beta)
    return _mean_from_cdf(mid, grid)


def readout_midband_mean_torch(beta, grid):
    """Batch midpoint mean using the same valid member as NumPy readouts."""
    import torch

    mid = midband_cdf_torch(beta)
    if not isinstance(grid, torch.Tensor):
        grid = torch.as_tensor(grid, dtype=mid.dtype, device=mid.device)
    else:
        grid = grid.to(dtype=mid.dtype, device=mid.device)
    pmf = torch.diff(mid, dim=-1, prepend=torch.zeros_like(mid[..., :1])).clamp(min=0)
    return (grid * pmf).sum(-1) / pmf.sum(-1).clamp(min=1e-12)


def readout_extremal_midpoint(beta, grid):
    """Midpoint of the median interval [median(upper CDF), median(lower CDF)]."""
    lo, up = bands.split(np.asarray(beta, dtype=np.float64))
    lo_f = np.maximum.accumulate(lo.copy())
    lo_f[-1] = 1.0
    q_hi = _quantile_from_cdf(up, grid, 0.5)      # upper CDF -> smaller quantile
    q_lo = _quantile_from_cdf(lo_f, grid, 0.5)
    return 0.5 * (q_hi + q_lo)


def readout_upper_extremal_mean(beta, grid):
    _, up = bands.split(np.asarray(beta, dtype=np.float64))
    return _mean_from_cdf(up, grid)


READOUTS = {
    "midband_median": readout_midband_median,
    "midband_mean": readout_midband_mean,
    "extremal_midpoint": readout_extremal_midpoint,
    "upper_extremal_mean": readout_upper_extremal_mean,
}


def all_readouts(beta, grid) -> dict:
    return {k: f(beta, grid) for k, f in READOUTS.items()}
