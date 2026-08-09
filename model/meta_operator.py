"""Unchanged CSMO and Band operator implementation.

This file consolidates the previously separate operator modules. The simplex,
positive ridge, Band polytope, K, mesh, parameter names, and computations are
unchanged.
"""
from __future__ import annotations

# Source section: convexoperator
"""The Convex Sieve-Mixture Operator (CSMO) -- the primary trainable model.

Core innovation, stated in 01_derivation/core_innovation.md:

    F_{theta,omega}(z) = sum_a  pi_theta(z)_a * sum_nu phi^a_nu(P_a z) omega^{(a)}_nu

with  omega^{(a)}_nu in Delta_m,  phi^a_nu >= 0,  sum_nu phi^a_nu = 1,
      pi_theta(z) in Delta_{A-1}.

CLAIM LEVELS
    F(z) in Delta_m for every parameter value and every input : THEOREM_EXACT
        (nested convex combination of simplex points; same argument the frozen
         theory uses in Meta-Learning sec.2)
    CSMO contains the exact frozen H_N as a special case (A=1, P=id)  : exact
        parameter embedding implemented by CSMO.embed_frozen_hn
    The FULL class (A>1, trained gate) inherits neither eps_approx nor Gamma_N :
        THEOREM_COMPATIBLE ONLY.  No proof is manufactured.
"""


import hashlib
import json
from itertools import product

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from .config import MetaSieveConfig as DeploymentConfig, DEFAULT_VIEWS, Z_NAMES
from .runtime import DEFAULT_DTYPE, assert_finite, require_cuda

DTYPE = DEFAULT_DTYPE
MAX_EXACT_HN_NODES = 1_000_000


def assert_unit_statistic_domain(z: torch.Tensor, name: str = "statistic z") -> None:
    """Reject off-coverage statistics instead of silently projecting them."""
    assert_finite(z, name)
    if bool(((z < 0.0) | (z > 1.0)).any().item()):
        raise ValueError(f"{name} is outside the declared [0,1] domain")


def context_index(z: torch.Tensor, cfg: DeploymentConfig) -> torch.Tensor:
    """Map declared statistic coordinates to the finite deployment context."""
    if z.ndim != 2 or z.shape[1] != cfg.d_z:
        raise ValueError("z does not match the declared statistic dimension")
    assert_unit_statistic_domain(z)

    def digit(column: int, edges: tuple) -> torch.Tensor:
        boundaries = torch.as_tensor(edges, dtype=z.dtype, device=z.device)
        return torch.bucketize(z[:, column].contiguous(), boundaries)

    y = digit(12, cfg.kappa_edges_y)
    mass = digit(16, cfg.kappa_edges_mass)
    continuous = digit(26, cfg.kappa_edges_context_cont)
    n_mass = len(cfg.kappa_edges_mass) + 1
    n_continuous = len(cfg.kappa_edges_context_cont) + 1
    return (y * n_mass + mass) * n_continuous + continuous


def exact_hn_node_count(res_N: int, coords) -> int:
    if int(res_N) < 1:
        raise ValueError("res_N must be at least one")
    return (int(res_N) + 1) ** len(tuple(coords))


def assert_exact_hn_node_budget(cfg: DeploymentConfig, res_N: int, coords) -> None:
    """Reject impractical literal H_N constructions before allocating GPU memory."""
    n_nodes = exact_hn_node_count(res_N, coords)
    if n_nodes > MAX_EXACT_HN_NODES:
        node_bytes = n_nodes * cfg.n_coef * torch.empty((), dtype=DTYPE).element_size()
        raise ValueError(
            f"literal H_N needs {n_nodes:,} nodes ({node_bytes / 2**30:.2f} GiB for "
            f"node coefficients), above the {MAX_EXACT_HN_NODES:,}-node safety budget"
        )


def project_simplex_torch(v: torch.Tensor) -> torch.Tensor:
    """Euclidean projection of each row of v onto Delta_m. Differentiable a.e."""
    n = v.shape[-1]
    u, _ = torch.sort(v, dim=-1, descending=True)
    css = torch.cumsum(u, dim=-1) - 1.0
    ind = torch.arange(1, n + 1, dtype=v.dtype, device=v.device)
    cond = (u - css / ind) > 0
    rho = cond.to(v.dtype).cumsum(dim=-1).argmax(dim=-1, keepdim=True)
    theta = torch.gather(css, -1, rho) / (rho.to(v.dtype) + 1.0)
    return torch.clamp(v - theta, min=0.0)


# ----------------------------------------------------------------------
# exact piecewise-multilinear sieve on a frozen coordinate view
# ----------------------------------------------------------------------
class MultilinearSieve(nn.Module):
    """G_a(z) = sum_nu phi_nu(P_a z) omega_nu,  omega_nu in Delta_m.

    phi are the tensor-product hat functions on a regular grid of resolution r:
    nonnegative and summing to one -- exactly the basis of the frozen H_N.
    """

    def __init__(self, view_idx, res: int, n_coef: int, init_logit_scale: float = 0.0,
                 seed: int = 0, node_param: str = "softmax", device=None,
                 dtype=DTYPE):
        super().__init__()
        device = require_cuda(device)
        self.register_buffer("view_idx", torch.as_tensor(list(view_idx), dtype=torch.long,
                                                          device=device))
        self.d_a = len(view_idx)
        self.res = int(res)
        self.n_coef = int(n_coef)
        self.node_param = node_param
        side = self.res + 1
        self.n_nodes = side ** self.d_a if self.d_a > 0 else 1
        strides = [side ** (self.d_a - 1 - d) for d in range(self.d_a)]
        self.register_buffer("strides", torch.as_tensor(strides, dtype=torch.long,
                                                         device=device)
                             if self.d_a > 0 else torch.zeros(0, dtype=torch.long,
                                                              device=device))
        g = torch.Generator(device=device).manual_seed(seed)
        if node_param == "softmax":
            init = init_logit_scale * torch.randn(self.n_nodes, n_coef,
                                                  generator=g, dtype=dtype, device=device)
        elif node_param == "projected":
            init = torch.full((self.n_nodes, n_coef), 1.0 / n_coef, dtype=dtype,
                              device=device)
            init = init + init_logit_scale * torch.randn(self.n_nodes, n_coef,
                                                         generator=g, dtype=dtype,
                                                         device=device)
        else:
            raise ValueError(node_param)
        self.node_logits = nn.Parameter(init)
        # enumerate the 2^{d_a} corners once
        self._corners = list(product([0, 1], repeat=self.d_a)) if self.d_a > 0 else [()]

    def omega(self) -> torch.Tensor:
        """(n_nodes, m+1), each row exactly in Delta_m.

        'softmax'  -> image is the relative interior of Delta_m: the realised
                      class is a dense subset of H_N (each member IS in H_N).
        'projected'-> Euclidean projection onto Delta_m: the realised class is
                      EXACTLY H_N, boundary included.
        """
        if self.node_param == "softmax":
            return torch.softmax(self.node_logits, dim=-1)
        return project_simplex_torch(self.node_logits)

    def _forward_res1(self, z: torch.Tensor, om: torch.Tensor) -> torch.Tensor:
        """Exact tensor-product contraction, valid when res == 1.

        Mathematically identical to the corner sum, but O(d_a * 2^{d_a}) instead
        of O(2^{d_a} * d_a) Python iterations -- this is what makes the exact
        H_N on the full 14-dimensional Z trainable at r_N = 1.
        """
        B = z.shape[0]
        K = self.n_coef
        f = z.index_select(-1, self.view_idx)                           # (B, d_a)
        W = om.reshape(2, -1, K)
        f0 = f[:, 0].view(B, 1, 1)
        W = W[0].unsqueeze(0) * (1.0 - f0) + W[1].unsqueeze(0) * f0    # (B, rest, K)
        for d in range(1, self.d_a):
            W = W.reshape(B, 2, -1, K)
            fd = f[:, d].view(B, 1, 1)
            W = W[:, 0] * (1.0 - fd) + W[:, 1] * fd
        return W.reshape(B, K)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        assert_unit_statistic_domain(
            z.index_select(-1, self.view_idx), "sieve input coordinates"
        )
        om = self.omega()
        assert_finite(om, "sieve simplex nodes")
        if self.d_a == 0:
            return om[0].unsqueeze(0).expand(z.shape[0], -1)
        if self.res == 1 and self.d_a >= 6:
            return self._forward_res1(z, om)
        t = z.index_select(-1, self.view_idx) * self.res              # (B, d_a)
        c = torch.floor(t).clamp(0, self.res - 1)
        f = (t - c).clamp(0.0, 1.0)
        c = c.to(torch.long)
        out = torch.zeros(z.shape[0], self.n_coef, dtype=z.dtype, device=z.device)
        for corner in self._corners:
            wgt = torch.ones(z.shape[0], dtype=z.dtype, device=z.device)
            idx = torch.zeros(z.shape[0], dtype=torch.long, device=z.device)
            for d, bit in enumerate(corner):
                wgt = wgt * (f[:, d] if bit else (1.0 - f[:, d]))
                idx = idx + (c[:, d] + bit) * self.strides[d]
            out = out + wgt.unsqueeze(-1) * om.index_select(0, idx)
        return out


# ----------------------------------------------------------------------
# gate
# ----------------------------------------------------------------------
class Gate(nn.Module):
    """pi_theta : Z -> Delta_{A-1}. Deep. THEOREM_COMPATIBLE.

    Two exact simplex parameterisations are available (ENGINEERING_CHOICE; both
    give pi in Delta_{A-1} exactly, so F(z) in Delta_m is structural either way):

    'softmax'  pi = softmax(logits / T)
               components compete multiplicatively -> observed MIXTURE LOCK-IN
               (see 09_handoff/FAILURE_LEDGER.md F-2).

    'stick'    stick-breaking: lambda_a = sigmoid(u_a) independently per view,
               pi_a = lambda_a * prod_{b<a}(1 - lambda_b),  pi_A = prod(1-lambda_b).
               Equivalent to a convex residual cascade
                   p_a = (1-lambda_a) p_{a-1} + lambda_a G_a(z),
               so each view refines the running estimate instead of competing
               for a shared softmax budget.

    'normalized'  pi_a = softplus(u_a) / sum_b softplus(u_b). This is
                  permutation-equivariant in branch order and has no privileged
                  early position; it is provided for the branch-order audit.

    With n_out == 1 BOTH return exactly 1.0, so the H_N containment holds for
    either parameterisation.
    """

    def __init__(self, d_in: int, n_out: int, width: int, depth: int,
                 input_mask=None, seed: int = 0, floor: float = 0.0,
                 param: str = "normalized", device=None, dtype=DTYPE):
        super().__init__()
        device = require_cuda(device)
        if param not in ("softmax", "stick", "normalized"):
            raise ValueError(f"unknown gate parameterization {param!r}")
        self.floor = float(floor)
        self.param = param
        self.n_out = int(n_out)
        if self.n_out < 1:
            raise ValueError("gate requires at least one output")
        if self.n_out == 1:
            self.net = None
        else:
            devices = [device.index if device.index is not None
                       else torch.cuda.current_device()]
            with torch.random.fork_rng(devices=devices):
                torch.manual_seed(seed)
                layers = []
                d = d_in
                for _ in range(depth - 1):
                    layers += [nn.Linear(d, width, dtype=dtype, device=device), nn.SiLU()]
                    d = width
                n_head = n_out if param in ("softmax", "normalized") else n_out - 1
                last = nn.Linear(d, n_head, dtype=dtype, device=device)
                # Keep all branches reachable at initialisation.
                nn.init.normal_(last.weight, std=1e-3)
                if param in ("softmax", "normalized"):
                    nn.init.zeros_(last.bias)
                else:
                    with torch.no_grad():
                        for a in range(n_head):
                            lam = 1.0 / (n_out - a)
                            last.bias[a] = float(np.log(lam / (1.0 - lam)))
                layers += [last]
                self.net = nn.Sequential(*layers)
        # Optimisation-schedule temperature. NOT part of the declared model:
        # training ends with temperature == 1.0, so the deployed map is exactly
        # pi = softmax(logits). See 04_training/meta_training_algorithm.md.
        self.register_buffer("temperature", torch.tensor(1.0, dtype=dtype, device=device))
        if input_mask is None:
            m = torch.ones(d_in, dtype=dtype, device=device)
        else:
            m = torch.zeros(d_in, dtype=dtype, device=device)
            m[list(input_mask)] = 1.0
        self.register_buffer("input_mask", m)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        assert_finite(z, "gate input z")
        if self.n_out == 1:
            return torch.ones(z.shape[0], 1, dtype=z.dtype, device=z.device)
        h = self.net(z * self.input_mask)
        assert_finite(h, "gate logits")
        if self.param == "softmax":
            pi = torch.softmax(h / self.temperature, dim=-1)
        elif self.param == "normalized":
            score = F.softplus(h / self.temperature) + torch.finfo(h.dtype).eps
            pi = score / score.sum(dim=-1, keepdim=True)
        else:
            lam = torch.sigmoid(h / self.temperature)          # (B, A-1)
            rest = torch.cumprod(1.0 - lam, dim=-1)            # (B, A-1)
            first = lam[:, :1]
            mid = lam[:, 1:] * rest[:, :-1]
            lastp = rest[:, -1:]
            pi = torch.cat([first, mid, lastp], dim=-1)
        if self.floor > 0.0:
            # (1-eps)*pi + eps*uniform: still EXACTLY a probability vector
            # (convex combination of two probability vectors), but every sieve
            # keeps a nonzero gradient path, which prevents mixture lock-in.
            pi = (1.0 - self.floor) * pi + self.floor / self.n_out
        return pi


class FixedGate(nn.Module):
    """A frozen simplex gate for uniform/random routing controls."""

    def __init__(self, n_out: int, weights=None, seed: int = 0, device=None,
                 dtype=DTYPE):
        super().__init__()
        device = require_cuda(device)
        uniform = weights is None
        if uniform:
            weights = torch.full((n_out,), 1.0 / n_out, dtype=dtype, device=device)
        else:
            weights = torch.as_tensor(weights, dtype=dtype, device=device)
        if weights.shape != (n_out,) or not bool(torch.isfinite(weights).all().item()):
            raise ValueError("fixed gate weights must be a finite vector of length n_out")
        if not uniform:
            weights = project_simplex_torch(weights)
        self.register_buffer("weights", weights)
        self.register_buffer("temperature", torch.tensor(1.0, dtype=dtype, device=device))

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.weights.unsqueeze(0).expand(z.shape[0], -1)


# ----------------------------------------------------------------------
# default view bank (frozen coordinate selections)
# ----------------------------------------------------------------------


class CSMO(nn.Module):
    """The primary model: F(z) in Delta_m by construction."""

    def __init__(self, cfg: DeploymentConfig, views=None, head: str = "sieve",
                 gate_input_mask=None, seed: int = 0, node_param: str = "softmax",
                 res: int = None, init_scale: float = 0.5, gate_param: str = "normalized",
                 gate_mode: str = "learned", fixed_gate_weights=None, device=None,
                 dtype=DTYPE):
        super().__init__()
        device = require_cuda(device)
        self.cfg = cfg
        self.head = head
        res = cfg.view_res if res is None else int(res)
        views = DEFAULT_VIEWS[: cfg.n_views] if views is None else list(views)
        self.views = views
        if head == "sieve":
            # NOTE: a nonzero init scale is REQUIRED. With identical sieves the
            # gate gradient is exactly zero (dF/dpi = 0 when all G_a coincide),
            # which would make the gate an unlearnable component.
            self.sieves = nn.ModuleList([
                MultilinearSieve(v, res, cfg.n_coef, init_scale, seed + 17 * i,
                                 node_param=node_param, device=device, dtype=dtype)
                for i, v in enumerate(views)])
            if gate_mode == "learned":
                self.gate = Gate(cfg.d_z, len(views), cfg.gate_width, cfg.gate_depth,
                                 gate_input_mask, seed, floor=cfg.gate_floor,
                                 param=gate_param, device=device, dtype=dtype)
            elif gate_mode in ("uniform", "fixed"):
                self.gate = FixedGate(len(views), fixed_gate_weights, seed,
                                      device=device, dtype=dtype)
            else:
                raise ValueError(f"unknown gate_mode={gate_mode}")
        elif head == "direct":
            # Candidate B ablation: deep simplex regressor (no sieve bank)
            self.sieves = nn.ModuleList()
            # Candidate B is defined with a SOFTMAX link (candidate_realizations.md);
            # keep it that way so the ablation isolates the sieve bank.
            self.gate = Gate(cfg.d_z, cfg.n_coef, cfg.gate_width, cfg.gate_depth,
                             gate_input_mask, seed, floor=0.0, param="softmax",
                             device=device, dtype=dtype)
        else:
            raise ValueError(head)

    # ------------------------------------------------------------------
    def gate_weights(self, z: torch.Tensor) -> torch.Tensor:
        return self.gate(z)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """z: (B, d_z) -> p: (B, m+1) exactly in Delta_m."""
        if not z.is_cuda:
            raise RuntimeError("CSMO requires CUDA input")
        if z.ndim != 2 or z.shape[1] != self.cfg.d_z:
            raise ValueError("z does not match the declared statistic dimension")
        assert_unit_statistic_domain(z)
        if self.head == "direct":
            p = self.gate(z)
            assert_finite(p, "CSMO output")
            return p
        pi = self.gate(z)                                   # (B, A) in Delta_{A-1}
        p = torch.zeros(z.shape[0], self.cfg.n_coef, dtype=z.dtype, device=z.device)
        for a, sv in enumerate(self.sieves):
            p = p + pi[:, a].unsqueeze(-1) * sv(z)          # convex comb of Delta_m pts
        assert_finite(p, "CSMO output")
        return p

    @classmethod
    def exact_hn(cls, cfg: DeploymentConfig, res_N: int, omega=None, seed: int = 0,
                 device=None, dtype=DTYPE):
        """Construct the literal A=1, P=id, projected-node H_N restriction."""
        assert_exact_hn_node_budget(cfg, res_N, range(cfg.d_z))
        mdl = cls(cfg, views=[tuple(range(cfg.d_z))], head="sieve", seed=seed,
                  node_param="projected", res=res_N, init_scale=0.0,
                  gate_mode="uniform", device=device, dtype=dtype)
        if omega is not None:
            omega = torch.as_tensor(omega, dtype=dtype,
                                    device=mdl.sieves[0].node_logits.device)
            expected = mdl.sieves[0].node_logits.shape
            if omega.shape != expected:
                raise ValueError(f"omega shape {tuple(omega.shape)} != {tuple(expected)}")
            if bool((omega < 0).any().item()) or not torch.allclose(
                    omega.sum(dim=-1), torch.ones(omega.shape[0], dtype=dtype,
                                                    device=omega.device)):
                raise ValueError("omega must belong to (Delta_m)^nu_N")
            with torch.no_grad():
                mdl.sieves[0].node_logits.copy_(omega)
        return mdl

    @classmethod
    def embed_frozen_hn(cls, reference):
        """The explicit embedding iota(omega) from FrozenHN into CSMO."""
        if reference.sieve.node_param != "projected":
            raise ValueError("only projected-node FrozenHN has an exact CSMO embedding")
        device = reference.sieve.node_logits.device
        dtype = reference.sieve.node_logits.dtype
        mdl = cls(reference.cfg, views=[tuple(reference.coords)], head="sieve",
                  node_param="projected", res=reference.res_N, init_scale=0.0,
                  gate_mode="uniform", device=device, dtype=dtype)
        with torch.no_grad():
            # Copy the raw representative so both implementations apply the
            # same simplex projection exactly once.  Mathematically this is
            # iota(omega)=(A=1,P=id,omega); the raw copy also makes the two
            # floating-point evaluation graphs bitwise identical.
            mdl.sieves[0].node_logits.copy_(reference.sieve.node_logits)
        return mdl

    # ------------------------------------------------------------------
    def n_params(self) -> dict:
        node = sum(s.node_logits.numel() for s in self.sieves)
        gate = sum(p.numel() for p in self.gate.parameters())
        return {"node_params": node, "gate_params": gate, "total": node + gate}


# ----------------------------------------------------------------------
# collapse baselines expressed as parameter restrictions of CSMO
# ----------------------------------------------------------------------
def constant_p_model(cfg, seed=0, device=None, dtype=DTYPE):
    """Pooled supervised regression: F(z) = const. CSMO with one 0-dim view."""
    mdl = CSMO(cfg, views=[()], head="sieve", seed=seed, device=device, dtype=dtype)
    # freeze gate to the constant 1 (A=1 -> softmax over one logit is exactly 1)
    for p in mdl.gate.parameters():
        p.requires_grad_(False)
    return mdl


def support_mean_calibration_model(cfg, seed=0, device=None, dtype=DTYPE):
    """F depends on the support only through z13 (mean support label)."""
    return CSMO(cfg, views=[(12,)], head="sieve", gate_input_mask=[12], seed=seed,
                device=device, dtype=dtype)


def query_only_model(cfg, seed=0, device=None, dtype=DTYPE):
    """F sees the query/declared blocks only. Support-blind by construction."""
    return CSMO(cfg, views=[(0, 1), (2, 3), (26, 27)], head="sieve",
                gate_input_mask=[0, 1, 2, 3, 26, 27], seed=seed,
                device=device, dtype=dtype)


def deep_direct_model(cfg, seed=0, device=None, dtype=DTYPE):
    """Candidate B: deep simplex regressor, no sieve bank."""
    return CSMO(cfg, head="direct", seed=seed, device=device, dtype=dtype)


def scale_only_model(cfg, seed=0, device=None, dtype=DTYPE):
    """F depends only on support-label spread z14 (scale-only correction)."""
    return CSMO(cfg, views=[(13,)], head="sieve", gate_input_mask=[13], seed=seed,
                device=device, dtype=dtype)


def uniform_mixture_model(cfg, seed=0, views=None, device=None, dtype=DTYPE):
    """Same sieve bank as CSMO with an exactly uniform frozen gate."""
    return CSMO(cfg, views=views, seed=seed, gate_mode="uniform",
                device=device, dtype=dtype)


def frozen_random_gate_model(cfg, seed=0, views=None, device=None, dtype=DTYPE):
    """Same sieve bank as CSMO with a fixed seeded nonuniform routing vector."""
    device = require_cuda(device)
    n_views = len(DEFAULT_VIEWS[: cfg.n_views] if views is None else views)
    gen = torch.Generator(device=device).manual_seed(seed + 1009)
    weights = torch.rand(n_views, generator=gen, dtype=dtype, device=device)
    return CSMO(cfg, views=views, seed=seed, gate_mode="fixed",
                fixed_gate_weights=weights, device=device, dtype=dtype)


def capacity_matched_single_view_model(cfg, view=DEFAULT_VIEWS[0], seed=0,
                                       device=None, dtype=DTYPE):
    """Single-view sieve with node count matched to the default full CSMO."""
    target = CSMO(cfg, seed=seed, device=device, dtype=dtype).n_params()["total"]
    res = max(1, int(round((target / cfg.n_coef) ** 0.5)) - 1)
    return CSMO(cfg, views=[view], seed=seed, res=res, device=device, dtype=dtype)


# Source section: bandoperator
"""The frozen deployment object B(z) and the set-valued output K(B(z)p).

Refactored from the audited implementation: the deployment no longer owns the
statistic (that is now the biological frontend's job); it owns only the frozen
band machinery. The mathematics is unchanged.

FROZEN (theorem-exact):
    B(z) = [ b_pop_{kappa(z)} | beta_1 | ... | beta_m ], anchors z-INDEPENDENT,
    every column in the valid band polytope, assembly p -> B(z)p linear and
    convex-closed, kappa_B = sup_z ||B(z)||_op finite, and
    A(F,z) = K(B(z)F(z)) is the SOLE operator.
"""

import numpy as np
import torch
import torch.nn as nn

from . import bands
from .runtime import require_cuda, to_numpy

DEFAULT_DTYPE = torch.float64

# fixed anchor lattice (data-free). Centres are spread across V=[0,1] so the
# midband readout can resolve the full affinity range; the previous 3-centre
# lattice (0.20/0.50/0.80) pinned every prediction to ~3 levels regardless of z.
# Six sharp logistic CDFs at evenly spaced centres plus one broad uniform anchor.
ANCHOR_SPEC = [
    ("logistic", 0.15, 0.08, 0.12),
    ("logistic", 0.30, 0.08, 0.12),
    ("logistic", 0.45, 0.08, 0.12),
    ("logistic", 0.60, 0.08, 0.12),
    ("logistic", 0.75, 0.08, 0.12),
    ("logistic", 0.90, 0.08, 0.12),
    ("uniform", 0.0, 0.0, 0.90),
]


def build_anchors(cfg, device=None, dtype=DEFAULT_DTYPE) -> torch.Tensor:
    grid = cfg.grid()
    out = []
    for kind, c, s, w in ANCHOR_SPEC[: cfg.m]:
        u = (bands.logistic_shape(grid, c, s) if kind == "logistic"
             else np.linspace(0.0, 1.0, cfg.n_grid))
        beta = bands.band_from_shape(u, w)
        bands.assert_valid(beta, name=f"anchor({kind},{c},{s},{w})")
        out.append(beta)
    return torch.as_tensor(np.stack(out, 0), dtype=dtype, device=device)


def build_population_bands(cfg, ctx_src, Y_src, device=None,
                           dtype=DEFAULT_DTYPE) -> torch.Tensor:
    """b_pop_c from SOURCE labels only."""
    grid = cfg.grid()
    ctx = to_numpy(ctx_src) if torch.is_tensor(ctx_src) else np.asarray(ctx_src)
    Y = to_numpy(Y_src) if torch.is_tensor(Y_src) else np.asarray(Y_src)
    n_context_cont = len(cfg.kappa_edges_context_cont) + 1

    def band_from_samples(samples):
        ecdf = (samples[:, None] <= grid[None, :]).mean(axis=0)
        eps = max(np.sqrt(np.log(2.0 / cfg.dkw_alpha) / (2.0 * len(samples))),
                  cfg.dkw_eps_min)
        return bands.band_from_ecdf(ecdf, eps)

    def data_free_band():
        # Deliberately broad when an endpoint has no sufficient source evidence.
        return bands.join(np.zeros(cfg.n_grid), np.ones(cfg.n_grid))

    out = np.zeros((cfg.n_context, cfg.band_dim))
    for c in range(cfg.n_context):
        sel = Y[ctx == c]
        context_cont_bin = c % n_context_cont
        endpoint = Y[(ctx % n_context_cont) == context_cont_bin]
        if len(sel) >= cfg.dkw_n_min:
            beta = band_from_samples(sel)
        elif len(endpoint) >= cfg.dkw_n_min:
            # Sparse joint context: borrow only from the same declared endpoint.
            beta = band_from_samples(endpoint)
        else:
            beta = data_free_band()
        bands.assert_valid(beta, name=f"b_pop[{c}]")
        out[c] = beta
    return torch.as_tensor(out, dtype=dtype, device=device)


class BandOperator(nn.Module):
    """Frozen B(.) plus the assembly and the emitted law class."""

    def __init__(self, cfg, anchors: torch.Tensor, b_pop: torch.Tensor,
                 device=None, dtype=DEFAULT_DTYPE):
        super().__init__()
        device = require_cuda(device)
        self.cfg = cfg
        Bt = torch.zeros(cfg.n_context, cfg.band_dim, cfg.n_coef,
                         dtype=dtype, device=device)
        Bt[:, :, 0] = b_pop.to(device=device, dtype=dtype)
        Bt[:, :, 1:] = anchors.to(device=device, dtype=dtype).T.unsqueeze(0)
        self.register_buffer("B_table", Bt)
        self.register_buffer("grid_t", torch.as_tensor(cfg.grid(), dtype=dtype,
                                                       device=device))
        for b in self.buffers():
            b.requires_grad_(False)
        self._audit_columns()

    # -- audits ----------------------------------------------------------
    def _audit_columns(self):
        Bt = to_numpy(self.B_table)
        for c in range(self.cfg.n_context):
            for k in range(self.cfg.n_coef):
                bands.assert_valid(Bt[c, :, k], name=f"B_table[{c}][:,{k}]")
        sep = self.column_separation()
        if sep < 0.05:
            raise AssertionError(f"degenerate deployment: column separation {sep:.4f}")

    def column_separation(self) -> float:
        Bt = to_numpy(self.B_table)
        worst = np.inf
        for c in range(self.cfg.n_context):
            cols = Bt[c].T
            for i in range(len(cols)):
                for j in range(i + 1, len(cols)):
                    worst = min(worst, float(np.max(np.abs(cols[i] - cols[j]))))
        return worst

    def assembly_norm(self) -> float:
        """kappa_B = sup_z ||B(z)||_op, EXACT (finite context alphabet)."""
        return float(np.max(np.linalg.norm(to_numpy(self.B_table), axis=2)))

    # -- the rule --------------------------------------------------------
    def context(self, z: torch.Tensor) -> torch.Tensor:
        return context_index(z, self.cfg)

    def B_of_z(self, z: torch.Tensor) -> torch.Tensor:
        return self.B_table[self.context(z)]

    def assemble(self, z: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
        """beta = B(z) p. No clipping: validity holds by convexity."""
        return bands.assemble(self.B_of_z(z), p)

    def constants(self) -> dict:
        cfg = self.cfg
        kB = self.assembly_norm()
        L_base = cfg.L_lip * kB
        return {
            "h": cfg.h, "mu": cfg.mu, "m": cfg.m, "D_V_val": cfg.D_V_val,
            "kappa_B": kB, "D_V": cfg.D_V_val * kB,
            "L_bar": cfg.L_bar, "L_lip": cfg.L_lip,
            "L_base_upper": L_base, "diam_Delta_m": cfg.diam_simplex,
            "L_p_star_upper": L_base + cfg.mu * cfg.diam_simplex,
            "Lambda": cfg.L_lip * kB + cfg.mu * cfg.diam_simplex,
            "two_h_floor": 2.0 * cfg.h, "n_context": cfg.n_context,
            "column_separation": self.column_separation(),
        }


def build_band_operator(cfg, z_source: torch.Tensor, Y_source: torch.Tensor,
                        device=None, dtype=DEFAULT_DTYPE) -> BandOperator:
    """Stage 4. SOURCE statistics and SOURCE labels only."""
    device = require_cuda(device)
    with torch.no_grad():
        ctx = context_index(z_source.to(device), cfg)
    anchors = build_anchors(cfg, device=device, dtype=dtype)
    b_pop = build_population_bands(cfg, ctx, Y_source, device=device, dtype=dtype)
    return BandOperator(cfg, anchors, b_pop, device=device, dtype=dtype)


def _sha256_json(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _band_table_sha256(B_table) -> str:
    table = to_numpy(B_table) if torch.is_tensor(B_table) else np.asarray(B_table)
    digest = hashlib.sha256()
    digest.update(str(table.dtype).encode("ascii"))
    digest.update(np.asarray(table.shape, dtype=np.int64).tobytes())
    digest.update(np.ascontiguousarray(table).tobytes())
    return digest.hexdigest()


DEPLOYMENT_ARTIFACT_KEYS = (
    "state_schema_hash", "view_registry_hash", "context_registry_hash",
    "mechanism_schema_hash", "protein_bank_hash", "ligand_bank_hash",
    "pair_bank_hash", "archetype_manifest_hash",
)


def deployment_manifest(cfg: DeploymentConfig, B_table, *, frontend_hash: str,
                        source_manifest_hash: str, artifact_hashes: dict | None = None) -> dict:
    """Bind a frozen B(z) artifact to its statistic and source-data semantics."""
    artifacts = {} if artifact_hashes is None else dict(artifact_hashes)
    unknown = set(artifacts) - set(DEPLOYMENT_ARTIFACT_KEYS)
    if unknown:
        raise ValueError(f"unknown deployment artifact hashes: {sorted(unknown)}")
    missing = sorted(key for key in DEPLOYMENT_ARTIFACT_KEYS if not artifacts.get(key))
    if missing:
        raise ValueError(f"deployment artifact hashes are incomplete: {missing}")
    if not frontend_hash or not source_manifest_hash:
        raise ValueError("frontend and source manifest hashes must be nonempty")
    manifest = {
        "version": 2,
        "math_config_hash": cfg.fingerprint(),
        "state_schema_hash": artifacts["state_schema_hash"],
        "view_registry_hash": artifacts["view_registry_hash"],
        "context_registry_hash": artifacts["context_registry_hash"],
        "mechanism_schema_hash": artifacts["mechanism_schema_hash"],
        "protein_bank_hash": artifacts["protein_bank_hash"],
        "ligand_bank_hash": artifacts["ligand_bank_hash"],
        "pair_bank_hash": artifacts["pair_bank_hash"],
        "frontend_hash": str(frontend_hash),
        "source_manifest_hash": str(source_manifest_hash),
        "B_table_hash": _band_table_sha256(B_table),
        "archetype_manifest_hash": artifacts["archetype_manifest_hash"],
    }
    manifest["deployment_hash"] = _sha256_json(manifest)
    return manifest


def validate_deployment_manifest(manifest: dict, cfg: DeploymentConfig, B_table, *,
                                 frontend_hash: str, source_manifest_hash: str,
                                 artifact_hashes: dict | None = None) -> None:
    """Fail closed when a deployment artifact predates the current z semantics."""
    expected = deployment_manifest(
        cfg, B_table, frontend_hash=frontend_hash,
        source_manifest_hash=source_manifest_hash,
        artifact_hashes=artifact_hashes,
    )
    missing = sorted(set(expected) - set(manifest))
    if missing:
        raise ValueError(f"deployment manifest is incomplete: {missing}")
    mismatched = [key for key, value in expected.items() if manifest.get(key) != value]
    if mismatched:
        raise ValueError(f"deployment artifact contract mismatch: {mismatched}")
