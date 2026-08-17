"""The exact episodic A2 moment operator, as `NEXT_RESEARCH_PLAN_A2_MOMENT` §2.

The v2 representation probe trained a *zero-shot bilinear delta predictor*

    delta_hat(P, L_i, L_j) = s · ⟨ g(P), U (e_i − e_j) ⟩

and used its failure to reject A2. That was an over-reach. The predictor shares
A2's feature space but not its structure: it never forms a support moment, it
reads a ligand *pair* rather than a query against an episode, and it has no
shrinkage in k. A2's actual claim is that **support residuals identify a
direction in a learned coordinate system**, and only the operator below tests
it.

    z_i     = A_φ(e0(P, L_i))                       # learned projection, rank R
    r_i     = stopgrad(y_i − f0(P, L_i))            # label-locked residual
    c_S     = (1/k) Σ_i r_i z_i                     # first moment, order-free
    η(k)    = η_∞ · k / (k + λ)                     # shrinkage, exactly 0 at k=0
    δ_q     = η(k) · ⟨c_S, z_q⟩
    f       = f0 + δ

Trainable: `A_φ` (a single `Linear(D, R, bias=False)`), `η_∞` and `λ`. Nothing
else. No solver, no pseudoinverse, no inner loop, no test-time gradient, no
query label anywhere in the forward path.

Five structural properties follow from the algebra and are tested rather than
asserted (`tests/test_operator_contract.py`):

1. **exact k=0 identity** — `η(0) = 0`, so `f ≡ f0` bit-exactly;
2. **non-scalar k=1** — unlike A0, whose k=1 transport is provably a pure level
   shift (`sar_adaptation ≡ 0`), `δ_q = η(1)·r_1·⟨z_1, z_q⟩` varies with the
   query through `⟨z_1, z_q⟩`. This is the one structural advantage A2 has over
   the incumbent and it must be demonstrated, not assumed;
3. **support permutation invariance** — `c_S` is a mean;
4. **query permutation equivariance** — queries never interact;
5. **no query-label path** — `query_y` appears in no signature here.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class A2MomentOperator(nn.Module):
    """`f = f0 + η(k)·⟨mean_i r_i z_i, z_q⟩` over a frozen representation."""

    def __init__(self, width: int, rank: int, *, learn_projection: bool = True,
                 eta_init: float = 1.0, lam_init: float = 1.0) -> None:
        super().__init__()
        self.width, self.rank = int(width), int(rank)
        self.projection = nn.Linear(width, rank, bias=False)
        if not learn_projection:
            # The random-projection control: the coordinate system is frozen at
            # initialisation and only the two shrinkage scalars train.
            self.projection.weight.requires_grad_(False)
        # Softplus-parameterised so both stay positive without a constraint.
        self.log_eta = nn.Parameter(torch.tensor(float(eta_init)))
        self.log_lambda = nn.Parameter(torch.tensor(float(lam_init)))

    def shrinkage(self, support_count: int, reference: torch.Tensor) -> torch.Tensor:
        """`η(k) = η_∞ k/(k+λ)`. Exactly zero at k=0, for any parameter value."""
        count = reference.new_tensor(float(support_count))
        lam = F.softplus(self.log_lambda)
        return F.softplus(self.log_eta) * count / (count + lam)

    def coordinates(self, features: torch.Tensor) -> torch.Tensor:
        return self.projection(features)

    def forward(self, support_features: torch.Tensor,
                support_residual: torch.Tensor,
                query_features: torch.Tensor) -> torch.Tensor:
        """Return `δ_q`, the correction added to the frozen zero-shot endpoint.

        `support_features`  [k, D]      frozen `e0(P, L_i)`
        `support_residual`  [k]         `y_i − f0(P, L_i)`, already detached
        `query_features`    [Q, D]      frozen `e0(P, L_q)`
        """
        support_count = int(support_features.shape[0])
        if support_count == 0:
            # Not merely η(0)=0 — return a literal zero so the k=0 identity
            # cannot depend on floating-point luck in an empty reduction.
            return query_features.new_zeros(query_features.shape[0])
        z_support = self.coordinates(support_features)               # [k, R]
        z_query = self.coordinates(query_features)                   # [Q, R]
        moment = (support_residual.unsqueeze(-1) * z_support).mean(0)  # [R]
        return self.shrinkage(support_count, z_query) * (z_query @ moment)


class ScalarLevelOperator(nn.Module):
    """The baseline A2 must beat: the shrunken support-residual mean.

    `δ_q = η(k)·mean_i r_i` — identical for every query. This is what A0's
    transport degenerates to at k=1, and any "adaptation" that does not beat it
    is a level shift wearing a coordinate system.
    """

    def __init__(self, eta_init: float = 1.0, lam_init: float = 1.0) -> None:
        super().__init__()
        self.log_eta = nn.Parameter(torch.tensor(float(eta_init)))
        self.log_lambda = nn.Parameter(torch.tensor(float(lam_init)))

    def shrinkage(self, support_count: int, reference: torch.Tensor) -> torch.Tensor:
        count = reference.new_tensor(float(support_count))
        return (F.softplus(self.log_eta) * count
                / (count + F.softplus(self.log_lambda)))

    def forward(self, support_features: torch.Tensor,
                support_residual: torch.Tensor,
                query_features: torch.Tensor) -> torch.Tensor:
        if support_features.shape[0] == 0:
            return query_features.new_zeros(query_features.shape[0])
        level = self.shrinkage(int(support_features.shape[0]), query_features)
        return (level * support_residual.mean()).expand(query_features.shape[0])


class SharedMomentOperator(nn.Module):
    """The protein-independent control: one learned direction, not a moment.

    `δ_q = η(k) · (mean_i r_i) · ⟨c, z_q⟩`

    Identical machinery to `A2MomentOperator` — same projection, same
    shrinkage, same query-specific inner product — except that the direction
    `c` is a **learned constant** instead of being formed from the support
    features. It therefore uses the support *labels* (through their mean) but
    not the support *chemistry*.

    If this matches A2, then A2's moment is not identifying anything episode-
    specific: a fixed direction scaled by the support mean reproduces it, and
    the "coordinate system identified by k≤5 labels" reduces to a query-
    dependent rescaling of a level shift.
    """

    def __init__(self, width: int, rank: int, eta_init: float = 1.0,
                 lam_init: float = 1.0) -> None:
        super().__init__()
        self.projection = nn.Linear(width, rank, bias=False)
        self.direction = nn.Parameter(torch.randn(rank) / rank ** 0.5)
        self.log_eta = nn.Parameter(torch.tensor(float(eta_init)))
        self.log_lambda = nn.Parameter(torch.tensor(float(lam_init)))

    def shrinkage(self, support_count: int, reference: torch.Tensor) -> torch.Tensor:
        count = reference.new_tensor(float(support_count))
        return (F.softplus(self.log_eta) * count
                / (count + F.softplus(self.log_lambda)))

    def forward(self, support_features: torch.Tensor,
                support_residual: torch.Tensor,
                query_features: torch.Tensor) -> torch.Tensor:
        if support_features.shape[0] == 0:
            return query_features.new_zeros(query_features.shape[0])
        z_query = self.projection(query_features)
        scale = self.shrinkage(int(support_features.shape[0]), z_query)
        return scale * support_residual.mean() * (z_query @ self.direction)


def tanimoto_transport(support_fingerprint: torch.Tensor,
                       query_fingerprint: torch.Tensor,
                       support_residual: torch.Tensor,
                       similarity_scale: float = 8.0,
                       shrinkage: float = 1.8546) -> torch.Tensor:
    """The incumbent's fixed Morgan/Tanimoto residual transport, parameter-free.

    Reproduces `model.similarity_grammar.SimilarityTransport` at A0's recorded
    initialisation (`similarity_scale=8.0`, `log_shrinkage=1.8546`), which the
    R3R4 record shows barely moves during training. This is the comparator A2
    has to beat to be worth anything: it already exploits ligand-side SAR
    continuity and it has no learned coordinate system at all.
    """
    support_count = support_fingerprint.shape[0]
    if support_count == 0:
        return query_fingerprint.new_zeros(query_fingerprint.shape[0])
    intersection = query_fingerprint @ support_fingerprint.T
    union = (query_fingerprint.sum(-1)[:, None]
             + support_fingerprint.sum(-1)[None, :] - intersection)
    similarity = intersection / union.clamp_min(1e-6)
    weight = torch.softmax(similarity_scale * similarity, -1)
    strength = F.softplus(torch.tensor(shrinkage, dtype=weight.dtype))
    shrink = support_count / (support_count + strength)
    return shrink * (weight @ support_residual)
