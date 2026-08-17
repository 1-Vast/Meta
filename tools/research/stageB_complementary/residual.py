"""Residual decomposition for complementary partial meta-adaptation.

The Stage A inner loop and the incumbent Tanimoto transport were fitted to the
*same* support residual, so they explained the same information twice. That is
the measured reason Stage A's gain shrank as support grew instead of growing.
This module splits the support signal into three disjoint parts so a
meta-adapter can be trained on only the part nothing else explains:

    support residual  r_i = y_i - zero_shot_i
    (1) target level        L        = shrink * mean(r)
    (2) neighbourhood       t_i      = leave-one-out Tanimoto transport of r - L
    (3) complementary       c_i      = r_i - L - t_i        <- the adapter's target

The leave-one-out rule in (2) is load-bearing. A support item's own label is its
own nearest neighbour at Tanimoto 1.0, so an ordinary support-to-support
transport would predict `r_i` from `r_i` and drive `c_i` to zero — the adapter
would then be trained on numerical noise and would look inert for a reason that
has nothing to do with chemistry.

Conventions here match the incumbent: `shrink = n / (n + softplus(log_shrinkage))`
and transport weights are `softmax(similarity_scale * tanimoto)`.
"""
from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import Tensor

from model.similarity_grammar import tanimoto


def level_term(support_residual: Tensor, shrink: Tensor) -> Tensor:
    """One scalar per episode: the shrunken mean support residual. `[B, 1]`."""
    if support_residual.shape[-1] == 0:
        return support_residual.new_zeros(support_residual.shape[0], 1)
    return shrink * support_residual.mean(-1, keepdim=True)


def transport_weights(query_fingerprint: Tensor, support_fingerprint: Tensor,
                      scale: Tensor) -> Tensor:
    """`softmax_k(scale * Tanimoto(q, k))`, the incumbent's weighting. `[B,Q,K]`."""
    return torch.softmax(scale * tanimoto(query_fingerprint,
                                          support_fingerprint), -1)


def query_transport(query_fingerprint: Tensor, support_fingerprint: Tensor,
                    support_value: Tensor, scale: Tensor) -> Tensor:
    """Ordinary query-to-support transport of `support_value`. `[B, Q]`."""
    if support_value.shape[-1] == 0:
        return support_value.new_zeros(query_fingerprint.shape[0],
                                       query_fingerprint.shape[1])
    weight = transport_weights(query_fingerprint, support_fingerprint, scale)
    return torch.einsum("bqk,bk->bq", weight, support_value)


def leave_one_out_transport(support_fingerprint: Tensor, support_value: Tensor,
                            scale: Tensor) -> Tensor:
    """Predict each support item from the *others* only. `[B, K]`.

    Masking the diagonal before the softmax is what makes this leave-one-out:
    a support item is its own Tanimoto-1.0 neighbour, so including it would let
    the transport reproduce `r_i` from `r_i` and leave nothing complementary to
    learn. With `K == 1` there is no other item and the result is exactly zero.
    """
    count = support_value.shape[-1]
    if count <= 1:
        return torch.zeros_like(support_value)
    similarity = tanimoto(support_fingerprint, support_fingerprint)
    logits = scale * similarity
    eye = torch.eye(count, device=logits.device, dtype=torch.bool)
    logits = logits.masked_fill(eye, torch.finfo(logits.dtype).min)
    weight = torch.softmax(logits, -1)
    return torch.einsum("bjk,bk->bj", weight, support_value)


def decompose(support_y: Tensor, support_zero: Tensor,
              support_fingerprint: Tensor, shrink: Tensor,
              scale: Tensor) -> dict[str, Tensor]:
    """Split the support residual into level, neighbourhood and complement.

    Centering uses the **raw** mean, not the shrunken one. Shrinkage is a
    confidence discount applied when the correction is emitted; folding it into
    the decomposition would leave `centered = r - shrink*mean(r)`, which is not
    mean-zero and — at k = 1 — would leave a spurious `(1 - shrink) * r` for the
    adapter to fit. One support label cannot identify shape, so that residue
    would be pure noise dressed as signal.
    """
    residual = support_y - support_zero
    raw_level = residual.mean(-1, keepdim=True) if residual.shape[-1] else \
        residual.new_zeros(residual.shape[0], 1)
    centered = residual - raw_level
    neighbourhood = leave_one_out_transport(
        support_fingerprint, centered, scale)
    complementary = centered - neighbourhood
    return {"residual": residual, "raw_level": raw_level,
            "level": shrink * raw_level, "centered": centered,
            "neighbourhood": neighbourhood, "complementary": complementary}


def inner_target(mode: str, support_y: Tensor, support_zero: Tensor,
                 parts: dict[str, Tensor]) -> Tensor:
    """What the inner loop is asked to fit, per arm.

    ``H`` (naive hybrid, Stage A's behaviour) and ``M`` (meta only) fit the raw
    label, so they re-explain the level and the neighbourhood that the transport
    also explains. ``C`` fits the zero-shot endpoint plus only the complementary
    part, so each mechanism owns a disjoint share of the support signal.
    """
    if mode in {"H", "M"}:
        return support_y
    if mode == "C":
        return support_zero + parts["complementary"]
    raise ValueError(f"no inner target for arm {mode!r}")


def centered_shape(correction: Tensor) -> Tensor:
    """The part of a query correction that reorders queries. `[B, Q]`.

    Stage A called the weight-only update "pure shape". That was wrong: the
    weight update moves a query by `-2*lr*r_s*(h_q . h_s)`, whose within-episode
    mean is not zero, so it carries a level component too. Only the correction
    with its own episode mean removed is shape.
    """
    if correction.shape[-1] < 2:
        return torch.zeros_like(correction)
    return correction - correction.mean(-1, keepdim=True)


def conditioning_alpha(hidden: Tensor, inner_lr: float,
                       adapt_bias: bool = True) -> Tensor:
    """`alpha = 2 * lr * (||h||^2 + 1)` when the bias is adapted too.

    Stage A reported `2 * lr * ||h||^2`, which omits the bias. The bias gradient
    of a squared error is `2 * (p - y) * 1`, contributing exactly 1 to the same
    sum, so the published Stage A alpha was an underestimate.

    This governs the **support** residual only. It does not predict the query
    error: one step moves a query by

        delta_f(q) = -2 * lr * r_s * (h_q . h_s + 1)

    which depends on the query-support inner product and therefore varies across
    the panel. Support contraction and query MSE are different quantities and
    must not be conflated.
    """
    bias = 1.0 if adapt_bias else 0.0
    return 2.0 * inner_lr * (hidden.square().sum(-1) + bias)


def query_step_delta(hidden_query: Tensor, hidden_support: Tensor,
                     support_residual: Tensor, inner_lr: float,
                     adapt_bias: bool = True) -> Tensor:
    """The exact first inner-step effect on each query, at k = 1. `[B, Q]`.

    `delta_f(q) = -2 * lr * r_s * (h_q . h_s + bias)`. Exposed so the report can
    state the query effect directly instead of inferring it from the support
    contraction.
    """
    bias = 1.0 if adapt_bias else 0.0
    inner = torch.einsum("bqh,bkh->bqk", hidden_query, hidden_support) + bias
    return -2.0 * inner_lr * torch.einsum(
        "bqk,bk->bq", inner, support_residual)
