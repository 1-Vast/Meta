"""Regression tests for the eight Stage A corrections.

Each test names the defect it guards and, where the defect was a wrong formula
or a wrong anchor, asserts the corrected behaviour differs from the old one —
otherwise a "correction" that changed nothing would still pass.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.research.stageB_complementary.residual import (           # noqa: E402
    centered_shape, conditioning_alpha, decompose, inner_target,
    leave_one_out_transport, level_term, query_step_delta,
    query_transport, transport_weights,
)


def fingerprints(count: int, bits: int = 64, seed: int = 0) -> torch.Tensor:
    generator = torch.Generator().manual_seed(seed)
    return (torch.rand(1, count, bits, generator=generator,
                       dtype=torch.float64) > 0.7).double()


# --- correction 3: the k=1 conditioning formula ----------------------------

def test_alpha_includes_the_adapted_bias():
    """Stage A used `2*lr*||h||^2`, omitting the bias term."""
    hidden = torch.tensor([[[0.3, 0.4]]], dtype=torch.float64)   # ||h||^2 = 0.25
    lr = 0.1
    corrected = float(conditioning_alpha(hidden, lr, adapt_bias=True))
    old = float(2.0 * lr * hidden.square().sum(-1))
    assert corrected == pytest.approx(2.0 * lr * (0.25 + 1.0), abs=1e-12)
    assert old == pytest.approx(2.0 * lr * 0.25, abs=1e-12)
    assert corrected > old, "the correction must change the number"
    # The omitted term is exactly `2*lr`, independent of the activation.
    assert corrected - old == pytest.approx(2.0 * lr, abs=1e-12)


def test_alpha_without_bias_recovers_the_weight_only_case():
    hidden = torch.randn(1, 3, 5, generator=torch.Generator().manual_seed(1),
                         dtype=torch.float64)
    weight_only = conditioning_alpha(hidden, 0.1, adapt_bias=False)
    both = conditioning_alpha(hidden, 0.1, adapt_bias=True)
    assert torch.allclose(both - weight_only, torch.full_like(both, 0.2),
                          atol=1e-12)


def test_query_step_delta_is_not_the_support_contraction():
    """The correction that must replace 'support contraction predicts query MSE'.

    One inner step moves each query by `-2*lr*r_s*(h_q . h_s + 1)`, which varies
    across the panel. A single scalar alpha cannot describe it.
    """
    generator = torch.Generator().manual_seed(2)
    hidden_support = torch.randn(1, 1, 6, generator=generator, dtype=torch.float64)
    hidden_query = torch.randn(1, 5, 6, generator=generator, dtype=torch.float64)
    residual = torch.tensor([[0.7]], dtype=torch.float64)
    delta = query_step_delta(hidden_query, hidden_support, residual, 0.1)
    assert delta.shape == (1, 5)
    # Queries move by different amounts: the effect is not one scalar.
    assert float(delta.std()) > 1e-6
    # And it is not equal to the support contraction for any query.
    alpha = float(conditioning_alpha(hidden_support, 0.1))
    support_delta = -alpha * float(residual)
    assert not torch.allclose(delta, torch.full_like(delta, support_delta),
                              atol=1e-6)


# --- correction 4: shape is the centered correction -------------------------

def test_centered_shape_removes_the_level_component():
    correction = torch.tensor([[1.0, 2.0, 3.0, 4.0]], dtype=torch.float64)
    shape = centered_shape(correction)
    assert float(shape.mean()) == pytest.approx(0.0, abs=1e-12)
    assert torch.allclose(shape, correction - 2.5, atol=1e-12)


def test_a_weight_only_update_is_not_pure_shape():
    """The claim Stage A made, shown false on the exact algebra.

    A weight-only step moves query q by `-2*lr*r_s*(h_q . h_s)`. Its
    within-episode mean is `-2*lr*r_s*(mean_q h_q . h_s)`, which is zero only if
    the mean query activation happens to be orthogonal to the support
    activation. Generic activations are not.
    """
    generator = torch.Generator().manual_seed(3)
    hidden_support = torch.randn(1, 1, 8, generator=generator, dtype=torch.float64)
    hidden_query = torch.randn(1, 6, 8, generator=generator, dtype=torch.float64)
    residual = torch.tensor([[0.9]], dtype=torch.float64)
    weight_only = query_step_delta(hidden_query, hidden_support, residual,
                                   0.1, adapt_bias=False)
    assert abs(float(weight_only.mean())) > 1e-6, (
        "a weight-only update carries a level component")
    assert float(centered_shape(weight_only).abs().max()) > 1e-9


def test_a_bias_only_update_is_exactly_level():
    """The complementary half: bias-only really is pure level."""
    generator = torch.Generator().manual_seed(4)
    hidden_support = torch.zeros(1, 1, 8, dtype=torch.float64)
    hidden_query = torch.randn(1, 6, 8, generator=generator, dtype=torch.float64)
    residual = torch.tensor([[0.9]], dtype=torch.float64)
    bias_only = query_step_delta(hidden_query, hidden_support, residual, 0.1)
    assert float(bias_only.std()) == pytest.approx(0.0, abs=1e-12)
    assert float(centered_shape(bias_only).abs().max()) == pytest.approx(0.0, abs=1e-12)


# --- the leave-one-out rule -------------------------------------------------

def test_leave_one_out_never_uses_an_item_to_predict_itself():
    """Without diagonal masking the transport reproduces `r_i` from `r_i`."""
    support = fingerprints(4, seed=5)
    values = torch.tensor([[1.0, -2.0, 3.0, 0.5]], dtype=torch.float64)
    scale = torch.tensor(8.0, dtype=torch.float64)
    loo = leave_one_out_transport(support, values, scale)
    naive = query_transport(support, support, values, scale)
    assert not torch.allclose(loo, naive, atol=1e-6)
    # The naive version is dominated by the self-match, so it tracks the values.
    assert float(torch.corrcoef(torch.stack(
        (naive[0], values[0])))[0, 1]) > 0.9


def test_leave_one_out_is_exactly_zero_for_a_single_support():
    support = fingerprints(1, seed=6)
    values = torch.tensor([[2.5]], dtype=torch.float64)
    loo = leave_one_out_transport(support, values, torch.tensor(8.0, dtype=torch.float64))
    assert torch.equal(loo, torch.zeros_like(values))


def test_leave_one_out_weights_are_a_proper_distribution_off_diagonal():
    support = fingerprints(5, seed=7)
    values = torch.eye(5, dtype=torch.float64)[None, 0]
    scale = torch.tensor(8.0, dtype=torch.float64)
    # Predict a one-hot: the result must be the off-diagonal weight on item 0.
    loo = leave_one_out_transport(support, values, scale)
    assert float(loo[0, 0]) == pytest.approx(0.0, abs=1e-12), (
        "item 0 must not see its own value")
    assert (loo[0, 1:] > 0).all()


def test_the_decomposition_is_exact_and_disjoint():
    """level + neighbourhood + complementary must reconstruct the residual."""
    generator = torch.Generator().manual_seed(8)
    support_y = torch.randn(1, 5, generator=generator, dtype=torch.float64)
    support_zero = torch.randn(1, 5, generator=generator, dtype=torch.float64)
    support_fp = fingerprints(5, seed=9)
    shrink = torch.tensor(0.7, dtype=torch.float64)
    scale = torch.tensor(8.0, dtype=torch.float64)
    parts = decompose(support_y, support_zero, support_fp, shrink, scale)
    total = parts["raw_level"] + parts["neighbourhood"] + parts["complementary"]
    assert torch.allclose(total, parts["residual"], atol=1e-12)
    # Centering uses the raw mean, so the centered part is exactly mean-zero.
    assert float(parts["centered"].mean().abs()) < 1e-12


def test_the_complementary_residual_is_exactly_zero_at_k1():
    """One support label carries no shape; the decomposition must say so."""
    support_y = torch.tensor([[1.7]], dtype=torch.float64)
    support_zero = torch.tensor([[0.4]], dtype=torch.float64)
    parts = decompose(support_y, support_zero, fingerprints(1, seed=30),
                      torch.tensor(0.6, dtype=torch.float64),
                      torch.tensor(8.0, dtype=torch.float64))
    assert float(parts["centered"].abs().max()) == pytest.approx(0.0, abs=1e-12)
    assert float(parts["neighbourhood"].abs().max()) == pytest.approx(0.0, abs=1e-12)
    assert float(parts["complementary"].abs().max()) == pytest.approx(0.0, abs=1e-12)
    # The whole residual is level, and the reported level carries the shrinkage.
    assert float(parts["raw_level"]) == pytest.approx(1.3, abs=1e-12)
    assert float(parts["level"]) == pytest.approx(0.6 * 1.3, abs=1e-12)


def test_the_complementary_target_differs_from_the_raw_label():
    """`C` must not be a relabelling of `H`."""
    generator = torch.Generator().manual_seed(10)
    support_y = torch.randn(1, 5, generator=generator, dtype=torch.float64)
    support_zero = torch.randn(1, 5, generator=generator, dtype=torch.float64)
    parts = decompose(support_y, support_zero, fingerprints(5, seed=11),
                      torch.tensor(0.7, dtype=torch.float64),
                      torch.tensor(8.0, dtype=torch.float64))
    naive = inner_target("H", support_y, support_zero, parts)
    complementary = inner_target("C", support_y, support_zero, parts)
    assert torch.allclose(naive, support_y, atol=0.0)
    assert not torch.allclose(complementary, naive, atol=1e-6)
    # And the complementary target removes the level exactly.
    assert abs(float((complementary - support_zero).mean()
                     - parts["complementary"].mean())) < 1e-12


def test_level_term_is_zero_on_an_empty_support():
    empty = torch.zeros(1, 0, dtype=torch.float64)
    assert float(level_term(empty, torch.tensor(0.5, dtype=torch.float64))) == 0.0


def test_transport_weights_sum_to_one():
    weights = transport_weights(fingerprints(3, seed=12), fingerprints(4, seed=13),
                                torch.tensor(8.0, dtype=torch.float64))
    assert torch.allclose(weights.sum(-1), torch.ones(1, 3, dtype=torch.float64),
                          atol=1e-12)


# --- correction 1: the matched-wrong anchor --------------------------------

def matched_wrong(support_y: torch.Tensor, anchor: torch.Tensor) -> torch.Tensor:
    return support_y - 2.0 * (support_y - anchor)


def test_matched_wrong_anchored_to_pre_adaptation_is_arm_independent():
    """The defect: Stage A anchored A1 to post-adaptation and A0 to pre.

    Anchoring to the pre-adaptation support prediction makes the corruption
    magnitude a property of the episode and the shared initialization, so both
    arms receive an identically corrupted control.
    """
    support_y = torch.tensor([[1.0, 2.0, 3.0]], dtype=torch.float64)
    pre = torch.tensor([[0.5, 1.0, 2.0]], dtype=torch.float64)
    post = torch.tensor([[0.9, 1.9, 2.9]], dtype=torch.float64)   # adapted, closer
    from_pre = matched_wrong(support_y, pre)
    from_post = matched_wrong(support_y, post)
    assert not torch.allclose(from_pre, from_post, atol=1e-6)
    # Post-anchoring shrinks the corruption, which is what biased the control.
    assert float((from_post - support_y).abs().mean()) < \
        float((from_pre - support_y).abs().mean())


def test_matched_wrong_preserves_residual_magnitude_and_flips_its_sign():
    support_y = torch.tensor([[1.0, 2.0, 3.0]], dtype=torch.float64)
    pre = torch.tensor([[0.5, 1.0, 2.0]], dtype=torch.float64)
    corrupted = matched_wrong(support_y, pre)
    assert torch.allclose(corrupted - pre, -(support_y - pre), atol=1e-12)
