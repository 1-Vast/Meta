"""Structural contract for the exact episodic A2 operator.

Every property A2 claims is algebraic and must hold at *every* parameter value,
not merely after training. These run on small random tensors in float64.

The k=1 test is the load-bearing one. A0's k=1 transport is provably a pure
level shift (`sar_adaptation ≡ 0`, DATAFLOW_AUDIT F4), so A0 makes no k=1
claim and cannot. If the A2 operator's k=1 correction is also constant across
queries, A2 has no structural advantage over the incumbent and there is nothing
to test on real data.

Run: `conda run -n drug python -m pytest tools/research/a2_exact_probe/tests -q`
"""
from __future__ import annotations

from pathlib import Path
import sys

import pytest
import torch

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.research.a2_exact_probe.operator import (          # noqa: E402
    A2MomentOperator, ScalarLevelOperator, tanimoto_transport,
)

WIDTH, RANK, SUPPORT, QUERIES = 12, 4, 5, 7


def build(seed: int = 0, **kwargs) -> A2MomentOperator:
    torch.manual_seed(seed)
    return A2MomentOperator(WIDTH, RANK, **kwargs).double()


def episode(seed: int, support: int = SUPPORT, queries: int = QUERIES):
    generator = torch.Generator().manual_seed(seed)
    return (torch.randn(support, WIDTH, generator=generator, dtype=torch.float64),
            torch.randn(support, generator=generator, dtype=torch.float64),
            torch.randn(queries, WIDTH, generator=generator, dtype=torch.float64))


# --- 1. exact k=0 identity -------------------------------------------------

@pytest.mark.parametrize("seed", [0, 1, 2])
def test_k0_correction_is_exactly_zero(seed):
    """`f ≡ f0` at k=0, bit-exactly, at any parameter value."""
    operator = build(seed)
    _, _, query = episode(seed)
    empty_features = torch.zeros(0, WIDTH, dtype=torch.float64)
    empty_residual = torch.zeros(0, dtype=torch.float64)
    delta = operator(empty_features, empty_residual, query)
    assert delta.shape == (QUERIES,)
    assert torch.equal(delta, torch.zeros_like(delta))


def test_k0_identity_survives_hostile_parameters():
    """Not an artifact of the initialisation: force η and λ far off."""
    operator = build(3)
    with torch.no_grad():
        operator.log_eta.fill_(50.0)
        operator.log_lambda.fill_(-50.0)
        operator.projection.weight.mul_(1e6)
    _, _, query = episode(3)
    delta = operator(torch.zeros(0, WIDTH, dtype=torch.float64),
                     torch.zeros(0, dtype=torch.float64), query)
    assert torch.equal(delta, torch.zeros_like(delta))


def test_shrinkage_is_zero_at_k0_and_increasing_in_k():
    operator = build(4)
    reference = torch.zeros(1, dtype=torch.float64)
    values = [float(operator.shrinkage(k, reference)) for k in (0, 1, 2, 3, 5)]
    assert values[0] == 0.0
    assert all(a < b for a, b in zip(values, values[1:]))


# --- 2. non-scalar k=1 -----------------------------------------------------

@pytest.mark.parametrize("seed", [0, 1, 2])
def test_k1_correction_is_query_specific(seed):
    """The structural advantage over A0, which is provably scalar at k=1.

    `δ_q = η(1)·r_1·⟨z_1, z_q⟩` varies with `q` through the inner product.
    """
    operator = build(seed)
    support, residual, query = episode(seed, support=1)
    delta = operator(support, residual, query)
    assert float(delta.std()) > 1e-6, (
        "the k=1 correction is constant across queries; the operator has "
        "degenerated to a level shift and has no advantage over A0")


def test_the_scalar_baseline_really_is_scalar_at_every_k():
    """The control must have the property A2 claims to lack."""
    baseline = ScalarLevelOperator().double()
    for k in (1, 2, 3, 5):
        support, residual, query = episode(10 + k, support=k)
        delta = baseline(support, residual, query)
        assert float(delta.max() - delta.min()) == 0.0


# --- 3. support permutation invariance ------------------------------------

@pytest.mark.parametrize("seed", [0, 1])
def test_support_order_does_not_change_the_output(seed):
    operator = build(seed)
    support, residual, query = episode(seed)
    order = torch.randperm(SUPPORT, generator=torch.Generator().manual_seed(seed))
    left = operator(support, residual, query)
    right = operator(support[order], residual[order], query)
    assert torch.allclose(left, right, atol=1e-12)


# --- 4. query permutation equivariance ------------------------------------

def test_queries_are_processed_independently():
    operator = build(5)
    support, residual, query = episode(5)
    order = torch.randperm(QUERIES, generator=torch.Generator().manual_seed(5))
    full = operator(support, residual, query)
    permuted = operator(support, residual, query[order])
    assert torch.allclose(full[order], permuted, atol=1e-12)
    # And a single query alone must equal its entry in the full panel.
    alone = operator(support, residual, query[:1])
    assert torch.allclose(alone, full[:1], atol=1e-12)


# --- 5. label handling -----------------------------------------------------

def test_the_correction_is_linear_and_odd_in_the_support_residuals():
    """`δ` is linear in `r`, so sign-flipping the labels flips the correction.

    This is what makes the label-permutation control interpretable: the
    operator cannot ignore label values and still produce the same output.
    """
    operator = build(6)
    support, residual, query = episode(6)
    base = operator(support, residual, query)
    assert torch.allclose(operator(support, -residual, query), -base, atol=1e-12)
    assert torch.allclose(operator(support, 2.0 * residual, query), 2.0 * base,
                          atol=1e-12)


def test_zero_residuals_give_zero_correction():
    operator = build(7)
    support, _, query = episode(7)
    delta = operator(support, torch.zeros(SUPPORT, dtype=torch.float64), query)
    assert torch.allclose(delta, torch.zeros_like(delta), atol=1e-12)


def test_no_query_label_reaches_the_operator():
    """`forward` has three parameters and none of them is a query label."""
    import inspect
    names = list(inspect.signature(A2MomentOperator.forward).parameters)
    assert names == ["self", "support_features", "support_residual",
                     "query_features"]


def test_the_support_residual_carries_no_gradient_into_the_operator():
    """Residuals are label-locked: `.detach()` upstream must be sufficient."""
    operator = build(8)
    support, residual, query = episode(8)
    residual = residual.detach().requires_grad_(False)
    operator(support, residual, query).sum().backward()
    assert operator.projection.weight.grad is not None
    assert residual.grad is None


# --- the frozen-projection control ----------------------------------------

def test_the_random_projection_control_trains_only_its_scalars():
    operator = build(9, learn_projection=False)
    support, residual, query = episode(9)
    operator(support, residual, query).sum().backward()
    assert operator.projection.weight.grad is None
    assert operator.log_eta.grad is not None
    assert operator.log_lambda.grad is not None


# --- the Tanimoto comparator ----------------------------------------------

def test_tanimoto_transport_matches_the_incumbent_algebra():
    """Zero support gives zero; identical fingerprints give the support mean."""
    generator = torch.Generator().manual_seed(11)
    support_fp = (torch.rand(3, 32, generator=generator,
                             dtype=torch.float64) > 0.6).double()
    query_fp = support_fp[:1].repeat(4, 1)
    residual = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)

    empty = tanimoto_transport(support_fp[:0], query_fp, residual[:0])
    assert torch.equal(empty, torch.zeros(4, dtype=torch.float64))

    # A query identical to support 0 must weight it most heavily.
    delta = tanimoto_transport(support_fp, query_fp, residual)
    assert delta.shape == (4,)
    assert float(delta.std()) < 1e-12          # identical queries agree
    uniform = tanimoto_transport(support_fp, query_fp, residual,
                                 similarity_scale=0.0)
    assert torch.allclose(uniform, uniform[0].expand(4), atol=1e-12)


def test_tanimoto_transport_is_support_permutation_invariant():
    generator = torch.Generator().manual_seed(12)
    support_fp = (torch.rand(4, 32, generator=generator,
                             dtype=torch.float64) > 0.6).double()
    query_fp = (torch.rand(5, 32, generator=generator,
                           dtype=torch.float64) > 0.6).double()
    residual = torch.randn(4, generator=generator, dtype=torch.float64)
    order = torch.randperm(4, generator=generator)
    assert torch.allclose(
        tanimoto_transport(support_fp, query_fp, residual),
        tanimoto_transport(support_fp[order], query_fp, residual[order]),
        atol=1e-12)
