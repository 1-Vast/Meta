"""Structural guarantees of the A2S-TRACE transport family.

These are not measurements.  Each one is a property that follows from the
functional form, and the corresponding control in the mechanism report is called
*structural* only because this file passes.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from research.a2s import a2s_trace as trace


DEVICE = trace.DEVICE


def make_inputs(batch: int = 4, n_query: int = 7, k: int = 5, seed: int = 0):
    generator = torch.Generator(device="cpu").manual_seed(seed)
    query_bits = (torch.rand(batch, n_query, trace.MORGAN_BITS, generator=generator) < 0.05).float()
    support_bits = (torch.rand(batch, k, trace.MORGAN_BITS, generator=generator) < 0.05).float()
    # Guarantee a non-empty fingerprint so Tanimoto is well defined.
    query_bits[..., 0] = 1.0
    support_bits[..., 0] = 1.0
    query_desc = torch.randn(batch, n_query, 10, generator=generator)
    support_desc = torch.randn(batch, k, 10, generator=generator)
    protein = torch.randn(batch, 1280, generator=generator)
    residual = torch.randn(batch, k, generator=generator)
    return tuple(
        tensor.to(DEVICE)
        for tensor in (query_bits, support_bits, query_desc, support_desc, protein, residual)
    )


def trained_like(config: trace.TraceConfig, seed: int = 3) -> trace.Trace:
    """A model whose heads are away from initialisation, as after training."""

    torch.manual_seed(seed)
    model = trace.Trace(config).to(DEVICE)
    with torch.no_grad():
        for parameter in model.parameters():
            parameter.add_(0.05 * torch.randn(parameter.shape, device=parameter.device))
    return model


@pytest.fixture(scope="module")
def config() -> trace.TraceConfig:
    return trace.TraceConfig()


def test_residual_null_is_a_bitwise_no_op(config):
    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, _ = make_inputs()
    zero = torch.zeros(query_bits.shape[0], support_bits.shape[1], device=DEVICE)
    delta, _ = model(query_bits, support_bits, query_desc, support_desc, protein, zero)
    assert torch.equal(delta, torch.zeros_like(delta))


def test_weights_do_not_read_the_residual(config):
    """A derangement must leave every transport weight bit-identical."""

    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    permutation = torch.as_tensor([[1, 2, 3, 4, 0]] * residual.shape[0], device=DEVICE)
    deranged = torch.gather(residual, 1, permutation)
    _, first = model(query_bits, support_bits, query_desc, support_desc, protein, residual)
    _, second = model(query_bits, support_bits, query_desc, support_desc, protein, deranged)
    assert torch.equal(first["modulation"], second["modulation"])
    assert torch.equal(first["gate"], second["gate"])


def test_transport_is_bounded_by_an_observed_quantity(config):
    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    residual = 5.0 * residual
    delta, _ = model(query_bits, support_bits, query_desc, support_desc, protein, residual)
    bound = residual.abs().amax(-1, keepdim=True)
    assert bool((delta.abs() <= bound + 1e-6).all())


def test_transport_is_linear_in_the_residual(config):
    """Doubling the evidence doubles the action, up to the declared bound."""

    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    residual = 0.05 * residual  # far from the clamp
    single, _ = model(query_bits, support_bits, query_desc, support_desc, protein, residual)
    double, _ = model(query_bits, support_bits, query_desc, support_desc, protein, 2.0 * residual)
    assert torch.allclose(double, 2.0 * single, atol=1e-6)


def test_support_order_is_irrelevant(config):
    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    order = torch.as_tensor([4, 0, 3, 1, 2], device=DEVICE)
    delta, _ = model(query_bits, support_bits, query_desc, support_desc, protein, residual)
    permuted, _ = model(
        query_bits,
        support_bits[:, order],
        query_desc,
        support_desc[:, order],
        protein,
        residual[:, order],
    )
    assert torch.allclose(delta, permuted, atol=1e-5)


def test_queries_do_not_influence_each_other(config):
    """No candidate-set statistic enters, so a query subset is exactly stable."""

    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    full, _ = model(query_bits, support_bits, query_desc, support_desc, protein, residual)
    subset, _ = model(
        query_bits[:, :3], support_bits, query_desc[:, :3], support_desc, protein, residual
    )
    assert torch.allclose(full[:, :3], subset, atol=1e-6)


def test_nw_restriction_reproduces_the_analytic_smoother():
    config = trace.TraceConfig(
        weights="nw", whiten=False, global_scale=False, modulation=False, gate=False
    )
    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    delta, _ = model(query_bits, support_bits, query_desc, support_desc, protein, residual)
    reference = trace.analytic_delta(query_bits, support_bits, residual, estimator="nw")
    assert torch.allclose(delta, reference, atol=1e-6)


def test_krr_restriction_reproduces_the_analytic_smoother():
    config = trace.TraceConfig(
        weights="krr", whiten=True, global_scale=False, modulation=False, gate=False, ridge=0.03
    )
    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    residual = 0.05 * residual  # keep the declared bound inactive
    delta, _ = model(query_bits, support_bits, query_desc, support_desc, protein, residual)
    reference = trace.analytic_delta(
        query_bits, support_bits, residual, estimator="krr", ridge=0.03
    )
    assert torch.allclose(delta, reference, atol=1e-5)


def test_protein_zero_removes_the_protein_channel():
    config = trace.TraceConfig(protein="zero")
    model = trained_like(config)
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    delta, _ = model(query_bits, support_bits, query_desc, support_desc, protein, residual)
    shuffled, _ = model(
        query_bits, support_bits, query_desc, support_desc, protein.flip(0), residual
    )
    assert torch.equal(delta, shuffled)


def test_level_channel_is_rank_null():
    query_bits, support_bits, _, _, _, residual = make_inputs()
    delta = trace.analytic_delta(query_bits, support_bits, residual, estimator="level")
    assert torch.allclose(delta, delta[:, :1].expand_as(delta))


def test_derangement_has_no_fixed_point():
    rng = np.random.default_rng(0)
    for size in (2, 3, 5, 8):
        for _ in range(20):
            order = trace.derangement(size, rng)
            assert not np.any(order == np.arange(size))


def test_stratum_edges_are_exhaustive_and_ordered():
    values = np.asarray([0.0, 0.19, 0.2, 0.34, 0.35, 0.54, 0.55, 1.0])
    assigned = trace.stratum_of(values)
    assert list(assigned) == [
        "t00_20", "t00_20", "t20_35", "t20_35", "t35_55", "t35_55", "t55_100", "t55_100",
    ]


def test_synthetic_control_is_recoverable_by_its_own_oracle():
    """The injected world must actually be learnable in principle."""

    spec = trace.SyntheticSpec()
    query_bits, support_bits, query_desc, support_desc, protein, residual = make_inputs()
    data = {
        "query_bits": query_bits, "support_bits": support_bits,
        "query_desc": query_desc, "support_desc": support_desc,
        "protein": protein, "residual": residual,
        "base": torch.zeros(query_bits.shape[:2], device=DEVICE),
    }
    synthetic = trace.apply_synthetic(data, spec, seed=11)
    oracle = trace.synthetic_oracle_delta(synthetic, spec)
    label = synthetic["label"]
    isotropic = trace.analytic_delta(
        query_bits, support_bits, synthetic["residual"], estimator="krr", ridge=0.03
    )
    oracle_error = float(((oracle - label) ** 2).mean())
    isotropic_error = float(((isotropic - label) ** 2).mean())
    assert oracle_error < isotropic_error
