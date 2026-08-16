from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import torch

from research.a2s.a2s_interaction_state_gate import (
    InteractionChannels,
    active_dimensions,
    adapted_prediction,
    assert_source_roles,
    fixed_rotation,
    norm_matched_transplant,
    pairwise_proper_loss,
    solve_state,
)


def test_budget_dimensions_are_explicit() -> None:
    assert active_dimensions(1) == 0
    assert active_dimensions(3) == 1
    assert active_dimensions(5) == 2
    with pytest.raises(ValueError):
        active_dimensions(2)


def test_k1_and_residual_null_are_exact_noops() -> None:
    phi_s = torch.tensor([[0.2, -0.4]])
    phi_q = torch.tensor([[0.9, 0.3], [-0.2, 0.5]])
    base = torch.tensor([1.0, 2.0])
    prediction, state, _ = adapted_prediction(phi_s, phi_q, torch.tensor([3.2]), base, 1)
    assert torch.equal(prediction, base)
    assert torch.equal(state, torch.zeros_like(state))

    phi_s = torch.tensor([[0.1, 0.3], [0.4, -0.2], [-0.5, 0.8]])
    prediction, state, _ = adapted_prediction(
        phi_s, phi_q, torch.ones(3), base, 3
    )
    assert torch.equal(prediction, base)
    assert torch.equal(state, torch.zeros_like(state))


def test_noiseless_state_recovers_query_direction() -> None:
    phi_s = torch.tensor([[-1.0, 0.0], [0.0, 0.0], [1.0, 0.0]])
    residual = torch.tensor([-1.0, 0.0, 1.0])
    phi_q = torch.tensor([[-0.8, 0.0], [0.8, 0.0]])
    prediction, state, _ = adapted_prediction(
        phi_s, phi_q, residual, torch.zeros(2), 3, ridge=1e-6
    )
    assert state[0] > 0.99
    assert prediction[1] > prediction[0]


def test_dense_eb_prediction_is_rotation_invariant() -> None:
    generator = torch.Generator().manual_seed(17)
    phi_s = torch.randn(5, 2, generator=generator)
    phi_q = torch.randn(9, 2, generator=generator)
    residual = torch.randn(5, generator=generator)
    base = torch.randn(9, generator=generator)
    rotation = fixed_rotation(91).cpu()
    original, _, _ = adapted_prediction(phi_s, phi_q, residual, base, 5)
    rotated, _, _ = adapted_prediction(
        phi_s @ rotation, phi_q @ rotation, residual, base, 5
    )
    assert torch.allclose(original, rotated, atol=1e-6, rtol=1e-6)


def test_norm_matched_transplant_preserves_centered_norm() -> None:
    correct = torch.tensor([1.0, 3.0, -2.0, 0.5, 2.5])
    wrong = torch.tensor([-4.0, 0.0, 1.0, 2.0, 8.0])
    transplanted = norm_matched_transplant(correct, wrong)
    assert float(transplanted.mean()) == pytest.approx(0.0, abs=1e-6)
    assert float(torch.linalg.vector_norm(transplanted)) == pytest.approx(
        float(torch.linalg.vector_norm(correct - correct.mean())), rel=1e-6
    )


def test_pairwise_proper_loss_rewards_correct_order() -> None:
    label = torch.tensor([3.0, 2.0, 1.0])
    correct = pairwise_proper_loss(label, torch.tensor([2.0, 1.0, 0.0]))
    reversed_loss = pairwise_proper_loss(label, torch.tensor([0.0, 1.0, 2.0]))
    assert correct < reversed_loss


def test_segment_encoder_is_permutation_invariant() -> None:
    torch.manual_seed(9)
    model = InteractionChannels(12, "segment")
    ligand = torch.randn(7, 12)
    segments = torch.randn(32, 32)
    first = model(ligand, segments)
    second = model(ligand, segments[torch.randperm(32)])
    assert first.shape == (7, 2)
    assert torch.allclose(first, second, atol=1e-6, rtol=1e-6)


def test_analytic_state_solve_backpropagates_into_interaction_channels() -> None:
    torch.manual_seed(77)
    model = InteractionChannels(12, "segment")
    ligand = torch.randn(12, 12)
    segments = torch.randn(32, 32)
    phi = model(ligand, segments)
    prediction, state, _ = adapted_prediction(
        phi[:5], phi[5:], torch.randn(5), torch.randn(7), 5
    )
    loss = pairwise_proper_loss(torch.randn(7), prediction)
    loss.backward()

    assert state.grad_fn is not None
    assert model.output.weight.grad is not None
    assert model.segment.weight.grad is not None
    assert float(model.output.weight.grad.norm()) > 0.0
    assert float(model.segment.weight.grad.norm()) > 0.0


def test_source_role_guard_rejects_locked_role() -> None:
    assert_source_roles(pd.DataFrame({"role": ["fit", "probe"]}))
    with pytest.raises(AssertionError):
        assert_source_roles(pd.DataFrame({"role": ["fit", "probe", "locked"]}))
