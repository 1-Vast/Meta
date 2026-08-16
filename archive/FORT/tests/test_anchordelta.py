from __future__ import annotations

import torch
from torch import nn

from model.anchordelta import (
    AnchorDelta,
    EncodedAnchorDelta,
    anchorabsolute,
    aggregateanchors,
    anchordeltaloss,
)


class _FakeInteraction(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.protein = nn.Linear(3, 4, bias=False)
        self.ligand = nn.Linear(5, 4, bias=False)

    def encodeprotein(self, tokens: torch.Tensor) -> torch.Tensor:
        return self.protein(tokens).mean(dim=0)

    def pairfromprotein(self, protein: torch.Tensor, ligand: torch.Tensor) -> torch.Tensor:
        return self.ligand(ligand) + protein


def test_antisymmetric_and_diagonal_zero() -> None:
    torch.manual_seed(7)
    model = AnchorDelta(feature_dim=4, hidden_dim=12)
    left = torch.randn(6, 4)
    right = torch.randn(6, 4)
    forward = model(left, right)
    reverse = model(right, left)
    assert torch.allclose(forward, -reverse, atol=1e-6, rtol=0.0)
    assert torch.equal(model(left, left), torch.zeros(6))


def test_matrix_and_anchor_permutation_invariance() -> None:
    torch.manual_seed(11)
    model = AnchorDelta(feature_dim=3, hidden_dim=8)
    query = torch.randn(4, 3)
    anchor = torch.randn(5, 3)
    labels = torch.linspace(4.0, 8.0, 5)
    deltas = model.matrix(query, anchor)
    expected = torch.stack([model(query, row.expand_as(query)) for row in anchor], dim=1)
    assert torch.allclose(deltas, expected)
    estimate, _ = anchorabsolute(labels, deltas)
    order = torch.tensor([2, 4, 0, 3, 1])
    permuted, _ = anchorabsolute(labels[order], deltas[:, order])
    assert torch.allclose(estimate, permuted)


def test_encoded_wrapper_freezes_interaction() -> None:
    torch.manual_seed(13)
    interaction = _FakeInteraction()
    model = EncodedAnchorDelta(interaction, feature_dim=4, hidden_dim=8)
    tokens = torch.randn(7, 3)
    query = torch.randn(4, 5)
    anchor = torch.randn(2, 5)
    output = model(tokens, query, anchor)
    assert output.shape == (4, 2)
    assert not any(parameter.requires_grad for parameter in interaction.parameters())
    assert any(parameter.requires_grad for parameter in model.head.parameters())


def test_pair_loss_masks_noise_floor() -> None:
    prediction = torch.tensor([1.0, -0.5, 0.0], requires_grad=True)
    target = torch.tensor([1.2, -0.2, 0.01])
    output = anchordeltaloss(prediction, target, order_weight=0.5, noise_floor=0.05)
    assert output["order"].ndim == 0
    output["loss"].backward()
    assert prediction.grad is not None


def test_weighted_aggregation_rejects_invalid_weights() -> None:
    values = torch.tensor([[1.0, 2.0], [3.0, 5.0]])
    estimate, variance = aggregateanchors(values)
    assert torch.equal(estimate, torch.tensor([1.5, 4.0]))
    assert torch.equal(variance, torch.tensor([0.25, 1.0]))
    try:
        aggregateanchors(values, method="weighted", weights=torch.zeros_like(values))
    except ValueError as error:
        assert "positive" in str(error)
    else:
        raise AssertionError("zero weights must be rejected")
