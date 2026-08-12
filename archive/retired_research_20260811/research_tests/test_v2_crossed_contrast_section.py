"""Algebraic guarantees for the V2 crossed-contrast candidate.

These tests verify arithmetic and identifiability claims only.  They are not a
biological PASS and they authorize no training, no Gate and no production
migration.
"""
from __future__ import annotations

import math

import pytest

torch = pytest.importorskip("torch")

from research.meta_fewshot.v2_crossed_contrast_section import (  # noqa: E402
    AnisotropicPairCoordinate,
    CrossedContrastSection,
    Rectangle,
    calibration_leakage,
    centered_section,
    crossed_contrast_loss,
    measured_crossed_difference,
    predicted_crossed_difference,
    wrong_wrong_invariance,
)

DTYPE = torch.float64


def _model(ligand_dim=6, protein_dim=4, section_dim=2, ridge=1.0, seed=0):
    torch.manual_seed(seed)
    model = CrossedContrastSection(ligand_dim, protein_dim, section_dim, ridge)
    return model.to(DTYPE)


# --------------------------------------------------------------------------- #
# calibration separation
# --------------------------------------------------------------------------- #
def test_constant_support_shift_moves_prediction_by_exactly_that_constant():
    model = _model()
    torch.manual_seed(1)
    support_ligand = torch.randn(5, 6, dtype=DTYPE)
    query_ligand = torch.randn(7, 6, dtype=DTYPE)
    support_protein = torch.randn(1, 4, dtype=DTYPE).expand(5, -1)
    query_protein = support_protein[:1].expand(7, -1)
    support_y = torch.randn(5, dtype=DTYPE)
    delta = calibration_leakage(model, support_protein, support_ligand, support_y,
                                query_protein, query_ligand, shift=2.5)
    assert torch.allclose(delta, torch.full_like(delta, 2.5), atol=1e-9)


def test_section_label_degrees_of_freedom_never_exceed_k():
    torch.manual_seed(2)
    for k in (1, 2, 3, 5):
        coordinates = torch.randn(k, 5, dtype=DTYPE)
        output = centered_section(coordinates, torch.randn(k, dtype=DTYPE),
                                  torch.randn(3, 5, dtype=DTYPE), ridge=0.5)
        assert output.rank <= max(k - 1, 0)
        assert 1 + output.rank <= k


def test_centered_section_rejects_non_positive_ridge():
    coordinates = torch.randn(3, 2, dtype=DTYPE)
    with pytest.raises(ValueError):
        centered_section(coordinates, torch.zeros(3, dtype=DTYPE),
                         coordinates, ridge=0.0)


# --------------------------------------------------------------------------- #
# the centring-gauge collapse: why the minimal repair is not sufficient
# --------------------------------------------------------------------------- #
def _affine_coordinates(ligand: torch.Tensor, transform: torch.Tensor,
                        offset: torch.Tensor) -> torch.Tensor:
    return ligand @ transform.T + offset


def test_centered_section_is_blind_to_protein_when_the_metric_is_shared():
    """m(P,L) = A_P phi(L) + b_P with A_P orthogonal gives identical predictions.

    This is the executable form of the centring-gauge collapse: after centring,
    the protein survives only through G_P = A_P^T A_P.  Orthogonal A_P share one
    metric, so correct/correct and wrong/wrong coincide exactly.
    """
    torch.manual_seed(3)
    ligand_support = torch.randn(5, 4, dtype=DTYPE)
    ligand_query = torch.randn(6, 4, dtype=DTYPE)
    residual = torch.randn(5, dtype=DTYPE)
    first = torch.linalg.qr(torch.randn(4, 4, dtype=DTYPE)).Q[:, :3].T
    second = torch.linalg.qr(torch.randn(4, 4, dtype=DTYPE)).Q[:, :3].T
    offset_a = torch.randn(3, dtype=DTYPE)
    offset_b = torch.randn(3, dtype=DTYPE)

    left = centered_section(
        _affine_coordinates(ligand_support, first, offset_a), residual,
        _affine_coordinates(ligand_query, first, offset_a), ridge=0.7).prediction
    right = centered_section(
        _affine_coordinates(ligand_support, second, offset_b), residual,
        _affine_coordinates(ligand_query, second, offset_b), ridge=0.7).prediction
    assert torch.allclose(left, right, atol=1e-8)


def test_anisotropic_metric_restores_protein_dependence():
    """A protein-specific, non-orthogonal metric does change the prediction."""
    torch.manual_seed(4)
    ligand_support = torch.randn(5, 4, dtype=DTYPE)
    ligand_query = torch.randn(6, 4, dtype=DTYPE)
    residual = torch.randn(5, dtype=DTYPE)
    base = torch.linalg.qr(torch.randn(4, 4, dtype=DTYPE)).Q[:, :3].T
    stretched = torch.diag(torch.tensor([4.0, 1.0, 0.25], dtype=DTYPE)) @ base
    zero = torch.zeros(3, dtype=DTYPE)

    left = centered_section(_affine_coordinates(ligand_support, base, zero), residual,
                            _affine_coordinates(ligand_query, base, zero),
                            ridge=0.7).prediction
    right = centered_section(_affine_coordinates(ligand_support, stretched, zero),
                             residual,
                             _affine_coordinates(ligand_query, stretched, zero),
                             ridge=0.7).prediction
    assert not torch.allclose(left, right, atol=1e-6)


def test_wrong_wrong_invariance_is_a_reported_falsifier():
    model = _model(seed=5)
    torch.manual_seed(5)
    support_ligand = torch.randn(5, 6, dtype=DTYPE)
    query_ligand = torch.randn(4, 6, dtype=DTYPE)
    support_y = torch.randn(5, dtype=DTYPE)
    correct = torch.randn(1, 4, dtype=DTYPE)
    wrong = torch.randn(1, 4, dtype=DTYPE)
    gap = wrong_wrong_invariance(model, support_ligand, support_y, query_ligand,
                                 correct, wrong)
    assert gap.shape == (4,)
    assert torch.isfinite(gap).all()


# --------------------------------------------------------------------------- #
# crossed contrast
# --------------------------------------------------------------------------- #
def test_measured_crossed_difference_annihilates_additive_effects():
    alpha = torch.tensor([1.5, -0.25], dtype=DTYPE)
    beta = torch.tensor([0.75, 2.0], dtype=DTYPE)
    values = alpha[:, None] + beta[None, :]
    assert torch.allclose(measured_crossed_difference(values),
                          torch.zeros((), dtype=DTYPE), atol=1e-12)


def _rectangle(seed: int, ligand_dim=6, protein_dim=4, k=5) -> Rectangle:
    torch.manual_seed(seed)
    return Rectangle(
        protein_a=torch.randn(protein_dim, dtype=DTYPE),
        protein_b=torch.randn(protein_dim, dtype=DTYPE),
        ligand_1=torch.randn(ligand_dim, dtype=DTYPE),
        ligand_2=torch.randn(ligand_dim, dtype=DTYPE),
        values=torch.randn(2, 2, dtype=DTYPE),
        support_a_ligand=torch.randn(k, ligand_dim, dtype=DTYPE),
        support_a_y=torch.randn(k, dtype=DTYPE),
        support_b_ligand=torch.randn(k, ligand_dim, dtype=DTYPE),
        support_b_y=torch.randn(k, dtype=DTYPE),
    )


def test_ligand_prior_and_task_calibration_cannot_produce_a_crossed_prediction():
    """With the section switched off, the predicted second difference is zero.

    Large ridge drives the centred correction to zero, leaving exactly
    ``mu_L + b_t``.  The predicted crossed difference is then identically zero
    whatever the ligand prior is, so the X-CON objective has no gradient through
    either the ligand prior or the calibration channel.
    """
    model = _model(ridge=1e12, seed=6)
    rectangle = _rectangle(6)
    predicted = predicted_crossed_difference(model, rectangle)
    assert abs(float(predicted)) < 1e-6


def test_crossed_prediction_is_invariant_to_a_constant_population_shift():
    model = _model(seed=7)
    rectangle = _rectangle(7)
    before = predicted_crossed_difference(model, rectangle)
    with torch.no_grad():
        model.ligand_population.bias += 3.0
    after = predicted_crossed_difference(model, rectangle)
    assert torch.allclose(before, after, atol=1e-9)


def test_crossed_prediction_is_invariant_to_a_per_protein_label_offset():
    model = _model(seed=8)
    rectangle = _rectangle(8)
    before = predicted_crossed_difference(model, rectangle)
    shifted = Rectangle(
        **{**rectangle.__dict__, "support_a_y": rectangle.support_a_y + 1.75})
    after = predicted_crossed_difference(model, shifted)
    assert torch.allclose(before, after, atol=1e-9)


def test_crossed_contrast_loss_reaches_the_anisotropic_family():
    model = _model(seed=9)
    loss = crossed_contrast_loss(model, [_rectangle(9), _rectangle(10)])
    loss.backward()
    assert model.coordinate.family.grad is not None
    assert float(model.coordinate.family.grad.abs().sum()) > 0.0
    assert float(model.coordinate.selector.weight.grad.abs().sum()) > 0.0


def test_crossed_contrast_loss_requires_rectangles():
    with pytest.raises(ValueError):
        crossed_contrast_loss(_model(seed=11), [])


# --------------------------------------------------------------------------- #
# gauge invariance of the reported metric
# --------------------------------------------------------------------------- #
def test_reported_metric_is_invariant_to_the_latent_gauge():
    torch.manual_seed(12)
    coordinate = AnisotropicPairCoordinate(6, 3, 4, family_size=5).to(DTYPE)
    protein = torch.randn(2, 4, dtype=DTYPE)
    before = coordinate.metric(protein)
    rotation = torch.linalg.qr(torch.randn(3, 3, dtype=DTYPE)).Q
    with torch.no_grad():
        coordinate.base.copy_(coordinate.base @ rotation)
        coordinate.family.copy_(coordinate.family @ rotation)
    after = coordinate.metric(protein)
    assert torch.allclose(before, after, atol=1e-9)


def test_family_selector_is_bounded():
    torch.manual_seed(13)
    coordinate = AnisotropicPairCoordinate(6, 2, 4, family_size=8).to(DTYPE)
    gate = coordinate.gate(torch.randn(32, 4, dtype=DTYPE) * 50.0)
    assert float(gate.norm(dim=-1).max()) <= 1.0 + 1e-9


def test_section_dimension_respects_the_support_budget():
    with pytest.raises(ValueError):
        AnisotropicPairCoordinate(6, 6, 4)
    with pytest.raises(ValueError):
        AnisotropicPairCoordinate(6, 0, 4)


def test_metric_matches_an_explicit_product():
    torch.manual_seed(14)
    coordinate = AnisotropicPairCoordinate(5, 2, 3, family_size=4).to(DTYPE)
    protein = torch.randn(1, 3, dtype=DTYPE)
    maps = coordinate.maps(protein)[0]
    assert torch.allclose(coordinate.metric(protein)[0], maps @ maps.T, atol=1e-10)
    assert not math.isnan(float(maps.sum()))
