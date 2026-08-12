import numpy as np
import pytest
import torch

from research.meta_fewshot.r2_calibration_orthogonal_section import (
    calibration_orthogonal_prediction,
    ridge_support_weights,
)
from research.meta_fewshot.r2_e2_tbasis_decomposition import (
    additive_projection,
    bipartite_two_core,
)
from research.meta_fewshot.r2_e3_crossed_census import interaction_df
from research.meta_fewshot.r2_reserved_fiber_section import (
    ReservedFiberSection,
    reserved_exposure_from_frames,
)


def unregularized_kernel(support, query):
    return query @ torch.linalg.solve(support.T @ support, support.T)


def test_ridge_invariance_boundary_is_gl_only_without_regularization():
    torch.manual_seed(3)
    support = torch.randn(5, 2, dtype=torch.float64)
    query = torch.randn(4, 2, dtype=torch.float64)
    transform = torch.tensor([[2.0, 0.7], [0.0, 0.5]], dtype=torch.float64)
    assert torch.allclose(
        unregularized_kernel(support, query),
        unregularized_kernel(support @ transform, query @ transform),
        atol=1e-10,
    )
    regularized = ridge_support_weights(support, query, 1.0)
    transformed = ridge_support_weights(support @ transform, query @ transform, 1.0)
    assert not torch.allclose(regularized, transformed, atol=1e-5)


def test_positive_ridge_preserves_orthogonal_coordinate_change():
    torch.manual_seed(5)
    support = torch.randn(5, 2, dtype=torch.float64)
    query = torch.randn(4, 2, dtype=torch.float64)
    q, _ = torch.linalg.qr(torch.randn(2, 2, dtype=torch.float64))
    assert torch.allclose(
        ridge_support_weights(support, query, 0.7),
        ridge_support_weights(support @ q, query @ q, 0.7),
        atol=1e-10,
    )


def test_calibration_orthogonal_section_separates_constant_residual():
    torch.manual_seed(7)
    support = torch.randn(5, 2, dtype=torch.float64)
    query = torch.randn(3, 2, dtype=torch.float64)
    mu_s = torch.randn(5, dtype=torch.float64)
    mu_q = torch.randn(3, dtype=torch.float64)
    prediction, calibration, specific = calibration_orthogonal_prediction(
        mu_s, support, mu_s + 1.25, mu_q, query, 1.0)
    assert calibration.item() == pytest.approx(1.25)
    assert torch.allclose(specific, torch.zeros_like(specific), atol=1e-12)
    assert torch.allclose(prediction, mu_q + 1.25, atol=1e-12)


def test_calibration_orthogonal_section_is_coordinate_translation_invariant():
    torch.manual_seed(11)
    support = torch.randn(5, 2, dtype=torch.float64)
    query = torch.randn(3, 2, dtype=torch.float64)
    mu_s, y_s = torch.randn(5, dtype=torch.float64), torch.randn(5, dtype=torch.float64)
    mu_q = torch.randn(3, dtype=torch.float64)
    base = calibration_orthogonal_prediction(mu_s, support, y_s, mu_q, query, 0.5)[0]
    shift = torch.tensor([12.0, -3.0], dtype=torch.float64)
    moved = calibration_orthogonal_prediction(
        mu_s, support + shift, y_s, mu_q, query + shift, 0.5)[0]
    assert torch.allclose(base, moved, atol=1e-10)


def test_additive_projection_and_two_core_on_synthetic_design():
    protein = np.repeat(np.arange(3), 4)
    ligand = np.tile(np.arange(4), 3)
    phi = np.stack([2.0 * protein + ligand, -protein + 0.5 * ligand], axis=1)
    result = additive_projection(phi.astype(float), protein, ligand)
    assert result["interaction_residual_fraction"] < 1e-20
    assert bipartite_two_core(protein, ligand).all()


def test_interaction_df_is_computed_within_assay_blocks():
    protein = np.asarray(["p0", "p0", "p1", "p1"])
    ligand = np.asarray(["l0", "l1", "l0", "l1"])
    block = np.asarray(["assay", "assay", "assay", "assay"])
    assert interaction_df(protein, ligand, block) == (1, 1)


def test_nonconstant_partner_map_does_not_guarantee_reserved_exposure():
    adaptable_s = torch.tensor([[1.0], [2.0]], dtype=torch.float64)
    reserved_s = torch.tensor([[2.0], [4.0]], dtype=torch.float64)
    adaptable_q = torch.tensor([[3.0]], dtype=torch.float64)
    gram = adaptable_s.T @ adaptable_s
    transfer = torch.linalg.solve(
        gram + torch.eye(1, dtype=torch.float64), adaptable_s.T @ reserved_s)
    reserved_q = adaptable_q @ transfer
    exposure = reserved_exposure_from_frames(
        adaptable_s, reserved_s, adaptable_q, reserved_q, ridge=1.0)
    assert torch.allclose(exposure, torch.zeros_like(exposure), atol=1e-12)
    partner_difference = torch.tensor([4.0], dtype=torch.float64)
    assert torch.allclose(exposure @ partner_difference, torch.zeros(1, dtype=torch.float64))


def test_rfms_accepts_singleton_batched_protein_and_rejects_too_few_supports():
    torch.manual_seed(13)
    model = ReservedFiberSection(4, 3, d=3, d_support=2, ridge=1.0).double()
    support_ligand = torch.randn(2, 4, dtype=torch.float64)
    query_ligand = torch.randn(3, 4, dtype=torch.float64)
    support_y = torch.randn(2, dtype=torch.float64)
    prediction, residual = model.episode(
        support_ligand, support_y, torch.randn(1, 3, dtype=torch.float64), query_ligand)
    assert prediction.shape == (3,)
    assert residual.ndim == 0
    with pytest.raises(ValueError, match="support size"):
        model.episode(support_ligand[:1], support_y[:1], torch.randn(3), query_ligand)
