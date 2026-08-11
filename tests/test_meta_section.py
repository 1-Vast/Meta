import numpy as np
import pytest
import torch

from research.meta_fewshot.meta_section import IdentifiableMetaSection
from research.meta_fewshot.synthetic_meta_section import SyntheticConfig, run


def fixed_model(input_dim=6, section_dim=3, ridge=0.2):
    torch.manual_seed(4)
    model = IdentifiableMetaSection(input_dim, section_dim, ridge)
    with torch.no_grad():
        model.raw_basis.copy_(torch.eye(input_dim, section_dim, dtype=torch.float64))
    return model


def test_positive_ridge_and_dimensions_are_enforced():
    with pytest.raises(ValueError):
        IdentifiableMetaSection(6, 3, ridge=0)
    with pytest.raises(ValueError):
        IdentifiableMetaSection(6, 6)


def test_section_is_order_invariant_and_rank_bounded_by_support():
    model = fixed_model()
    phi = torch.tensor([[1., 0., 0., 0., 0., 0.],
                        [0., 1., 0., 0., 0., 0.]], dtype=torch.float64)
    residual = torch.tensor([2., -1.], dtype=torch.float64)
    query = torch.tensor([[1., 2., 0., 0., 0., 0.]], dtype=torch.float64)
    first = model.adapt(phi, residual)
    second = model.adapt(phi.flip(0), residual.flip(0))
    assert first.rank <= len(phi)
    assert torch.allclose(model.query(query, first)[0], model.query(query, second)[0])


def test_off_row_query_has_zero_correction_and_coverage():
    model = fixed_model()
    state = model.adapt(
        torch.tensor([[1., 0., 0., 0., 0., 0.]], dtype=torch.float64),
        torch.tensor([3.], dtype=torch.float64),
    )
    correction, coverage, _ = model.query(
        torch.tensor([[0., 1., 0., 0., 0., 0.]], dtype=torch.float64), state
    )
    assert correction.item() == pytest.approx(0.0, abs=1e-12)
    assert coverage.item() == pytest.approx(0.0, abs=1e-12)


def test_query_loss_reaches_shared_basis():
    model = fixed_model()
    support = torch.randn(3, 6, dtype=torch.float64)
    state = model.adapt(support, torch.randn(3, dtype=torch.float64))
    prediction, _, _ = model.query(torch.randn(4, 6, dtype=torch.float64), state)
    prediction.square().mean().backward()
    assert model.raw_basis.grad is not None
    assert torch.linalg.vector_norm(model.raw_basis.grad) > 0


def test_measurement_covariance_matches_monte_carlo():
    model = fixed_model(ridge=0.3)
    support = torch.tensor([[1., .2, 0., 0., 0., 0.],
                            [.1, 1., .3, 0., 0., 0.]], dtype=torch.float64)
    sigma = torch.diag(torch.tensor([0.04, 0.09], dtype=torch.float64))
    state = model.adapt(support, torch.zeros(2, dtype=torch.float64), sigma)
    query = torch.tensor([[.5, .8, .2, 0., 0., 0.]], dtype=torch.float64)
    _, _, analytic = model.query(query, state)
    rng = np.random.default_rng(9)
    noise = torch.from_numpy(rng.multivariate_normal(np.zeros(2), sigma.numpy(), 30000))
    coefficients = noise @ state.adaptation_map.T
    empirical = (coefficients @ model.coordinates(query).squeeze()).std(unbiased=True)
    assert empirical.item() == pytest.approx(analytic.item(), rel=0.025)


def test_bounded_support_noise_cannot_exceed_analytic_radius():
    model = fixed_model(ridge=0.3)
    support = torch.randn(3, 6, dtype=torch.float64)
    residual = torch.randn(3, dtype=torch.float64)
    query = torch.randn(5, 6, dtype=torch.float64)
    state = model.adapt(support, residual)
    bound = torch.tensor([0.1, 0.2, 0.15], dtype=torch.float64)
    radius = model.bounded_noise_radius(query, state, bound)
    for signs in torch.cartesian_prod(*[torch.tensor([-1., 1.])] * 3):
        changed = model.adapt(support, residual + signs * bound)
        baseline = model.query(query, state)[0]
        perturbed = model.query(query, changed)[0]
        assert torch.all((perturbed - baseline).abs() <= radius + 1e-12)


def test_small_synthetic_family_is_learnable(tmp_path):
    result = run(SyntheticConfig(
        source_tasks=32, evaluation_tasks=16, steps=220, tasks_per_step=4
    ), tmp_path / "result.json")
    assert result["checks"]["d_true_beats_d0"]
    assert result["checks"]["correct_beats_foreign"]
    assert result["query_labels_exposed_to_model"] is False
