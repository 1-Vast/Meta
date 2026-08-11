import numpy as np
import torch

from research.meta_fewshot.train_main_v0 import (
    MetaSectionRegressor,
    bootstrap_contrast,
    cluster_bootstrap_contrast,
    draw_episode,
)


def test_episode_does_not_accept_query_labels_and_backpropagates_through_solve():
    torch.manual_seed(1)
    model = MetaSectionRegressor(8, 3, 0.1)
    ligand = torch.randn(9, 8)
    family = torch.randn(9, 8)
    y = torch.randn(9)
    prediction = model.episode(ligand[:5], family[:5], y[:5], ligand[5:], family[5:])
    prediction.square().mean().backward()
    assert model.raw_basis.grad is not None
    assert torch.isfinite(model.raw_basis.grad).all()


def test_d0_prediction_is_support_independent():
    model = MetaSectionRegressor(4, 0, 0.1)
    ligand = torch.randn(8, 4)
    family = torch.randn(8, 4)
    first = model.episode(ligand[:5], family[:5], torch.zeros(5), ligand[5:], family[5:])
    second = model.episode(ligand[:5], family[:5], torch.ones(5), ligand[5:], family[5:])
    assert torch.equal(first, second)


def test_production_candidate_uses_uncentered_positive_dual_ridge():
    model = MetaSectionRegressor(4, 2, 0.7)
    with torch.no_grad():
        model.population.weight.zero_()
        model.population.bias.fill_(0.25)
        model.raw_basis.copy_(torch.eye(4, 2))
        model.population_coordinate.copy_(torch.tensor([0.2, -0.1]))
    support_ligand = torch.zeros(3, 4)
    query_ligand = torch.zeros(2, 4)
    support_family = torch.tensor([
        [2., 0., 0., 0.],
        [1., 1., 0., 0.],
        [3., 1., 0., 0.],
    ])
    query_family = torch.tensor([
        [1.5, 0.5, 0., 0.],
        [4.0, 1.0, 0., 0.],
    ])
    support_y = torch.tensor([1.4, -0.2, 2.1])

    support_population, matrix = model.components(support_ligand, support_family)
    query_population, query_matrix = model.components(query_ligand, query_family)
    residual = support_y - support_population
    identity = torch.eye(len(support_y))
    dual = torch.linalg.solve(matrix @ matrix.T + 0.7 * identity, residual)
    expected = query_population + query_matrix @ matrix.T @ dual

    actual = model.episode(
        support_ligand, support_family, support_y, query_ligand, query_family)
    assert torch.allclose(actual, expected)
    assert not torch.allclose(matrix.mean(0), torch.zeros(2))


def test_draw_episode_has_disjoint_support_and_query():
    support, query = draw_episode(np.arange(12), np.random.default_rng(1), 5)
    assert len(support) == 5
    assert set(support).isdisjoint(query)


def test_bootstrap_contrast_uses_target_level_paired_loss():
    losses = {
        ("correct", "a"): 1.0, ("control", "a"): 2.0,
        ("correct", "b"): 1.0, ("control", "b"): 3.0,
    }
    result = bootstrap_contrast(losses, "correct", "control", 999, 1)
    assert result["mean_mse_reduction"] == 1.5
    assert result["pass"]


def test_cluster_bootstrap_averages_targets_within_cluster_first():
    losses = {
        ("correct", "a"): 1.0, ("control", "a"): 2.0,
        ("correct", "b"): 1.0, ("control", "b"): 2.0,
        ("correct", "c"): 1.0, ("control", "c"): 5.0,
    }
    result = cluster_bootstrap_contrast(
        losses, {"a": "x", "b": "x", "c": "y"},
        "correct", "control", 999, 1)
    assert result["clusters"] == 2
    assert result["cluster_macro_mse_reduction"] == 2.5
