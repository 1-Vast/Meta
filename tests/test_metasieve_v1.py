import pytest
import torch

from model.metasieve_v1 import MetaSieveV1, TaskScheduler, uniform_label_noise


def test_v1_episode_is_uncentered_positive_dual_ridge():
    torch.manual_seed(4)
    model = MetaSieveV1(input_dim=4, section_dim=2, ridge=0.7,
                        dtype=torch.float64)
    support_ligand = torch.randn(3, 4, dtype=torch.float64)
    support_pair = torch.randn(3, 4, dtype=torch.float64)
    query_ligand = torch.randn(5, 4, dtype=torch.float64)
    query_pair = torch.randn(5, 4, dtype=torch.float64)
    support_y = torch.randn(3, dtype=torch.float64)
    predicted = model.episode(
        support_ligand, support_pair, support_y, query_ligand, query_pair)
    support_population, support_coordinates = model.components(
        support_ligand, support_pair)
    query_population, query_coordinates = model.components(query_ligand, query_pair)
    expected = query_population + query_coordinates @ support_coordinates.T @ \
        torch.linalg.solve(
            support_coordinates @ support_coordinates.T
            + 0.7 * torch.eye(3, dtype=torch.float64),
            support_y - support_population,
        )
    assert torch.allclose(predicted, expected, atol=1e-12)


def test_v1_meta_gradient_crosses_analytic_support_solve():
    torch.manual_seed(8)
    model = MetaSieveV1(input_dim=3, section_dim=2, ridge=1.0,
                        dtype=torch.float64)
    loss = model.episode(
        torch.randn(2, 3, dtype=torch.float64),
        torch.randn(2, 3, dtype=torch.float64),
        torch.randn(2, dtype=torch.float64),
        torch.randn(4, 3, dtype=torch.float64),
        torch.randn(4, 3, dtype=torch.float64),
    ).square().mean()
    loss.backward()
    assert model.population.weight.grad is not None
    assert model.raw_basis.grad is not None
    assert model.population_coordinate.grad is not None


def test_batched_episode_matches_individual_episodes_and_gradients():
    torch.manual_seed(21)
    model = MetaSieveV1(input_dim=5, section_dim=2, ridge=0.8,
                        dtype=torch.float64)
    support_ligand = torch.randn(4, 3, 5, dtype=torch.float64)
    support_pair = torch.randn(4, 3, 5, dtype=torch.float64)
    support_y = torch.randn(4, 3, dtype=torch.float64)
    query_ligand = torch.randn(4, 6, 5, dtype=torch.float64)
    query_pair = torch.randn(4, 6, 5, dtype=torch.float64)
    batched, _ = model.batched_episode(
        support_ligand, support_pair, support_y, query_ligand, query_pair)
    individual = torch.stack([
        model.episode(
            support_ligand[index], support_pair[index], support_y[index],
            query_ligand[index], query_pair[index])
        for index in range(4)
    ])
    assert torch.allclose(batched, individual, atol=1e-12)
    batched.square().mean().backward()
    batched_gradients = [parameter.grad.detach().clone()
                         for parameter in model.parameters()]
    model.zero_grad(set_to_none=True)
    individual.square().mean().backward()
    assert all(torch.allclose(parameter.grad, expected, atol=1e-11)
               for parameter, expected in zip(model.parameters(), batched_gradients))


def test_support_only_section_has_no_pair_dependent_zero_shot_path():
    torch.manual_seed(31)
    model = MetaSieveV1(
        input_dim=4, section_dim=2, ridge=1.0,
        support_only_section=True, dtype=torch.float64)
    ligand = torch.randn(5, 4, dtype=torch.float64)
    pair_left = torch.randn(5, 4, dtype=torch.float64)
    pair_right = torch.randn(5, 4, dtype=torch.float64)
    left, _ = model.components(ligand, pair_left)
    right, _ = model.components(ligand, pair_right)
    assert torch.equal(left, right)
    assert not hasattr(model, "population_coordinate")

    support_ligand = ligand[:2]
    support_pair = pair_left[:2]
    query_ligand = ligand[2:]
    query_pair = pair_left[2:]
    first = model.episode(
        support_ligand, support_pair, torch.tensor([0.0, 1.0]),
        query_ligand, query_pair)
    second = model.episode(
        support_ligand, support_pair, torch.tensor([1.0, 0.0]),
        query_ligand, query_pair)
    assert not torch.equal(first, second)


def test_narrow_population_mlp_preserves_episode_contract():
    torch.manual_seed(32)
    model = MetaSieveV1(
        input_dim=4, section_dim=2, ridge=0.5,
        support_only_section=True, population_hidden_dim=3,
        dtype=torch.float64)
    support_ligand = torch.randn(2, 4, dtype=torch.float64)
    support_pair = torch.randn(2, 4, dtype=torch.float64)
    query_ligand = torch.randn(3, 4, dtype=torch.float64)
    query_pair = torch.randn(3, 4, dtype=torch.float64)
    prediction = model.episode(
        support_ligand, support_pair, torch.randn(2, dtype=torch.float64),
        query_ligand, query_pair)
    assert prediction.shape == (3,)
    prediction.square().mean().backward()
    assert all(parameter.grad is not None for parameter in model.parameters())


def test_narrow_pair_map_remains_inside_orthonormal_meta_section():
    torch.manual_seed(33)
    model = MetaSieveV1(
        input_dim=6, section_dim=2, ridge=1.0,
        support_only_section=True, population_hidden_dim=4,
        pair_hidden_dim=3, dtype=torch.float64)
    basis = model.basis()
    assert basis.shape == (3, 2)
    assert torch.allclose(
        basis.T @ basis, torch.eye(2, dtype=torch.float64), atol=1e-12)
    loss = model.episode(
        torch.randn(2, 6, dtype=torch.float64),
        torch.randn(2, 6, dtype=torch.float64),
        torch.randn(2, dtype=torch.float64),
        torch.randn(3, 6, dtype=torch.float64),
        torch.randn(3, 6, dtype=torch.float64),
    ).square().mean()
    loss.backward()
    assert model.raw_basis.grad is not None
    assert all(parameter.grad is not None for parameter in model.pair_encoder.parameters())


def test_scheduler_samples_without_replacement_and_is_normalized():
    scheduler = TaskScheduler(hidden_dim=5, dtype=torch.float64)
    statistics = torch.tensor([
        [0.2, 0.8, 0.0], [1.0, -0.2, 0.2], [0.5, 0.3, 0.4],
        [1.5, -0.8, 0.6], [0.1, 0.9, 0.8], [0.7, 0.0, 1.0],
    ], dtype=torch.float64)
    probability = scheduler.probabilities(statistics)
    assert probability.sum().item() == pytest.approx(1.0)
    selected = scheduler.sample(
        statistics, 4, torch.Generator().manual_seed(3))
    assert len(torch.unique(selected)) == 4


def test_uniform_label_noise_is_deterministic_zero_mean_contract():
    labels = torch.zeros(20000)
    left = uniform_label_noise(
        labels, 0.2, generator=torch.Generator().manual_seed(11))
    right = uniform_label_noise(
        labels, 0.2, generator=torch.Generator().manual_seed(11))
    assert torch.equal(left, right)
    assert abs(float(left.mean())) < 0.005
    assert float(left.std()) == pytest.approx(0.2, abs=0.005)
    assert uniform_label_noise(
        labels, 0.0, generator=torch.Generator().manual_seed(1)) is labels
