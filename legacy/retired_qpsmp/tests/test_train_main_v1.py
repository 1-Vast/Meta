import numpy as np
import pytest
import torch

from model.metasieve_v1 import MetaSieveV1, uniform_label_noise
from research.meta_fewshot.train_main_v1 import (
    _flat_batched_gradients,
    _difference_mse,
    add_episodewise_label_noise,
    build_tasks,
    cluster_tasks,
    different_cluster_donors,
    draw_episode,
    draw_nested_episode,
    gradient_cosine,
    gradient_cosine_rows,
    pack_episode_indices,
    pretrain_population,
    residualize_pair_against_ligand,
    sample_cluster_balanced,
    V1TrainConfig,
)


def test_population_pretraining_is_source_only_and_reduces_source_loss():
    torch.manual_seed(7)
    model = MetaSieveV1(
        input_dim=4, section_dim=2, support_only_section=True,
        population_hidden_dim=8)
    cells = [
        {"split": "meta_train" if index < 32 else "meta_val"}
        for index in range(40)]
    ligand = torch.randn(40, 4)
    target = ligand[:, 0] - 0.5 * ligand[:, 1]
    tensors = {"ligand": ligand, "y": target}
    before = torch.nn.functional.mse_loss(
        model.population(ligand[:32]).squeeze(-1), target[:32])
    detail = pretrain_population(
        model, cells, tensors, seed=11,
        config=V1TrainConfig(
            population_pretrain_steps=100,
            population_pretrain_batch_size=16,
            population_pretrain_learning_rate=0.01))
    after = torch.nn.functional.mse_loss(
        model.population(ligand[:32]).squeeze(-1), target[:32])
    assert detail["steps"] == 100
    assert after < before


def test_pair_residualization_fits_source_ligand_projection_only():
    ligand = torch.tensor([
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [2.0, -1.0],
    ])
    weights = torch.tensor([[2.0, -1.0], [0.5, 3.0]])
    pair = ligand @ weights
    residual = residualize_pair_against_ligand(
        pair, ligand, torch.tensor([0, 1, 2]), ridge=1e-6)
    assert torch.allclose(residual[:3], torch.zeros(3, 2), atol=5e-6)
    assert torch.isfinite(residual).all()


def test_v1_task_builder_supports_all_cold_target_shots():
    cells = []
    for target, count, cluster in (("a", 4, "x"), ("b", 6, "y"), ("c", 8, "z")):
        for index in range(count):
            cells.append({"split": "meta_test", "target_id": target,
                          "protein_group_40": cluster, "cell_id": f"{target}-{index}"})
    assert set(build_tasks(cells, "meta_test", 1, 3)) == {"a", "b", "c"}
    assert set(build_tasks(cells, "meta_test", 3, 3)) == {"b", "c"}
    assert set(build_tasks(cells, "meta_test", 5, 3)) == {"c"}
    tasks = build_tasks(cells, "meta_test", 1, 3)
    assert cluster_tasks(cells, tasks) == {"x": ["a"], "y": ["b"], "z": ["c"]}


def test_cluster_sampling_and_episode_draw_are_deterministic_and_disjoint():
    clusters = {f"c{i}": [f"t{i}"] for i in range(10)}
    left = sample_cluster_balanced(clusters, np.random.default_rng(5), 6)
    right = sample_cluster_balanced(clusters, np.random.default_rng(5), 6)
    assert left == right
    assert len(set(left)) == 6
    support, query = draw_episode(np.arange(12), np.random.default_rng(8), 3, 5)
    assert len(support) == 3 and len(query) == 5
    assert not set(support) & set(query)


def test_cold_target_k_sweep_uses_nested_support_and_shared_query():
    indices = np.arange(20)
    episodes = {
        k: draw_nested_episode(indices, np.random.default_rng(17), k, 5, 6)
        for k in (1, 2, 3, 5)}
    support5, query5 = episodes[5]
    for k, (support, query) in episodes.items():
        assert np.array_equal(support, support5[:k])
        assert np.array_equal(query, query5)


def test_foreign_support_donor_is_from_a_different_protein_cluster():
    cells = [
        {"target_id": "a", "protein_group_40": "g1"},
        {"target_id": "b", "protein_group_40": "g1"},
        {"target_id": "c", "protein_group_40": "g2"},
    ]
    donors = different_cluster_donors(cells, ["a", "b", "c"])
    cluster = {row["target_id"]: row["protein_group_40"] for row in cells}
    assert all(cluster[target] != cluster[donor] for target, donor in donors.items())


def test_gradient_cosine_and_difference_metric_contracts():
    assert gradient_cosine(torch.tensor([1., 0.]), torch.tensor([2., 0.])) == 1
    assert gradient_cosine(torch.zeros(2), torch.ones(2)) == 0
    y = np.array([1.0, 2.0, 4.0])
    assert _difference_mse(y, y + 10.0) == pytest.approx(0.0)
    assert _difference_mse(y[:1], y[:1]) != _difference_mse(y[:1], y[:1])


def test_batched_gradient_and_episode_packing_match_scalar_contracts():
    parameter = torch.nn.Parameter(torch.tensor([0.5, -1.0], dtype=torch.float64))
    weights = torch.tensor([[1.0, 2.0], [-3.0, 0.5], [0.2, -0.7]],
                           dtype=torch.float64)
    losses = ((weights * parameter).sum(dim=1)).square()
    batched = _flat_batched_gradients(
        losses, (parameter,), retain_graph=True)
    scalar = torch.stack([
        torch.autograd.grad(loss, parameter, retain_graph=True)[0]
        for loss in losses
    ])
    assert torch.allclose(batched, scalar, atol=1e-12)
    assert torch.allclose(
        gradient_cosine_rows(batched, scalar), torch.ones_like(losses))
    tiny = torch.tensor([[1e-40, 0.0], [0.0, 0.0]], dtype=torch.float64)
    reference = torch.tensor([[2e-40, 0.0], [1.0, 1.0]], dtype=torch.float64)
    assert torch.equal(
        gradient_cosine_rows(tiny, reference),
        torch.tensor([1.0, 0.0], dtype=torch.float64))

    episodes = [
        (np.array([1, 2]), np.array([3, 4, 5])),
        (np.array([6, 7]), np.array([8])),
    ]
    support, query, mask = pack_episode_indices(episodes, "cpu")
    assert support.tolist() == [[1, 2], [6, 7]]
    assert query.tolist() == [[3, 4, 5], [8, 8, 8]]
    assert mask.tolist() == [[True, True, True], [True, False, False]]


@pytest.mark.parametrize(
    "device", ["cpu"] + (["cuda"] if torch.cuda.is_available() else []))
def test_episodewise_noise_preserves_rng_mapping_when_model_is_batched(device):
    labels = torch.arange(45, dtype=torch.float32, device=device).reshape(15, 3)
    expected_generator = torch.Generator(device=device).manual_seed(901)
    actual_generator = torch.Generator(device=device).manual_seed(901)
    expected = torch.stack([
        uniform_label_noise(row, 0.2, generator=expected_generator)
        for row in labels
    ])
    actual = add_episodewise_label_noise(
        labels, 0.2, generator=actual_generator)
    assert torch.equal(actual, expected)
    assert torch.equal(actual_generator.get_state(), expected_generator.get_state())
