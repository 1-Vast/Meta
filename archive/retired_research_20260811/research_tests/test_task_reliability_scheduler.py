import numpy as np
import torch

from model.metasieve_v1 import MetaSieveV1
from research.meta_fewshot.audit_task_reliability_scheduler import (
    episode_view,
    permuted_labels,
    task_covariates,
)
from research.meta_fewshot.task_reliability_scheduler import (
    apply_cross_fitted,
    component_bootstrap,
    cross_fitted_ridge,
    permute_informative_rows,
    stable_fold_assignments,
)


def test_group_folds_are_deterministic_and_never_split_a_component():
    groups = np.asarray(["a", "a", "b", "c", "c", "d", "e", "f"])
    left = stable_fold_assignments(groups, folds=3)
    right = stable_fold_assignments(groups, folds=3)
    assert np.array_equal(left, right)
    for group in np.unique(groups):
        assert len(np.unique(left[groups == group])) == 1


def test_cross_fitted_scorer_only_predicts_rows_held_from_its_group_fit():
    rng = np.random.default_rng(7)
    groups = np.repeat([f"g{i}" for i in range(15)], 2)
    features = rng.normal(size=(len(groups), 4))
    target = features[:, 0] - 0.3 * features[:, 1]
    prediction, folds = cross_fitted_ridge(features, target, groups, folds=5)
    assert np.isfinite(prediction).all()
    assert sorted(np.concatenate([fold.test_indices for fold in folds]).tolist()) == list(
        range(len(groups)))
    transformed = features.copy()
    transformed[:, 0] *= -1
    assert apply_cross_fitted(folds, transformed).shape == prediction.shape


def test_equal_capacity_null_only_permutes_informative_columns_within_strata():
    features = np.arange(80, dtype=float).reshape(10, 8)
    size = np.repeat([1.0, 2.0], 5)
    permuted = permute_informative_rows(features, size, seed=11)
    assert np.array_equal(permuted[:, 2:], features[:, 2:])
    for value in np.unique(size):
        rows = size == value
        assert sorted(map(tuple, permuted[rows, :2])) == sorted(
            map(tuple, features[rows, :2]))


def test_component_bootstrap_keeps_destructive_controls_evaluation_only():
    rng = np.random.default_rng(13)
    groups = np.repeat([f"g{i}" for i in range(20)], 2)
    utility = rng.normal(size=len(groups))
    score = utility + 0.05 * rng.normal(size=len(groups))
    null = rng.normal(size=len(groups))
    controls = {
        "intercept_only": rng.normal(size=len(groups)),
        "label_permutation": rng.normal(size=len(groups)),
        "ligand_only": rng.normal(size=len(groups)),
        "protein_shuffle": rng.normal(size=len(groups)),
        "wrong_support": rng.normal(size=len(groups)),
    }
    control_utility = {name: utility - 0.1 for name in controls}
    result = component_bootstrap(
        score, utility, null, controls, control_utility,
        rng.normal(size=(len(groups), 3)), groups, seed=17, draws=99)
    assert result["components"] == 20
    assert result["observed"]["matched_null_advantage"] > 0
    assert set(result["observed"]["score_correlation_losses"]) == set(controls)


def test_empty_scaffold_is_missing_and_never_a_shared_familiarity_value():
    cells = [
        {"ligand_id": "known", "pK": 7.0, "replicate_count": 1,
         "panel_count": 1},
        {"ligand_id": "empty", "pK": 8.0, "replicate_count": 1,
         "panel_count": 1},
    ]
    tasks = {"a": np.asarray([0]), "b": np.asarray([1])}
    values, names, metadata = task_covariates(
        cells, tasks, {"ring"}, {"known": "ring", "empty": ""})
    overlap = names.index("scaffold_overlap_meta_val_nonempty")
    missing = names.index("scaffold_missing_fraction")
    assert values[0, overlap] == 1.0
    assert values[1, overlap] == 1.0  # imputed only for nuisance regression
    assert values[:, missing].tolist() == [0.0, 1.0]
    assert metadata["all_scaffolds_missing_tasks"] == 1
    assert metadata["enters_scheduler_scorer"] is False


def test_label_override_handles_chunks_with_different_query_widths():
    model = MetaSieveV1(
        input_dim=4, section_dim=1, ridge=1.0,
        support_only_section=True)
    tensors = {
        "y": torch.linspace(0.0, 1.0, 8),
        "ligand": torch.randn(8, 4),
        "correct": torch.randn(8, 4),
    }
    episodes = [
        (np.asarray([0]), np.asarray([1])),
        (np.asarray([2]), np.asarray([3, 4, 5])),
    ]
    support_y, query_y = permuted_labels(episodes, tensors)
    difficulty, alignment, gradients = episode_view(
        model, tensors, episodes, family="correct",
        support_override=support_y, query_override=query_y, chunk_size=1)
    assert difficulty.shape == (2,)
    assert alignment.shape == (2,)
    assert gradients.shape[0] == 2
