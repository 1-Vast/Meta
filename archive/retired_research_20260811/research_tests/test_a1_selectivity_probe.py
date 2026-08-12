import numpy as np

from research.meta_fewshot.a1_selectivity_probe import (
    bootstrap_contrast,
    contrast_sensitivity,
    component_folds,
    component_losses,
    score,
)


def test_component_losses_average_rows_within_group_then_component():
    rows = [
        {"group": "g1", "component": "c1"},
        {"group": "g1", "component": "c1"},
        {"group": "g2", "component": "c1"},
        {"group": "g3", "component": "c2"},
    ]
    losses = component_losses(np.zeros(4), np.asarray([1.0, 1.0, 2.0, 3.0]), rows,
                              np.arange(4))
    assert losses == {"c1": 2.5, "c2": 9.0}


def test_component_folds_never_split_a_component():
    rows = [{"group": f"g{i}", "component": f"c{i // 2}"} for i in range(10)]
    folds = component_folds(rows, folds=3)
    flattened = [component for fold in folds for component in fold]
    assert len(flattened) == len(set(flattened)) == 5


def test_bootstrap_uses_component_paired_losses():
    result = bootstrap_contrast({"a": 2.0, "b": 4.0}, {"a": 1.0, "b": 1.0}, draws=999)
    assert result["mean_loss_reduction"] == 2.0
    assert result["pass"]


def test_pair_sign_ignores_truth_ties_and_scores_prediction_ties_as_half():
    rows = [
        {"group": "g", "component": "c"},
        {"group": "g", "component": "c"},
        {"group": "g", "component": "c"},
    ]
    result = score(np.zeros(3), np.asarray([-1.0, -1.0, 2.0]), rows)
    assert result["component_macro_pair_sign_accuracy"] == 0.5
    assert result["component_macro_group_pearson"] is None
    assert result["pair_sign_counts"] == {
        "comparable_truth_pairs": 2,
        "truth_tied_pairs_exact": 1,
        "prediction_tied_pairs_among_comparable": 2,
        "truth_tie_tolerance": 0.0,
    }


def test_contrast_sensitivity_reports_giant_and_leave_one_out():
    result = contrast_sensitivity(
        {"giant": 3.0, "small": 1.0}, {"giant": 1.0, "small": 2.0},
        {"giant": 9, "small": 1},
    )
    assert result["components_favoring_correct"] == 1
    assert result["components_favoring_control"] == 1
    assert result["group_weighted_mean"] == 1.7
    assert result["leave_giant_out_component_macro_mean"] == -1.0
