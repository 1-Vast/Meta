import numpy as np
import torch

from research.e0_identifiability.run_objective_design_solver import (
    _loss_and_gradient, _pinv_solution, _svd_audit,
)


def test_old_rank_objective_conflicts_with_exact_residual_teacher():
    features = np.asarray([[0.0], [1.0], [2.0]], dtype=np.float64)
    weights = np.asarray([1.0], dtype=np.float64)
    residuals = features[:, 0]
    baselines = np.asarray([3.0, 0.0, 0.0])
    labels = baselines + residuals
    left = np.asarray([0, 0, 1])
    right = np.asarray([1, 2, 2])
    _, point_gradient = _loss_and_gradient(
        "point_huber", features, weights, residuals, baselines, labels, left, right)
    _, old_gradient = _loss_and_gradient(
        "old_logistic", features, weights, residuals, baselines, labels, left, right)
    _, difference_gradient = _loss_and_gradient(
        "residual_difference_huber", features, weights, residuals, baselines,
        labels, left, right)
    assert point_gradient < 1e-12
    assert difference_gradient < 1e-12
    assert old_gradient > 1e-3


def test_pinv_recovers_exact_predictions_under_parameter_nonuniqueness():
    features = np.asarray([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    target = np.asarray([2.0, 4.0, 6.0])
    solution = _pinv_solution(features, target, 1e-10)
    assert np.allclose(features @ solution, target, atol=1e-10)
    assert np.allclose(solution, [1.0, 1.0], atol=1e-10)


def test_svd_audit_separates_rank_from_holdout_transport():
    train = np.asarray([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    holdout = np.asarray([[4.0, 4.0]])
    teacher = np.asarray([1.0, 1.0])
    audit, _, _, _ = _svd_audit(train, holdout, teacher)
    assert audit["primary_rank"] == 1
    assert audit["holdout_row_space_coverage"]["minimum"] > 1 - 1e-12
    assert audit["teacher_specific_unseen_contribution"]["maximum_absolute"] < 1e-12
