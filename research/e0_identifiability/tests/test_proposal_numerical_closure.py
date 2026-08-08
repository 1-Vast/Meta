import numpy as np

from research.e0_identifiability.run_proposal_numerical_closure import (
    _augmented_least_squares, _corrected_objective,
)


def test_augmented_solver_recovers_point_and_pair_differences():
    features = np.asarray([
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [2.0, -1.0],
    ])
    teacher = np.asarray([0.75, -1.25])
    targets = features @ teacher
    left = np.asarray([0, 0, 1, 2])
    right = np.asarray([1, 2, 3, 3])
    solution = _augmented_least_squares(features, targets, left, right)
    audit = _corrected_objective(
        features, targets, solution, left, right)
    assert np.allclose(solution, teacher, atol=1e-12)
    assert audit["train_rmse"] < 1e-12
    assert audit["full_gradient_l2"] < 1e-12


def test_augmented_solver_handles_parameter_nonuniqueness():
    features = np.asarray([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    targets = np.asarray([2.0, 4.0, 6.0])
    left = np.asarray([0, 0, 1])
    right = np.asarray([1, 2, 2])
    solution = _augmented_least_squares(features, targets, left, right)
    assert np.allclose(features @ solution, targets, atol=1e-12)

