import numpy as np

from .component_independent import _task_location


def test_task_location_is_protein_independent_and_permutation_invariant():
    residual = np.asarray([0.0, 0.2, 0.4, -0.05, 0.6])
    permutation = np.asarray([3, 0, 4, 1, 2])
    reference = _task_location(residual, 10.0)
    assert reference == _task_location(residual[permutation], 10.0)
    # There is deliberately no protein/surface argument.
    assert reference == _task_location(residual.copy(), 10.0)


def test_task_location_is_bounded():
    assert _task_location(np.full(5, 100.0), 0.0) == 0.5
    assert _task_location(np.full(5, -100.0), 0.0) == -0.5
