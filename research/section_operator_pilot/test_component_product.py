import numpy as np

from .component_product import _location


def test_location_is_label_permutation_invariant():
    surface = np.asarray([-0.2, 0.1, 0.3, -0.1, 0.5])
    residual = np.asarray([0.0, 0.2, 0.4, -0.05, 0.6])
    permutation = np.asarray([3, 0, 4, 1, 2])
    assert _location(surface, residual, 10.0) == _location(
        surface, residual[permutation], 10.0)


def test_location_is_bounded():
    surface = np.zeros(5)
    assert _location(surface, np.full(5, 100.0), 0.0) == 0.5
    assert _location(surface, np.full(5, -100.0), 0.0) == -0.5
