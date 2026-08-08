import numpy as np

from .gauge_orthogonal import _curve, _solve_orthogonal


def test_orthogonal_section_is_pair_permutation_invariant():
    surface = np.asarray([-0.2, 0.1, 0.3, -0.1, 0.5])
    residual = np.asarray([0.0, 0.2, 0.4, -0.05, 0.6])
    permutation = np.asarray([3, 0, 4, 1, 2])
    reference = _solve_orthogonal(surface, residual, 10.0, 1.0)
    permuted = _solve_orthogonal(surface[permutation], residual[permutation], 10.0, 1.0)
    np.testing.assert_allclose(reference, permuted, atol=1e-12)


def test_coordinates_obey_frozen_bounds():
    surface = np.asarray([-1.0, 0.0, 1.0])
    residual = np.asarray([100.0, -100.0, 100.0])
    coordinate = _solve_orthogonal(surface, residual, 0.01, 0.01)
    assert -0.5 <= coordinate[0] <= 0.5
    assert 0.0 <= coordinate[1] <= 2.0
    assert np.all(np.isfinite(_curve(surface, coordinate)))
