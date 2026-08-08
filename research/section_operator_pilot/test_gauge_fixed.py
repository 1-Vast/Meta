import numpy as np

from .gauge_fixed import _curve, _solve_gain


def test_gain_section_is_pair_permutation_invariant():
    surface = np.asarray([-0.2, 0.1, 0.3, -0.1, 0.5])
    residual = np.asarray([0.0, 0.2, 0.4, -0.05, 0.6])
    permutation = np.asarray([3, 0, 4, 1, 2])
    reference = _solve_gain(surface, residual, 10.0, 1.0)
    permuted = _solve_gain(
        surface[permutation], residual[permutation], 10.0, 1.0)
    np.testing.assert_allclose(reference, permuted, atol=1e-12)


def test_large_penalties_recover_frozen_prior():
    surface = np.asarray([-1.0, 0.0, 1.0])
    residual = np.asarray([1.0, -1.0, 0.5])
    coefficient = _solve_gain(surface, residual, 1e12, 1e12)
    np.testing.assert_allclose(coefficient, [0.0, 1.0], atol=1e-10)
    np.testing.assert_allclose(_curve(surface, coefficient), surface, atol=1e-10)
