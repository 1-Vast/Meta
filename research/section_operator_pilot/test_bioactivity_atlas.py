import numpy as np

from .bioactivity_atlas import (
    _anchor_weights,
    _apply_section,
    _atlas_surfaces,
    _solve_section,
)


def test_anchor_weights_are_sparse_and_normalized():
    source = [np.arange(20, dtype=float)[:, None]] * 3
    query = [np.asarray([[2.2], [15.7]])] * 3
    weights = _anchor_weights(query, source, [1.0, 1.0, 1.0], n_anchor=4)
    for value in weights:
        np.testing.assert_allclose(value.sum(axis=1), 1.0)
        assert np.all(np.count_nonzero(value, axis=1) == 4)


def test_atlas_surfaces_are_task_centered():
    rng = np.random.default_rng(3)
    profile = rng.normal(size=(11, 7))
    source = [rng.normal(size=(7, 3)) for _ in range(3)]
    query = [rng.normal(size=(5, 3)) for _ in range(3)]
    surface = _atlas_surfaces(profile, query, source, [1.0, 1.0, 1.0])
    assert surface.shape == (11, 5, 3)
    np.testing.assert_allclose(surface.mean(axis=0), 0.0, atol=1e-12)


def test_ridge_section_is_support_permutation_invariant():
    rng = np.random.default_rng(4)
    surface = rng.normal(size=(5, 3))
    residual = rng.normal(size=5)
    permutation = np.asarray([2, 4, 0, 3, 1])
    coefficient = _solve_section(surface, residual, 10.0, 0.1)
    permuted = _solve_section(
        surface[permutation], residual[permutation], 10.0, 0.1)
    np.testing.assert_allclose(coefficient, permuted, atol=1e-12)
    np.testing.assert_allclose(
        _apply_section(surface, coefficient),
        _apply_section(surface, permuted), atol=1e-12)
