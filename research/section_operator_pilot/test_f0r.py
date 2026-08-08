import numpy as np

from .f0r_support_design import _design, _query_indices, _solve


def test_support_and_query_are_scaffold_disjoint():
    scaffold = np.asarray(["a", "a", "b", "c", "d", "e"], dtype=object)
    factor = np.arange(12, dtype=float).reshape(6, 2)
    y = np.ones((6, 1))
    support = np.asarray([0, 2])
    query = _query_indices(y, 0, support, scaffold)
    assert 0 not in query and 1 not in query and 2 not in query
    assert set(query) == {3, 4, 5}


def test_random_design_uses_distinct_scaffolds():
    scaffold = np.asarray(["a", "a", "b", "c", "d"], dtype=object)
    factor = np.arange(10, dtype=float).reshape(5, 2)
    support = _design(scaffold, factor, np.arange(5), 3, "random", 7)
    assert len({scaffold[index] for index in support}) == 3


def test_d_optimal_design_is_deterministic_and_distinct():
    scaffold = np.asarray(["a", "a", "b", "c", "d"], dtype=object)
    factor = np.asarray([[0, 0], [5, 5], [1, 0], [0, 1], [1, 1]], dtype=float)
    left = _design(scaffold, factor, np.arange(5), 3, "d_optimal", 1)
    right = _design(scaffold, factor, np.arange(5), 3, "d_optimal", 999)
    assert np.array_equal(left, right)
    assert len({scaffold[index] for index in left}) == 3


def test_positive_ridge_section_is_finite():
    factor = np.asarray([[0, 0], [1, 0], [0, 1]], dtype=float)
    coefficient = _solve(factor, np.asarray([0.2, 0.4, 0.3]),
                         np.zeros(2), 1.0, 1.0)
    assert coefficient.shape == (3,)
    assert np.isfinite(coefficient).all()
