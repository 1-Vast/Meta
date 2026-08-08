import numpy as np

from .counterfactual_design import _counterfactual_support


def test_counterfactual_design_is_deterministic_and_scaffold_unique():
    rng = np.random.default_rng(4)
    correct = rng.normal(size=(12, 3))
    wrong = rng.normal(size=(12, 3))
    scaffold = np.asarray([f"s{i // 2}" for i in range(12)], dtype=object)
    first = _counterfactual_support(scaffold, correct, wrong, np.arange(12), 5)
    second = _counterfactual_support(scaffold, correct, wrong, np.arange(12), 5)
    np.testing.assert_array_equal(first, second)
    assert len({scaffold[index] for index in first}) == 5


def test_counterfactual_design_requires_enough_scaffolds():
    surface = np.zeros((5, 3))
    scaffold = np.asarray(["a", "a", "b", "b", "c"], dtype=object)
    assert _counterfactual_support(scaffold, surface, surface, np.arange(5), 4) is None
