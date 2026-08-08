import numpy as np
import pandas as pd
import torch

from .conditional_bilinear import (
    SectionConditionedBilinear,
    _d_optimal,
    _nearest_nonself,
    _subpocket_weights,
)


def test_section_is_support_permutation_invariant():
    torch.manual_seed(7)
    model = SectionConditionedBilinear([5, 4, 3], [3, 4, 5])
    support = torch.randn(2, 5, 3)
    residual = torch.randn(2, 5)
    permutation = torch.tensor([3, 0, 4, 1, 2])
    original = model.adapt(support, residual)
    permuted = model.adapt(support[:, permutation], residual[:, permutation])
    torch.testing.assert_close(original, permuted, atol=1e-7, rtol=1e-6)


def test_d_optimal_selects_distinct_scaffolds():
    scaffold = np.asarray(["a", "a", "b", "c", "d", "e"], dtype=object)
    factor = np.asarray([[0, 0], [1, 1], [1, 0], [0, 1], [-1, 0], [0, -1]], dtype=float)
    selected = _d_optimal(scaffold, factor, np.arange(6), 4, "d_optimal", 7)
    assert len(selected) == 4
    assert len({scaffold[index] for index in selected}) == 4


def test_nearest_nonself_never_returns_self():
    views = [np.asarray([[0.0], [0.1], [2.0]], dtype=np.float32)] * 3
    nearest = _nearest_nonself(views)
    assert np.all(nearest != np.arange(3))
    assert nearest.tolist()[:2] == [1, 0]


def test_kissim_masks_are_bounded_and_distinct(tmp_path):
    path = tmp_path / "min_max_distances_fine.csv"
    columns = [f"position_{index}" for index in range(85)]
    rows = []
    for offset, subpocket in enumerate(
            ("hinge_region", "dfg_region", "front_pocket")):
        low = np.linspace(0.5 + offset, 8.5 - offset, 85)
        high = low + 1.0
        rows.extend([
            {"subpocket": subpocket, "min_max": "min",
             **dict(zip(columns, low))},
            {"subpocket": subpocket, "min_max": "max",
             **dict(zip(columns, high))},
        ])
    pd.DataFrame(rows).to_csv(path, index=False)
    weights = _subpocket_weights(path)
    assert set(weights) == {"hinge_region", "dfg_region", "front_pocket"}
    assert all(value.shape == (85,) for value in weights.values())
    assert all(np.isclose(value.max(), 1.0) and value.min() > 0.0
               for value in weights.values())
    assert not np.allclose(weights["hinge_region"], weights["dfg_region"])
