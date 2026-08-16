import numpy as np

from research.meta_fewshot.pcsar_oracle_gate import (
    apply_basis,
    component_contrast,
    draw_episode,
    fit_pca_basis,
    fit_ridge,
    predict,
)


def test_draw_episode_is_deterministic_and_disjoint():
    indices = np.arange(12)

    support_a, query_a = draw_episode(indices, seed=7, k=5, max_query=3)
    support_b, query_b = draw_episode(indices, seed=7, k=5, max_query=3)

    assert np.array_equal(support_a, support_b)
    assert np.array_equal(query_a, query_b)
    assert len(set(support_a) & set(query_a)) == 0
    assert len(query_a) == 3


def test_fit_ridge_predicts_simple_linear_signal():
    x = np.asarray([[0.0], [1.0], [2.0], [3.0]], dtype=np.float64)
    y = np.asarray([0.0, 1.0, 2.0, 3.0], dtype=np.float64)

    model = fit_ridge(x, y, ridge=0.01)
    prediction = predict(model, x)

    assert np.square(y - prediction).mean() < 1e-3


def test_fit_pca_basis_reduces_to_requested_dimension():
    x = np.asarray([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=np.float64)

    model = fit_pca_basis(x, dim=2)
    projected = apply_basis(model, x)

    assert projected.shape == (3, 2)


def test_component_contrast_uses_target_macro_direction():
    rows = []
    for target in ("t1", "t2", "t3"):
        rows.append({"target_id": target, "arm": "ORACLE", "squared_error": 0.25})
        rows.append({"target_id": target, "arm": "GLOBAL", "squared_error": 1.0})

    result = component_contrast(rows, "ORACLE", "GLOBAL", draws=999, seed=1)

    assert result["targets"] == 3
    assert result["target_macro_reduction"] == 0.75
    assert result["pass"] is True
