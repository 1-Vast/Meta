import numpy as np

from research.crossed_interaction.train_sardelta_edge_cq_observable import (
    cell_edge_feature,
    edge_summary,
)


def test_edge_summary_reports_distribution_statistics():
    summary = edge_summary([1.0, 2.0, 3.0])

    assert summary.shape == (8,)
    assert summary[0] == 2.0
    assert summary[2] == 1.0
    assert summary[3] == 3.0


def test_cell_edge_feature_returns_zero_without_panel_neighbors():
    target = np.asarray([1.0, 2.0], dtype=np.float64)
    ligand = np.asarray([3.0, 5.0], dtype=np.float64)
    model = {
        "mean": np.zeros(4, dtype=np.float64),
        "scale": np.ones(4, dtype=np.float64),
        "weights": np.ones(4, dtype=np.float64),
        "y_mean": 0.0,
    }

    feature = cell_edge_feature(
        target, ligand, [ligand], model, feature_mode="delta_target")

    assert feature.shape == (8,)
    assert np.allclose(feature, 0.0)


def test_cell_edge_feature_uses_panel_neighbor_predictions():
    target = np.asarray([1.0, 2.0], dtype=np.float64)
    ligand = np.asarray([3.0, 5.0], dtype=np.float64)
    neighbor = np.asarray([1.0, 4.0], dtype=np.float64)
    model = {
        "mean": np.zeros(4, dtype=np.float64),
        "scale": np.ones(4, dtype=np.float64),
        "weights": np.ones(4, dtype=np.float64),
        "y_mean": 0.0,
    }

    feature = cell_edge_feature(
        target, ligand, [neighbor], model, feature_mode="delta_target")

    assert feature[0] == 6.0
    assert feature[1] == 0.0
