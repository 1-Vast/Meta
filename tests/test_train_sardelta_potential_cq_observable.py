import numpy as np

from research.crossed_interaction.train_sardelta_potential_cq_observable import (
    sardelta_potential,
)


def test_sardelta_potential_returns_zero_without_neighbors():
    target = np.asarray([1.0, 2.0], dtype=np.float64)
    ligand = np.asarray([0.5, 0.25], dtype=np.float64)
    model = {
        "mean": np.zeros(4, dtype=np.float64),
        "scale": np.ones(4, dtype=np.float64),
        "weights": np.ones(4, dtype=np.float64),
        "y_mean": 0.0,
    }

    value = sardelta_potential(
        target, ligand, [], model, feature_mode="delta_target")

    assert value.shape == (1,)
    assert value[0] == 0.0


def test_sardelta_potential_uses_target_and_neighbor_delta():
    target = np.asarray([1.0, 2.0], dtype=np.float64)
    ligand = np.asarray([3.0, 5.0], dtype=np.float64)
    neighbors = [np.asarray([1.0, 4.0], dtype=np.float64)]
    model = {
        "mean": np.zeros(4, dtype=np.float64),
        "scale": np.ones(4, dtype=np.float64),
        "weights": np.ones(4, dtype=np.float64),
        "y_mean": 0.0,
    }

    value = sardelta_potential(
        target, ligand, neighbors, model, feature_mode="delta_target")

    assert value[0] == 6.0
