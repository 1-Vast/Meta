import numpy as np

from research.crossed_interaction.generate_tbasis_features import (
    apply_frozen_calibration,
    ligand_channels_smiles,
)


def test_smiles_channels_follow_heavy_atom_order():
    channels = ligand_channels_smiles("CC(=O)O")
    assert channels.shape == (4, 8)
    assert np.all(channels.sum(axis=1) > 0)


def test_frozen_calibration_shape_and_order():
    raw = np.arange(2 * 3 * 2, dtype=np.float64).reshape(2, 3, 2)
    calibration = {
        "coef": np.eye(2),
        "intercept": np.zeros(2),
        "mean": np.zeros(12),
        "scale": np.ones(12),
        "active": np.ones(12, dtype=bool),
    }
    assert np.array_equal(apply_frozen_calibration(raw, calibration), raw.reshape(-1))
