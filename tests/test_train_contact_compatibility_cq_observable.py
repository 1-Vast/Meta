import numpy as np

from research.crossed_interaction.train_contact_compatibility_cq_observable import (
    contact_compatibility_descriptor,
)


def test_contact_compatibility_descriptor_has_fixed_shape():
    blocks = np.asarray([
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
    ], dtype=np.float64)
    mask = np.asarray([1, 1, 1], dtype=bool)
    ligand = np.asarray([1.0, 0.0], dtype=np.float64)
    prior = {
        "mean": np.zeros(4, dtype=np.float64),
        "scale": np.ones(4, dtype=np.float64),
        "coefficients": np.asarray([0.0, 5.0, 0.0, -5.0, 0.0], dtype=np.float64),
        "prior_feature_mode": "ligand_conditioned",
    }

    descriptor = contact_compatibility_descriptor(blocks, mask, ligand, prior)

    assert descriptor.shape == (17,)
    assert np.isfinite(descriptor).all()
    assert descriptor[3] >= descriptor[0]


def test_contact_compatibility_descriptor_changes_with_ligand():
    blocks = np.asarray([
        [1.0, 0.0],
        [0.0, 1.0],
    ], dtype=np.float64)
    mask = np.asarray([1, 1], dtype=bool)
    prior = {
        "mean": np.zeros(4, dtype=np.float64),
        "scale": np.ones(4, dtype=np.float64),
        "coefficients": np.asarray([0.0, 10.0, -10.0, -10.0, 10.0], dtype=np.float64),
        "prior_feature_mode": "ligand_conditioned",
    }

    first = contact_compatibility_descriptor(
        blocks, mask, np.asarray([1.0, 0.0]), prior)
    second = contact_compatibility_descriptor(
        blocks, mask, np.asarray([0.0, 1.0]), prior)

    assert not np.array_equal(first, second)


def test_contact_compatibility_descriptor_zeroes_empty_mask():
    blocks = np.ones((3, 2), dtype=np.float64)
    mask = np.zeros(3, dtype=bool)
    ligand = np.ones(2, dtype=np.float64)
    prior = {
        "mean": np.zeros(4, dtype=np.float64),
        "scale": np.ones(4, dtype=np.float64),
        "coefficients": np.zeros(5, dtype=np.float64),
        "prior_feature_mode": "ligand_conditioned",
    }

    descriptor = contact_compatibility_descriptor(blocks, mask, ligand, prior)

    assert np.array_equal(descriptor, np.zeros(17, dtype=np.float64))
