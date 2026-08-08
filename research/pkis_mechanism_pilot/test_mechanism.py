import numpy as np

from research.pkis_mechanism_pilot.mechanism import (
    bounded_biological_z,
    channel_contributions,
    deterministic_derangement,
    double_center,
    fit_channel_bounds,
    ligand_from_smiles,
    ordered_anchor_simplex,
    pair_feature_matrix,
    protein_pair_properties,
    TargetRecord,
)


def test_ligand_statistic_is_fixed_shape_and_finite():
    record = ligand_from_smiles("x", "CCOc1ccccc1N")
    assert record.model_features.shape == (1036,)
    assert record.pharmacophore_shells.shape == (7, 4)
    assert np.isfinite(record.model_features).all()
    assert np.isfinite(record.pharmacophore_shells).all()


def test_pair_tensor_and_five_bounded_coordinates():
    ligand = ligand_from_smiles("x", "CCOc1ccccc1N")
    pocket = "A" * 85
    features = pair_feature_matrix(
        ligand.pharmacophore_shells[None], protein_pair_properties(pocket)[None]
    )
    assert features.shape == (1, 1, 5 * 4 * 85)
    coefficients = np.linspace(-0.5, 0.5, features.shape[-1])
    contribution = channel_contributions(features, coefficients)
    bounds = fit_channel_bounds(np.concatenate([contribution - 1, contribution + 1], axis=0))
    z = bounded_biological_z(contribution, bounds)
    assert z.shape == (1, 1, 5)
    assert np.isfinite(z).all()
    assert np.all((z >= 0.0) & (z <= 1.0))


def test_double_center_removes_both_additive_axes():
    matrix = np.arange(30, dtype=float).reshape(5, 6) ** 1.3
    centered = double_center(matrix)
    assert np.max(np.abs(centered.mean(axis=0))) < 1e-12
    assert np.max(np.abs(centered.mean(axis=1))) < 1e-12


def test_derangement_has_no_fixed_points_and_prefers_group():
    records = [
        TargetRecord(i, f"n{i}", f"g{i}", "f", "A" if i < 3 else "B", "u", "A" * 85)
        for i in range(6)
    ]
    mapping = deterministic_derangement(records)
    assert np.all(mapping != np.arange(6))
    assert all(records[i].group == records[int(mapping[i])].group for i in range(6))


def test_law_interface_is_exact_simplex_without_scalar_output():
    weights = ordered_anchor_simplex(np.asarray([0.0, 0.2, 0.5, 1.0]),
                                     np.asarray([0.0, 0.2, 0.5, 1.0]))
    assert weights.shape == (4, 8)
    assert np.min(weights) >= 0.0
    assert np.max(np.abs(weights.sum(axis=1) - 1.0)) < 1e-12
    assert weights[-1, 7] == 1.0

