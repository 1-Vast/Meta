import numpy as np

from research.e0_identifiability.run_tbasis_radial import (
    bin_rbf_expectation,
    radial_basis,
    slot_composition,
)


def test_radial_basis_has_fixed_cutoff():
    value = radial_basis(np.asarray([2.0, 10.0, 12.0]))
    assert value.shape == (3, 6)
    assert value[0].max() > 0
    assert np.all(value[1:] == 0)


def test_bin_moments_match_geometry_contract():
    value = bin_rbf_expectation(points=512)
    assert value.shape == (5, 6)
    assert np.isfinite(value).all()
    assert np.all(value[-1] == 0)


def test_slot_composition_is_normalized_on_occupied_slots():
    value = slot_composition("ACDEFGHIKLMNPQRSTVWY" * 8)
    occupied = value.sum(axis=1) > 0
    assert value.shape == (128, 6)
    assert np.allclose(value[occupied].sum(axis=1), 1.0)
