from research.e0_identifiability.run_tdir_pilot import (
    _collect_ligand_serials,
    _complex_weights,
    _cross_split_overlap,
)


def test_collect_ligand_serials_filters_protein_indices():
    from collections import namedtuple

    Charge = namedtuple("Charge", ["atoms_orig_idx", "type"])
    value = {
        "atoms_orig_idx": [11, 12],
        "protein_orig_idx": 900,
        "nested": Charge(atoms_orig_idx=[13], type="negative"),
    }
    assert _collect_ligand_serials(value, {11, 12, 13}) == {11, 12, 13}


def test_complex_weights_give_each_complex_equal_mass():
    import numpy as np

    identifiers = np.asarray(["a", "a", "a", "b"])
    weights = _complex_weights(identifiers)
    assert np.isclose(weights[identifiers == "a"].sum(), weights[identifiers == "b"].sum())
    assert np.isclose(weights.mean(), 1.0)


def test_cross_split_overlap_reports_shared_value():
    records = [
        {"source_split": "train", "value": "x"},
        {"source_split": "test", "value": "x"},
        {"source_split": "test", "value": "y"},
    ]
    assert _cross_split_overlap(records, "value") == [
        {"value": "x", "splits": ["test", "train"]}
    ]
