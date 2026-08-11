import numpy as np
import pytest

from scripts.cache_r0b_exact_geometry import (
    _validate_geometry_arrays,
    pack_exact_geometry_rows,
    unpack_exact_geometry_row,
)


def _row(entry, atoms, residues):
    metadata = {"source_entry_id": entry, "atoms": atoms, "residues": residues}
    distance = np.arange(atoms * residues, dtype=np.float32).reshape(atoms, residues)
    arrays = {
        "distance_angstrom": distance,
        "distance_bin": np.digitize(
            distance, [4.0, 6.0, 8.0, 10.0]).astype(np.uint8),
        "slot_of_residue": np.arange(residues, dtype=np.uint8),
    }
    return metadata, arrays


def test_exact_geometry_ragged_shard_round_trip():
    rows = [_row("a", 2, 3), _row("b", 3, 4)]
    payload = pack_exact_geometry_rows(rows)
    assert payload["pair_offsets"].tolist() == [0, 6, 18]
    assert payload["residue_offsets"].tolist() == [0, 3, 7]
    for index, expected in enumerate(rows):
        actual = unpack_exact_geometry_row(payload, index)
        assert np.array_equal(actual["distance_angstrom"],
                              expected[1]["distance_angstrom"])
        assert np.array_equal(actual["distance_bin"], expected[1]["distance_bin"])
        assert np.array_equal(actual["slot_of_residue"],
                              expected[1]["slot_of_residue"])


def test_exact_geometry_unpack_rejects_bad_offsets():
    payload = pack_exact_geometry_rows([_row("a", 2, 3)])
    payload["pair_offsets"][-1] -= 1
    with pytest.raises(ValueError, match="pair offsets"):
        unpack_exact_geometry_row(payload, 0)


def test_exact_geometry_rejects_shape_and_bin_contract_mismatch():
    metadata, arrays = _row("a", 2, 3)
    broken = dict(arrays)
    broken["distance_angstrom"] = arrays["distance_angstrom"][:, :2]
    with pytest.raises(ValueError, match="dimensions"):
        _validate_geometry_arrays(metadata, broken)

    broken = dict(arrays)
    broken["distance_bin"] = arrays["distance_bin"].copy()
    broken["distance_bin"][0, 0] = 4
    with pytest.raises(ValueError, match="distance bins"):
        _validate_geometry_arrays(metadata, broken)
