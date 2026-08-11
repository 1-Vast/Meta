import numpy as np
import pytest

from scripts.cache_exact_structure_proteins import (
    pack_projected_rows,
    sequence_key,
    unpack_projected_row,
)


def test_ragged_exact_projected_rows_round_trip():
    rows = [
        ("a", np.arange(12, dtype=np.float16).reshape(3, 4),
         np.zeros((128, 4), dtype=np.float16), np.ones(128, dtype=np.uint8)),
        ("b", np.arange(20, dtype=np.float16).reshape(5, 4),
         np.ones((128, 4), dtype=np.float16), np.ones(128, dtype=np.uint8)),
    ]
    packed = pack_projected_rows(rows)
    assert packed["offsets"].tolist() == [0, 3, 8]
    first = unpack_projected_row(packed, 0)
    second = unpack_projected_row(packed, 1)
    assert np.array_equal(first["exact_projected"], rows[0][1])
    assert np.array_equal(second["exact_projected"], rows[1][1])
    assert np.array_equal(second["slot_projected"], rows[1][2])


def test_exact_cache_contract_rejects_duplicate_keys_or_bad_slots():
    exact = np.zeros((2, 3), dtype=np.float16)
    slots = np.zeros((128, 3), dtype=np.float16)
    mask = np.ones(128, dtype=np.uint8)
    with pytest.raises(ValueError, match="duplicate"):
        pack_projected_rows([("x", exact, slots, mask), ("x", exact, slots, mask)])
    with pytest.raises(ValueError, match="128-slot"):
        pack_projected_rows([("x", exact, slots[:127], mask)])


def test_sequence_key_is_content_addressed_and_nonempty():
    assert sequence_key("ACDE") == sequence_key("ACDE")
    assert sequence_key("ACDE") != sequence_key("ACDF")
    with pytest.raises(ValueError, match="empty"):
        sequence_key("")
