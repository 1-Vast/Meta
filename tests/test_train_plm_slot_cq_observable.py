import json

import numpy as np

from research.crossed_interaction.train_plm_slot_cq_observable import (
    load_protein_bank,
    protein_slot_descriptor,
)


def test_protein_slot_descriptor_has_segment_block_shape():
    residues = np.arange(8 * 6, dtype=np.float64).reshape(8, 6)
    mask = np.asarray([1, 1, 0, 0, 1, 1, 1, 1], dtype=bool)
    descriptor = protein_slot_descriptor(
        residues, mask, slot_segments=4, hidden_blocks=3)
    assert descriptor.shape == (12,)
    assert np.isfinite(descriptor).all()


def test_protein_slot_descriptor_zeroes_empty_segments():
    residues = np.ones((8, 4), dtype=np.float64)
    mask = np.asarray([1, 1, 0, 0, 1, 1, 1, 1], dtype=bool)
    descriptor = protein_slot_descriptor(
        residues, mask, slot_segments=4, hidden_blocks=2)
    assert np.array_equal(descriptor[2:4], np.zeros(2))


def test_protein_slot_descriptor_is_region_sensitive():
    first = np.zeros((8, 4), dtype=np.float64)
    second = np.zeros((8, 4), dtype=np.float64)
    first[:2, :] = 1.0
    second[-2:, :] = 1.0
    mask = np.ones(8, dtype=bool)
    assert not np.array_equal(
        protein_slot_descriptor(first, mask, slot_segments=4, hidden_blocks=2),
        protein_slot_descriptor(second, mask, slot_segments=4, hidden_blocks=2))


def test_load_protein_bank_reads_required_keys(tmp_path):
    bank = tmp_path
    np.savez_compressed(
        bank / "shard_000000.npz",
        keys=np.asarray(["p1", "p2"]),
        pooled=np.zeros((2, 4), dtype=np.float16),
        residues=np.zeros((2, 3, 4), dtype=np.float16),
        mask=np.ones((2, 3), dtype=np.uint8),
    )
    manifest = {"shards": [{"path": "shard_000000.npz"}], "model_id": "test"}
    (bank / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    rows, loaded_manifest = load_protein_bank(bank, {"p2"})
    assert sorted(rows) == ["p2"]
    assert rows["p2"]["residues"].shape == (3, 4)
    assert loaded_manifest["model_id"] == "test"
