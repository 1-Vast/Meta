import json

import numpy as np

from research.crossed_interaction.train_structure_pocket_prior_cq_observable import (
    contact_prior_features,
    fit_contact_prior_from_arrays,
    fit_structure_pocket_prior,
    load_structure_contact_labels,
    pocket_weighted_descriptor,
)


def test_fit_contact_prior_learns_separable_slot_signal():
    blocks = np.asarray([
        [0.0, 0.0],
        [0.1, 0.0],
        [1.0, 1.0],
        [1.1, 1.0],
    ], dtype=np.float64)
    labels = np.asarray([0.0, 0.0, 1.0, 1.0], dtype=np.float64)

    prior = fit_contact_prior_from_arrays(blocks, labels, ridge=0.01)
    scores = np.column_stack([
        np.ones(blocks.shape[0]),
        (blocks - prior["mean"]) / prior["scale"],
    ]) @ prior["coefficients"]

    assert prior["positive_contact_rate"] == 0.5
    assert scores[-1] > scores[0]


def test_pocket_weighted_descriptor_respects_mask_and_prior():
    blocks = np.asarray([
        [10.0, 0.0],
        [0.0, 20.0],
        [100.0, 100.0],
    ], dtype=np.float64)
    mask = np.asarray([1, 1, 0], dtype=bool)
    prior = {
        "mean": np.zeros(2, dtype=np.float64),
        "scale": np.ones(2, dtype=np.float64),
        "coefficients": np.asarray([0.0, -10.0, 10.0], dtype=np.float64),
    }

    descriptor = pocket_weighted_descriptor(blocks, mask, prior, mode="structure_prior")

    assert descriptor[1] > descriptor[0]
    assert descriptor[0] < 1.0
    assert descriptor[1] > 19.0


def test_ligand_conditioned_pocket_descriptor_changes_with_ligand():
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

    first = pocket_weighted_descriptor(
        blocks, mask, prior, ligand=np.asarray([1.0, 0.0]), mode="structure_prior")
    second = pocket_weighted_descriptor(
        blocks, mask, prior, ligand=np.asarray([0.0, 1.0]), mode="structure_prior")

    assert first[0] > 0.99
    assert second[1] > 0.99


def test_contact_prior_features_ligand_conditioned_outer_shape():
    blocks = np.ones((3, 2), dtype=np.float64)
    ligand = np.arange(4, dtype=np.float64)
    features = contact_prior_features(
        blocks, ligand, prior_feature_mode="ligand_conditioned")
    assert features.shape == (3, 8)


def test_load_structure_contact_labels_or_merges_same_sequence(tmp_path):
    supervision = tmp_path / "supervision"
    supervision.mkdir()
    records = tmp_path / "complexes.jsonl"
    records.write_text(
        '{"source_entry_id":"e1","sequence_sha256":"s1"}\n'
        '{"source_entry_id":"e2","sequence_sha256":"s1"}\n',
        encoding="utf-8")
    (supervision / "manifest.json").write_text("{}", encoding="utf-8")
    (supervision / "pairs.jsonl").write_text(
        '{"source_entry_id":"e1","shard":"shard_000000.npz","shard_index":0}\n'
        '{"source_entry_id":"e2","shard":"shard_000000.npz","shard_index":1}\n',
        encoding="utf-8")
    contact = np.zeros((2, 2, 4), dtype=np.uint8)
    contact[0, 0, 1] = 1
    contact[1, 1, 3] = 1
    residue_mask = np.ones((2, 4), dtype=np.uint8)
    np.savez_compressed(
        supervision / "shard_000000.npz",
        contact=contact,
        distance_bin=np.zeros((2, 2, 4), dtype=np.uint8),
        atom_mask=np.ones((2, 2), dtype=np.uint8),
        residue_mask=residue_mask,
    )

    labels, metadata = load_structure_contact_labels(supervision, records)

    assert labels["s1"].tolist() == [False, True, False, True]
    assert metadata["structure_pairs_used"] == 2


def test_fit_structure_pocket_prior_uses_external_bank(tmp_path):
    supervision = tmp_path / "supervision"
    protein_bank = tmp_path / "protein_bank"
    supervision.mkdir()
    protein_bank.mkdir()
    records = tmp_path / "complexes.jsonl"
    records.write_text(
        '{"source_entry_id":"e1","sequence_sha256":"s1"}\n'
        '{"source_entry_id":"e2","sequence_sha256":"s2"}\n',
        encoding="utf-8")
    (supervision / "manifest.json").write_text("{}", encoding="utf-8")
    (supervision / "pairs.jsonl").write_text(
        '{"source_entry_id":"e1","shard":"shard_000000.npz","shard_index":0}\n'
        '{"source_entry_id":"e2","shard":"shard_000000.npz","shard_index":1}\n',
        encoding="utf-8")
    contact = np.zeros((2, 2, 4), dtype=np.uint8)
    contact[0, 0, 0] = 1
    contact[1, 0, 3] = 1
    np.savez_compressed(
        supervision / "shard_000000.npz",
        contact=contact,
        distance_bin=np.zeros((2, 2, 4), dtype=np.uint8),
        atom_mask=np.ones((2, 2), dtype=np.uint8),
        residue_mask=np.ones((2, 4), dtype=np.uint8),
    )
    np.savez_compressed(
        protein_bank / "shard_000000.npz",
        keys=np.asarray(["s1", "s2"]),
        pooled=np.zeros((2, 4), dtype=np.float16),
        residues=np.asarray([
            [[1, 0, 1, 0], [0, 1, 0, 1], [0, 1, 0, 1], [0, 1, 0, 1]],
            [[0, 1, 0, 1], [0, 1, 0, 1], [0, 1, 0, 1], [1, 0, 1, 0]],
        ], dtype=np.float16),
        mask=np.ones((2, 4), dtype=np.uint8),
    )
    (protein_bank / "manifest.json").write_text(
        json.dumps({"shards": [{"path": "shard_000000.npz"}], "model_id": "test"}),
        encoding="utf-8")

    prior, metadata = fit_structure_pocket_prior(
        supervision, records, protein_bank,
        hidden_blocks=2, pocket_ridge=0.1, max_records=None)

    assert prior["slot_rows"] == 8
    assert metadata["structure_sequences_used_for_prior"] == 2
    assert metadata["hidden_blocks"] == 2


def test_fit_structure_pocket_prior_can_use_ligand_conditioning(tmp_path):
    supervision = tmp_path / "supervision"
    protein_bank = tmp_path / "protein_bank"
    supervision.mkdir()
    protein_bank.mkdir()
    records = tmp_path / "complexes.jsonl"
    records.write_text(
        '{"source_entry_id":"e1","sequence_sha256":"s1","canonical_smiles":"CCO"}\n'
        '{"source_entry_id":"e2","sequence_sha256":"s2","canonical_smiles":"CCN"}\n',
        encoding="utf-8")
    (supervision / "manifest.json").write_text("{}", encoding="utf-8")
    (supervision / "pairs.jsonl").write_text(
        '{"source_entry_id":"e1","shard":"shard_000000.npz","shard_index":0}\n'
        '{"source_entry_id":"e2","shard":"shard_000000.npz","shard_index":1}\n',
        encoding="utf-8")
    contact = np.zeros((2, 2, 4), dtype=np.uint8)
    contact[0, 0, 0] = 1
    contact[1, 0, 3] = 1
    np.savez_compressed(
        supervision / "shard_000000.npz",
        contact=contact,
        distance_bin=np.zeros((2, 2, 4), dtype=np.uint8),
        atom_mask=np.ones((2, 2), dtype=np.uint8),
        residue_mask=np.ones((2, 4), dtype=np.uint8),
    )
    np.savez_compressed(
        protein_bank / "shard_000000.npz",
        keys=np.asarray(["s1", "s2"]),
        pooled=np.zeros((2, 4), dtype=np.float16),
        residues=np.ones((2, 4, 4), dtype=np.float16),
        mask=np.ones((2, 4), dtype=np.uint8),
    )
    (protein_bank / "manifest.json").write_text(
        json.dumps({"shards": [{"path": "shard_000000.npz"}], "model_id": "test"}),
        encoding="utf-8")

    prior, metadata = fit_structure_pocket_prior(
        supervision, records, protein_bank,
        hidden_blocks=2, pocket_ridge=0.1, max_records=None,
        prior_feature_mode="ligand_conditioned")

    assert prior["prior_feature_mode"] == "ligand_conditioned"
    assert prior["mean"].shape[0] == 80
    assert metadata["prior_feature_mode"] == "ligand_conditioned"
