import json

import numpy as np

from research.crossed_interaction.train_pocket_prototype_cq_observable import (
    build_pocket_prototypes,
    farthest_first_indices,
    ligand_prototype_similarity,
    load_structure_pocket_examples,
    protein_prototype_similarity,
)


def test_farthest_first_indices_returns_spread_subset():
    features = np.asarray([
        [1.0, 0.0],
        [0.9, 0.1],
        [-1.0, 0.0],
        [0.0, 1.0],
    ], dtype=np.float64)
    indices = farthest_first_indices(features, count=2, seed=1)
    assert indices.shape == (2,)
    assert len(set(indices.tolist())) == 2


def test_protein_prototype_similarity_uses_best_active_slot():
    blocks = np.asarray([
        [1.0, 0.0],
        [0.0, 1.0],
        [10.0, 10.0],
    ], dtype=np.float64)
    mask = np.asarray([1, 1, 0], dtype=bool)
    prototypes = np.asarray([
        [1.0, 0.0],
        [0.0, 1.0],
    ], dtype=np.float64)
    similarity = protein_prototype_similarity(blocks, mask, prototypes)
    assert np.allclose(similarity, [1.0, 1.0])


def test_ligand_prototype_similarity_is_cosine():
    ligand = np.asarray([1.0, 0.0], dtype=np.float64)
    prototypes = np.asarray([
        [1.0, 0.0],
        [0.0, 1.0],
    ], dtype=np.float64)
    similarity = ligand_prototype_similarity(ligand, prototypes)
    assert np.allclose(similarity, [1.0, 0.0])


def test_load_structure_pocket_examples_reads_contact_and_ligand(tmp_path):
    supervision = tmp_path / "supervision"
    supervision.mkdir()
    records = tmp_path / "complexes.jsonl"
    records.write_text(
        '{"source_entry_id":"e1","sequence_sha256":"s1","canonical_smiles":"CCO"}\n',
        encoding="utf-8")
    (supervision / "manifest.json").write_text("{}", encoding="utf-8")
    (supervision / "pairs.jsonl").write_text(
        '{"source_entry_id":"e1","shard":"shard_000000.npz","shard_index":0}\n',
        encoding="utf-8")
    contact = np.zeros((1, 2, 4), dtype=np.uint8)
    contact[0, 0, 2] = 1
    np.savez_compressed(
        supervision / "shard_000000.npz",
        contact=contact,
        distance_bin=np.zeros((1, 2, 4), dtype=np.uint8),
        atom_mask=np.ones((1, 2), dtype=np.uint8),
        residue_mask=np.ones((1, 4), dtype=np.uint8),
    )

    examples, metadata = load_structure_pocket_examples(supervision, records)

    assert len(examples) == 1
    assert examples[0]["label"].tolist() == [False, False, True, False]
    assert examples[0]["ligand"].shape == (40,)
    assert metadata["structure_pairs_used"] == 1


def test_build_pocket_prototypes_selects_external_pockets(tmp_path):
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
        residues=np.asarray([
            [[1, 0, 1, 0], [0, 1, 0, 1], [0, 1, 0, 1], [0, 1, 0, 1]],
            [[0, 1, 0, 1], [0, 1, 0, 1], [0, 1, 0, 1], [1, 0, 1, 0]],
        ], dtype=np.float16),
        mask=np.ones((2, 4), dtype=np.uint8),
    )
    (protein_bank / "manifest.json").write_text(
        json.dumps({"shards": [{"path": "shard_000000.npz"}], "model_id": "test"}),
        encoding="utf-8")

    prototypes, metadata = build_pocket_prototypes(
        supervision, records, protein_bank,
        hidden_blocks=2, prototype_count=2, max_records=None,
        ligand_weight=1.0, seed=1)

    assert prototypes["pockets"].shape == (2, 2)
    assert prototypes["ligands"].shape == (2, 40)
    assert metadata["candidate_pockets"] == 2
