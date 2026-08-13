from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import torch

from scripts.qpsmp_data import EpisodeSpec, QPSMPData
from scripts.train_qpsmp import component_target_mean, evaluate


ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"


def governed_data():
    return QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK)


def test_governed_banks_cover_main_v0_and_hard_components_are_disjoint():
    data = governed_data()

    assert len(data.protein_bank) == 499
    assert len(data.ligand_bank) == 9880
    split_components = [set(data.components[split])
                        for split in ("meta_train", "meta_val", "meta_test")]
    assert split_components[0].isdisjoint(split_components[1])
    assert split_components[0].isdisjoint(split_components[2])
    assert split_components[1].isdisjoint(split_components[2])


def test_fixed_episode_bank_is_reproducible_and_support_query_disjoint():
    data = governed_data()
    left = data.fixed_episode_bank("meta_val", 5, 3, 1, 17)
    right = data.fixed_episode_bank("meta_val", 5, 3, 1, 17)

    assert left == right
    assert {spec.target for spec in left} == {
        target for target, indices in data.tasks["meta_val"].items() if len(indices) >= 6}
    assert all(set(spec.support).isdisjoint(spec.query) for spec in left)
    assert all(spec.target != spec.donor_target for spec in left)


def test_nested_episode_banks_share_query_and_use_support_prefixes():
    data = governed_data()
    banks = data.fixed_nested_episode_banks(
        "meta_val", (1, 2, 3, 5), 3, 1, 29, 1)

    assert {key: len(value) for key, value in banks.items()} == {
        1: len(banks[5]), 2: len(banks[5]), 3: len(banks[5]), 5: len(banks[5])}
    for episodes in zip(*(banks[size] for size in (1, 2, 3, 5))):
        largest = episodes[-1]
        assert all(item.query == largest.query for item in episodes)
        assert all(item.support == largest.support[:size]
                   for size, item in zip((1, 2, 3, 5), episodes))
        assert not set(largest.support) & set(largest.query)
        support_ligands = {data.cells[index]["ligand_id"] for index in largest.support}
        query_ligands = {data.cells[index]["ligand_id"] for index in largest.query}
        assert support_ligands.isdisjoint(query_ligands)


def test_materialized_episode_uses_raw_biological_banks_and_pki_labels():
    data = governed_data()
    spec = data.fixed_episode_bank("meta_val", 2, 2, 1, 19)[0]
    episode = data.materialize(spec)

    assert episode.protein_pooled.shape == (640,)
    assert episode.protein_tokens.shape == (128, 640)
    assert episode.support_atoms.shape == (2, 128, 32)
    assert episode.support_bonds.shape == (2, 128, 128, 12)
    assert episode.query_atoms.shape == (2, 128, 32)
    assert torch.equal(
        episode.support_y,
        torch.tensor([data.cells[index]["pK"] for index in spec.support], dtype=torch.float32))
    ligand = data.cells[spec.support[0]]["ligand_id"]
    atoms, bonds, mask = data.ligand_bank.get(ligand)
    assert torch.equal(episode.support_atoms[0], torch.from_numpy(atoms))
    assert torch.equal(episode.support_bonds[0], torch.from_numpy(bonds))
    assert torch.equal(episode.support_mask[0], torch.from_numpy(mask).bool())


def test_materialize_rejects_support_query_overlap():
    data = governed_data()
    spec = data.fixed_episode_bank("meta_val", 2, 2, 1, 23)[0]
    invalid = replace(spec, query=(spec.support[0],))

    try:
        data.materialize(invalid)
    except ValueError as error:
        assert "overlap" in str(error)
    else:
        raise AssertionError("overlapping support/query was accepted")


def test_materialize_rejects_ligand_identity_overlap():
    data = governed_data()
    spec = data.fixed_episode_bank("meta_val", 2, 2, 1, 31)[0]
    query_index = spec.query[0]
    original = data.cells[query_index]["ligand_id"]
    data.cells[query_index]["ligand_id"] = data.cells[spec.support[0]]["ligand_id"]
    try:
        with pytest.raises(ValueError, match="ligand identities overlap"):
            data.materialize(spec)
    finally:
        data.cells[query_index]["ligand_id"] = original


def test_component_target_mean_does_not_weight_large_targets_or_components_more():
    rows = [
        {"component": "a", "target": "a1", "loss": 0.0},
        {"component": "a", "target": "a1", "loss": 0.0},
        {"component": "a", "target": "a2", "loss": 2.0},
        {"component": "b", "target": "b1", "loss": 10.0},
    ]

    assert np.isclose(component_target_mean(rows, "loss"), 5.5)
