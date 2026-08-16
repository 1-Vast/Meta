from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import torch

from research.a2s.a2s_cfes_semantic_gate import (
    CONTACT_TYPES,
    SAFE_SPLIT_COLUMNS,
    AdditiveContactModel,
    BilinearResidual,
    NoCrossResidual,
    assert_allowed_splits,
    assert_safe_columns,
    contact_profile,
    dual_cold_indices,
    matched_donors,
    normalized_state_pool,
    parameter_count,
    parse_fasta,
    pocket_features,
)


def test_firewall_rejects_affinity_and_test_rows() -> None:
    assert "system_proper_pocket_num_residues" in SAFE_SPLIT_COLUMNS
    assert "system_proper_num_pocket_residues" not in SAFE_SPLIT_COLUMNS
    assert_safe_columns(("system_id",), ("system_id",))
    with pytest.raises(AssertionError):
        assert_safe_columns(("ligand_binding_affinity",), ("ligand_binding_affinity",))
    assert_allowed_splits(pd.DataFrame({"split": ["train", "val"]}))
    with pytest.raises(AssertionError):
        assert_allowed_splits(pd.DataFrame({"split": ["train", "test"]}))


def test_fasta_parser_retains_accessions_and_sequences() -> None:
    text = ">sp|P00390|ONE description\nACDE\nFG\n>tr|Q00001-2|TWO\nHIKL\n"
    assert parse_fasta(text) == {"P00390": "ACDEFG", "Q00001-2": "HIKL"}


def test_contact_profile_uses_registered_types_only() -> None:
    profile = contact_profile(
        [
            "1.A_4_type:hydrogen_bonds__protisdon:True",
            "1.A_5_type:hydrogen_bonds__protisdon:False",
            "1.A_8_type:pi_stacks__stack_type:T",
            "1.A_9_type:unregistered_contact",
        ]
    )
    assert profile.shape == (len(CONTACT_TYPES),)
    assert profile[CONTACT_TYPES.index("hydrogen_bonds")] == pytest.approx(np.log1p(2))
    assert profile[CONTACT_TYPES.index("pi_stacks")] == pytest.approx(np.log1p(1))
    assert np.count_nonzero(profile) == 2


def test_pocket_mapping_preserves_composition_but_randomizes_position_chemistry() -> None:
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 2
    residues = [f"1.A_{index + 1}_{index}" for index in range(20)]
    original = pocket_features(residues, sequence)
    randomized = pocket_features(residues, sequence, randomize_seed=17)
    assert original is not None and randomized is not None
    original_vector, count, fraction = original
    randomized_vector = randomized[0]
    assert count == 20
    assert fraction == 1.0
    assert np.allclose(original_vector[:26], randomized_vector[:26])
    assert not np.allclose(original_vector[34:-1], randomized_vector[34:-1])


def test_dual_cold_purge_removes_every_registered_overlap() -> None:
    columns = {
        "system_pocket_UniProt": ["A", "B", "C", "A"],
        "ligand_rdkit_canonical_smiles": ["x", "y", "z", "q"],
        "scaffold": ["sx", "sy", "sz", "sq"],
        "entry_pdb_id": ["p1", "p2", "p3", "p4"],
        "system_id_no_biounit": ["n1", "n2", "n3", "n4"],
        "cluster": ["c1", "c2", "c3", "c4"],
        "cluster_for_val_split": ["v1", "v2", "v3", "v4"],
        "uniqueness": ["u1", "u2", "u3", "u4"],
        "split": ["train", "train", "train", "val"],
    }
    train, audit, summary = dual_cold_indices(pd.DataFrame(columns))
    assert train.tolist() == [1, 2]
    assert audit.tolist() == [3]
    assert not any(summary["overlap_after"].values())


def test_cross_and_no_cross_residuals_are_parameter_matched() -> None:
    cross = BilinearResidual(42, 83, 8)
    additive_capacity = NoCrossResidual(42, 83, 8)
    assert parameter_count(cross) == parameter_count(additive_capacity)
    ligand = torch.randn(5, 42)
    pocket = torch.randn(5, 83)
    assert cross(ligand, pocket).shape == (5, 8)
    assert torch.equal(cross(ligand, pocket), torch.zeros(5, 8))


def test_additive_model_has_no_pairwise_multiplication() -> None:
    torch.manual_seed(11)
    model = AdditiveContactModel(4, 3, 2)
    ligand = torch.randn(1, 4)
    pocket_a = torch.randn(1, 3)
    pocket_b = torch.randn(1, 3)
    difference_a = model(ligand, pocket_a) - model(torch.zeros_like(ligand), pocket_a)
    difference_b = model(ligand, pocket_b) - model(torch.zeros_like(ligand), pocket_b)
    assert torch.allclose(difference_a, difference_b, atol=1e-6)


def test_matched_donors_are_physically_distinct_and_deterministic() -> None:
    values = np.asarray([10, 11, 20, 21, 30, 31], dtype=float)
    identities = ["a", "b", "c", "d", "e", "f"]
    first = matched_donors(values, identities, seed=91, bins=3)
    second = matched_donors(values, identities, seed=91, bins=3)
    assert np.array_equal(first, second)
    assert np.all(first != np.arange(len(values)))
    assert all(identities[donor] != identities[row] for row, donor in enumerate(first))


def test_normalized_state_pool_is_order_and_duplication_invariant() -> None:
    states = np.asarray([[1.0, 2.0], [3.0, 4.0], [-1.0, 0.0]])
    expected = normalized_state_pool(states)
    assert np.array_equal(expected, normalized_state_pool(states[[2, 0, 1]]))
    assert np.array_equal(
        expected, normalized_state_pool(np.concatenate((states, states[[1]]), axis=0))
    )
