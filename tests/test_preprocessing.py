import json

import pytest
import torch

from scripts.build_ligand_bank import ATOM_FEAT_DIM, BOND_FEAT_DIM, featurize_smiles
from scripts.build_protein_bank import build_protein_bank
from scripts.data_contract import (canonical_sequence, canonical_smiles, read_jsonl,
                                   write_jsonl)


def test_canonical_contract_rejects_invalid_biological_inputs():
    assert canonical_sequence("ac d\n") == "ACD"
    assert canonical_smiles(" CCO ") == "CCO"
    with pytest.raises(ValueError, match="invalid protein"):
        canonical_sequence("ACD-*")
    with pytest.raises(ValueError, match="invalid SMILES"):
        canonical_smiles("not-a-smiles")


def test_fixed_schema_ligand_graph_is_complete():
    graph = featurize_smiles("CC(=O)N[C@@H](C)C(=O)O", max_atoms=32)
    assert tuple(graph["X"].shape) == (32, ATOM_FEAT_DIM)
    assert tuple(graph["A"].shape) == (32, 32, BOND_FEAT_DIM)
    assert graph["mask"].sum() > 0


def test_protein_bank_requires_an_explicit_provider(tmp_path):
    rows = [{"target_key": "p1", "sequence": "ACDE"},
            {"target_key": "p2", "sequence": "FGHI"}]
    rows_path = tmp_path / "rows.jsonl"
    write_jsonl(rows_path, rows)

    def provider(sequences):
        return {key: (torch.ones(4), torch.ones(2, 4)) for key in sequences}

    output = tmp_path / "proteins.pt"
    summary = build_protein_bank(rows_path, output, provider, provider_metadata={
        "provider": "test.fixture", "model_id": "fixture", "model_revision": "sha256:test",
        "tokenizer_revision": "sha256:test", "pooling": "fixture_mean",
        "slot_policy": "fixture_two_slots", "dtype": "float32",
    })
    assert summary["targets"] == 2
    assert summary["schema"] == "MetaSieve.ProteinBank.v3"
    assert summary["tensor_hashes"]["p1"]["residues"]
    assert output.with_suffix(".manifest.json").is_file()
    assert set(torch.load(output, weights_only=False)) == {"p1", "p2"}
