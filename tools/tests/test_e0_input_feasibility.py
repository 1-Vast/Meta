import json

from scripts.audit_e0_input_feasibility import audit_e0_inputs


def _write_jsonl(path, rows):
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")


def test_audit_drops_affinity_values_and_counts_contract_failures(tmp_path):
    task = "t"
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    import hashlib
    protein = hashlib.sha256(sequence.encode()).hexdigest()
    splits = tmp_path / "splits.jsonl"
    rows = tmp_path / "rows.jsonl"
    _write_jsonl(splits, [{
        "task_id": task, "outer_oof_fold": 0, "closure_component_id": "c",
        "protein_sequence_sha256": protein,
    }])
    _write_jsonl(rows, [{
        "task_id": task, "ligand_connectivity_key": "ligand", "canonical_smiles": "CCO",
        "standard_inchi_key": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "protein_sequence": sequence, "protein_sequence_sha256": protein,
        "p_affinity": "DO_NOT_READ", "standard_value": "DO_NOT_READ",
    }])
    result = audit_e0_inputs(rows, splits)
    assert result["affinity_value_fields_materialized"] is False
    assert result["counts"]["governed_rows"] == 1
    assert result["counts"]["valid_rows"] == 1
    assert result["counts"]["unique_valid_ligands"] == 1
    assert result["post_contract_floor"]["pass"] is False
