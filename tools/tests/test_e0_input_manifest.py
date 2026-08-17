import hashlib
import json

from scripts.build_e0_input_manifest import build_e0_input_manifest


def _write_jsonl(path, rows):
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")


def test_manifest_is_label_free_and_uses_stereo_state_key(tmp_path):
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    protein = hashlib.sha256(sequence.encode()).hexdigest()
    splits = tmp_path / "splits.jsonl"
    rows = tmp_path / "rows.jsonl"
    output = tmp_path / "output"
    _write_jsonl(splits, [{"task_id": "t", "outer_oof_fold": 2,
                           "closure_component_id": "c",
                           "protein_sequence_sha256": protein}])
    values = []
    for index in range(20):
        values.append({
            "activity_id": index, "task_id": "t", "endpoint_family": "Ki",
            "ligand_connectivity_key": f"connectivity-{index}",
            "standard_inchi_key": f"state-{index}", "canonical_smiles": "CCO",
            "protein_sequence": sequence, "protein_sequence_sha256": protein,
            "p_affinity": "DO_NOT_READ",
        })
    _write_jsonl(rows, values)
    result = build_e0_input_manifest(rows, splits, output)
    assert result["affinity_values_present"] is False
    assert result["counts"] == {"rows": 20, "tasks": 1, "ligand_states": 20,
                                "proteins": 1}
    text = (output / "rows.label_blind.jsonl").read_text(encoding="utf-8")
    assert "p_affinity" not in text
    assert "state-0" in text
