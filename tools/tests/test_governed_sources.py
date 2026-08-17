import json

from scripts.governed_sources import collect_governed_sources


def _write_jsonl(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")


def _accepted(identifier, protein, sequence, eligible=True):
    return {
        "source_row_id": identifier,
        "model_eligible": eligible,
        "standard_relation": "=",
        "endpoint_family": "Ki",
        "p_value": 7.0,
        "ligand_connectivity_key": f"ligand-{identifier}",
        "canonical_smiles": "CCO",
        "protein_sequence": sequence,
        "protein_sequence_sha256": protein,
        "task_keys": {"protein_assay_context": f"{protein}-assay"},
    }


def test_collects_all_roots_before_split_and_builds_homology_only(tmp_path):
    root_a, root_b = tmp_path / "a", tmp_path / "b"
    shared = _accepted("shared", "p1", "AAAA")
    _write_jsonl(root_a / "a.jsonl", [shared, _accepted("a", "p2", "AAAT")])
    _write_jsonl(root_b / "b.jsonl", [shared, _accepted("b", "p3", "CCCC"),
                                      _accepted("skip", "p4", "GGGG", False)])

    def identity(left, right):
        return (0.75 if {left, right} == {"AAAA", "AAAT"} else 0.0, 4, 0)

    output = tmp_path / "governed"
    manifest = collect_governed_sources([root_a, root_b], output, identity)
    identities = [json.loads(line) for line in (output / "identity.jsonl").read_text().splitlines()]
    assignments = [json.loads(line) for line in (output / "homology_assignments.jsonl").read_text().splitlines()]
    component = {row["protein_sequence_sha256"]: row["homology_component_id"] for row in assignments}

    assert len(identities) == 3
    assert manifest["repeated_identical_source_ids"] == 1
    assert manifest["skipped_ineligible_rows"] == 1
    assert manifest["rows_rejected_by_split_homology_scaffold"] == 0
    assert manifest["document_closure_applied"] is False
    assert component["p1"] == component["p2"] != component["p3"]
