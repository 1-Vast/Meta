import hashlib
import json

from scripts.census_source_affinity import (
    audit_closure_fields,
    audit_davis_protected_homology,
    census_rows,
    census_task_manifest,
    readiness,
    row_is_eligible,
)


def _row(compound, **updates):
    row = {
        "standard_relation": "=",
        "endpoint_family": "Ki",
        "standard_units": "nM",
        "standard_value": 10.0,
        "protein_sequence": "ACDEFG",
        "protein_sequence_sha256": "protein-a",
        "canonical_smiles": "CCO",
        "ligand_connectivity_key": compound,
        "assay_chembl_id": "ASSAY-A",
        "assay_context_sha256": "context-a",
        "document_chembl_id": "DOC-A",
    }
    row.update(updates)
    return row


def test_row_contract_rejects_censoring_missing_sequence_and_bad_units():
    assert row_is_eligible(_row("L1"))
    assert not row_is_eligible(_row("L1", standard_relation="<"))
    assert not row_is_eligible(_row("L1", protein_sequence=""))
    assert not row_is_eligible(_row("L1", standard_units="mg/L"))


def test_row_census_counts_exact_compounds_without_endpoint_pooling():
    rows = [_row(f"L{i}") for i in range(50)]
    rows += [_row(f"KD{i}", endpoint_family="Kd") for i in range(19)]
    rows += [_row("L0"), _row("CENSORED", standard_relation="<")]
    result = census_rows(rows)
    assert result["eligible_rows"] == 70
    assert result["task_count"] == 2
    assert result["tasks_at_threshold"] == {"20": 1, "32": 1, "50": 1}
    assert result["max_exact_compounds"] == 50
    assert result["document_fields"]["complete"]


def test_row_census_filters_to_allowed_sequences():
    rows = [_row(f"L{i}") for i in range(20)]
    excluded = [
        _row(f"X{i}", protein_sequence_sha256="protein-b") for i in range(20)
    ]
    result = census_rows(rows + excluded, {"protein-a"})
    assert result["tasks_at_threshold"]["20"] == 1
    assert result["eligible_rows_excluded_by_split_or_homology"] == 20


def test_manifest_census_is_explicitly_not_row_level_recomputed():
    records = []
    for count in (19, 20, 32, 50):
        records.append({
            "split": "source",
            "protein_sequence_sha256": f"P{count}",
            "endpoint_families": ["Ki"],
            "assays": [f"A{count}"],
            "task_id": f"T{count}",
            "unique_compounds": count,
            "documents": [f"D{count}"],
        })
    payload = {
        "comparison_is_source_only": True,
        "endpoint_pooling_allowed": False,
        "summaries": {"protein_assay_context": {"records": records}},
    }
    result = census_task_manifest(payload)
    assert result["evidence_level"] == "MANIFEST_ONLY"
    assert not result["row_level_recomputed"]
    assert result["tasks_at_threshold"] == {"20": 3, "32": 2, "50": 1}
    assert result["manifest_contract"]["required_task_fields_complete"]


def test_closure_audit_does_not_conflate_internal_split_with_davis_exclusion():
    result = audit_closure_fields({
        "identity_threshold": 0.4,
        "n_pairs_aligned": 10,
        "exhaustive_cross_split_pairs_at_or_above_0_40": 0,
        "homology_components_straddling": 0,
        "documents_straddling": 0,
    })
    assert result["internal_40pct_homology_closure_fields_complete"]
    assert result["document_closure_fields_complete"]
    assert not result["davis_protected_40pct_exclusion_documented"]
    state = readiness({"evidence_level": "MANIFEST_ONLY"}, result)
    assert not state["affinity_mechanism_pilot_ready"]
    assert len(state["blockers"]) == 2


def test_davis_homology_audit_uses_label_free_protected_split(tmp_path):
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    sequence_hash = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    governance = tmp_path / "governance.jsonl"
    governance.write_text(json.dumps({
        "target_key": sequence_hash, "split": "recipient", "rows": 1,
    }) + "\n")
    fasta = tmp_path / "protected.fasta"
    fasta.write_text(f">s_{sequence_hash[:24]}\n{sequence}\n")
    rows = [{
        "protein_sequence": sequence,
        "protein_sequence_sha256": sequence_hash,
    }]
    result = audit_davis_protected_homology(
        rows, {sequence_hash}, governance, fasta,
    )
    assert result["protected_davis_sequences"] == 1
    assert result["excluded_source_sequence_count"] == 1
    assert result["recipient_labels_read"] is False
