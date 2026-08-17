from scripts.audit_affinity_corpus import audit, build_split, verify_split


def _task(task_id, target, endpoint="Kd", documents=()):
    return {
        "task_id": task_id, "protein_sequence_sha256": target,
        "endpoint_family": endpoint, "document_ids": list(documents),
        "exact_compound_count": 20, "non_tied_pair_comparisons": 190,
        "eligible_e0_core": True,
    }


def test_audit_unions_homology_and_document_components():
    tasks = [
        _task("1", "A", documents=("D1",)),
        _task("2", "B", documents=("D1",)),
        _task("3", "C", endpoint="Ki"),
    ]
    homology = [
        {"protein_sequence_sha256": "A", "homology_component_id": "H1"},
        {"protein_sequence_sha256": "B", "homology_component_id": "H2"},
        {"protein_sequence_sha256": "C", "homology_component_id": "H2"},
    ]
    result = audit(tasks, homology)
    assert result["document_components"] == 2
    assert result["homology_components"] == 2
    assert result["target_document_closure_components"] == 1
    assert result["endpoint_distribution"] == {"Kd": 2, "Ki": 1}
    assert result["largest_component_task_fraction"] == 1.0
    split = build_split(tasks, homology)
    assert len({row["outer_oof_fold"] for row in split}) == 1
    assert len({row["closure_component_id"] for row in split}) == 1
    verification = verify_split(tasks, homology, split)
    assert verification["homology_components_straddling"] == 0
    assert verification["documents_straddling"] == 0
