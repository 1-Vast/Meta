from research.meta_fewshot.fs_corpus_rebuild import Components, audit_cells, freeze_split


def test_fs_audit_counts_distinct_ligands_and_closed_components():
    components = Components(["a", "b", "c"])
    components.union("a", "b")
    cells = [
        {"target": "a", "ligand": f"a{i}", "scaffold": f"sa{i}", "documents": {"d1"}}
        for i in range(8)
    ] + [
        {"target": "b", "ligand": f"b{i}", "scaffold": f"sb{i}", "documents": {"d2"}}
        for i in range(8)
    ] + [
        {"target": "c", "ligand": f"c{i}", "scaffold": f"sc{i}", "documents": {"d3"}}
        for i in range(3)
    ]
    result = audit_cells(cells, components, homology_candidates_count=2, homology_edges=1)
    assert result["targets_usable_at_k"]["5"] == 2
    assert result["components_with_k5_eligible_target"] == 1
    assert result["largest_component_target_share"] == 2 / 3


def test_frozen_split_keeps_only_largest_component_as_source():
    result = freeze_split({
        "component_summary": [
            {"component": "large", "targets": 20}, {"component": "small", "targets": 2}],
        "eligible_k5_by_component": {"large": 12, "small": 3},
    })
    assert result["source_component"] == "large"
    assert result["evaluation_k5_eligible_targets"] == 3
