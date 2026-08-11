from research.meta_fewshot.numeric_availability_audit import identity_pairs, summarize


def test_summarize_counts_distinct_ligands_per_frozen_side():
    cells = [
        {"target": target, "ligand": f"{target}-{index}"}
        for target in ("source", "eval")
        for index in range(8)
    ]
    pairs = {(cell["target"], cell["ligand"]) for cell in cells}
    result = summarize(
        cells,
        {"source": "root-source", "eval": "root-eval"},
        "root-source",
        {"root-eval"},
        pairs,
    )
    assert result["source_k5_targets"] == 1
    assert result["evaluation_k5_targets"] == 1
    assert result["admitted_pair_overlap"] == 16


def test_identity_pairs_filters_endpoint_and_counts_replicates(tmp_path):
    import gzip
    import json

    path = tmp_path / "labels.jsonl.gz"
    rows = [
        {"endpoint": "Ki", "target_id": "p", "ligand_id": "a", "pK": 1.0},
        {"endpoint": "Ki", "target_id": "p", "ligand_id": "a", "pK": 9.0},
        {"endpoint": "Kd", "target_id": "p", "ligand_id": "b", "pK": 2.0},
    ]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")
    pairs, counts = identity_pairs(path, endpoint="Ki")
    assert pairs == {("p", "a")}
    assert counts[("p", "a")] == 2
