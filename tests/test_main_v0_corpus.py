import json

from research.meta_fewshot.build_main_v0_corpus import (
    aggregate_observations,
    assign_clusters,
    parse_cdhit_clusters,
)


def test_two_stage_median_weights_panels_equally():
    metadata = {
        str(index): {"chain_count": 1, "target_sequence": "ACDE",
                     "ligand_smiles": "CCO"}
        for index in range(4)
    }
    labels = [
        {"endpoint": "Ki", "source_row_id": "0", "panel_id": "a",
         "target_id": "p", "ligand_id": "l", "pK": 1.0},
        {"endpoint": "Ki", "source_row_id": "1", "panel_id": "a",
         "target_id": "p", "ligand_id": "l", "pK": 3.0},
        {"endpoint": "Ki", "source_row_id": "2", "panel_id": "a",
         "target_id": "p", "ligand_id": "l", "pK": 100.0},
        {"endpoint": "Ki", "source_row_id": "3", "panel_id": "b",
         "target_id": "p", "ligand_id": "l", "pK": 5.0},
    ]
    cells, counts = aggregate_observations(metadata, labels)
    assert len(cells) == 1
    assert cells[0]["pK"] == 4.0  # median([median(1,3,100), median(5)])
    assert cells[0]["panel_count"] == 2
    assert counts["cross_panel_pairs"] == 1


def test_cluster_assignment_is_complete_deterministic_and_cluster_closed():
    mapping = {f"p{i}": f"c{i // 3}" for i in range(30)}
    first = assign_clusters(mapping, seed=7)
    second = assign_clusters(mapping, seed=7)
    assert first == second
    assert set(first) == set(mapping)
    for left, left_cluster in mapping.items():
        for right, right_cluster in mapping.items():
            if left_cluster == right_cluster:
                assert first[left] == first[right]


def test_parse_cdhit_clusters(tmp_path):
    left, right = "a" * 64, "b" * 64
    path = tmp_path / "x.clstr"
    path.write_text(
        f">Cluster 0\n0 10aa, >{left}... *\n"
        f">Cluster 1\n0 10aa, >{right}... *\n",
        encoding="utf-8",
    )
    parsed = parse_cdhit_clusters(path)
    assert set(parsed) == {left, right}
    assert parsed[left] != parsed[right]
