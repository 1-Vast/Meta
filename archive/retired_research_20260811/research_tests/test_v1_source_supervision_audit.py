from research.meta_fewshot.v1_source_supervision_audit import (
    admitted_source_rows,
    median_rows,
    summarize,
)


def test_admitted_source_rows_excludes_meta_test():
    cells = [
        {"split": "meta_train", "target_id": "t1", "ligand_id": "l1",
         "source_row_ids": ["a"]},
        {"split": "meta_val", "target_id": "t2", "ligand_id": "l2",
         "source_row_ids": ["b"]},
        {"split": "meta_test", "target_id": "t3", "ligand_id": "l3",
         "source_row_ids": ["c"]},
    ]
    allowed, forbidden = admitted_source_rows(cells)
    assert allowed == {("a", "t1", "l1"): "meta_train",
                       ("b", "t2", "l2"): "meta_val"}
    assert forbidden == {("c", "t3", "l3")}


def test_median_rows_aggregates_only_exact_group():
    rows = [
        {"split": "meta_train", "panel_id": "p", "target_id": "t",
         "ligand_id": "l", "pK": 6.0},
        {"split": "meta_train", "panel_id": "p", "target_id": "t",
         "ligand_id": "l", "pK": 8.0},
    ]
    result = median_rows(rows)
    assert result[0]["pK"] == 7.0
    assert result[0]["replicates"] == 2


def test_partner_group_requires_different_frozen_families():
    rows = [
        {"split": "meta_train", "panel_id": "p", "target_id": "a",
         "ligand_id": "x", "pK": 6.0, "replicates": 1},
        {"split": "meta_train", "panel_id": "p", "target_id": "b",
         "ligand_id": "x", "pK": 8.0, "replicates": 1},
        {"split": "meta_train", "panel_id": "p", "target_id": "a",
         "ligand_id": "y", "pK": 7.0, "replicates": 1},
    ]
    result = summarize(rows, {"a": "g1", "b": "g2"}, "meta_train")
    assert result["measured_partner_groups_cross_cdhit40"] == 1
    assert result["measured_cross_family_pairs"] == 1
    assert result["absolute_delta_pK"]["median"] == 2.0
