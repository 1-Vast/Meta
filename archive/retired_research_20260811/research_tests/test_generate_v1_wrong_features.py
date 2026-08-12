from research.meta_fewshot.generate_v1_wrong_features import matched_donor_map


def test_wrong_donors_are_source_only_and_cross_group():
    source = [
        {"target_id": "a", "protein_group_40": "g1"},
        {"target_id": "b", "protein_group_40": "g2"},
    ]
    validation = [{"target_id": "v", "protein_group_40": "g3"}]
    sequences = {"a": "A" * 10, "b": "C" * 30, "v": "A" * 11}
    result = matched_donor_map(source, validation, sequences)
    assert result == {"v": "a"}


def test_wrong_donor_never_uses_same_group():
    source = [
        {"target_id": "a", "protein_group_40": "same"},
        {"target_id": "b", "protein_group_40": "other"},
    ]
    validation = [{"target_id": "v", "protein_group_40": "same"}]
    sequences = {"a": "A" * 10, "b": "A" * 20, "v": "A" * 10}
    assert matched_donor_map(source, validation, sequences)["v"] == "b"
