from research.meta_fewshot.seal_v1_development import fixed_episodes


def test_fixed_episodes_keep_query_labels_out_of_model_input():
    cells = [
        {"cell_id": str(index), "target_id": "t", "protein_group_40": "g",
         "pK": float(index)} for index in range(8)
    ]
    episodes, truth = fixed_episodes(cells, draws=1)
    assert len(episodes) == 5
    assert all("query_pK" not in row for row in episodes)
    assert all(len(row["support_pK"]) == 5 for row in episodes)
    truth_ids = {row["cell_id"] for row in truth}
    query_ids = {cell_id for row in episodes for cell_id in row["query_cell_ids"]}
    assert truth_ids == query_ids


def test_fixed_episode_support_and_query_are_disjoint():
    cells = [
        {"cell_id": str(index), "target_id": "t", "protein_group_40": "g",
         "pK": float(index)} for index in range(10)
    ]
    episodes, _ = fixed_episodes(cells, draws=1)
    for row in episodes:
        assert set(row["support_cell_ids"]).isdisjoint(row["query_cell_ids"])
