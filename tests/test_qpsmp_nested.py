from types import SimpleNamespace

import numpy as np

from scripts.evaluate_qpsmp import (
    build_nested_manifest, episode_spec, manifest_payload,
    paired_component_bootstrap, paired_component_effects,
)


def synthetic_data():
    cells, tasks = [], {"meta_test": {}}
    components = {"meta_test": {"c1": ("t1",), "c2": ("t2",)}}
    for component, target in (("c1", "t1"), ("c2", "t2")):
        indices = []
        for index in range(10):
            indices.append(len(cells))
            cells.append({"cell_id": f"{target}-{index}", "target_id": target,
                          "protein_group_40": component, "split": "meta_test"})
        tasks["meta_test"][target] = np.asarray(indices)
    return SimpleNamespace(cells=cells, tasks=tasks, components=components)


def test_nested_manifest_uses_common_query_and_support_prefixes():
    data = synthetic_data()
    records = build_nested_manifest(data, "meta_test", (1, 2, 3, 5), 3, 2, 17)

    assert len(records) == 4
    assert {record.target for record in records} == {"t1", "t2"}
    for record in records:
        specs = [episode_spec(data, record, k) for k in (1, 2, 3, 5)]
        assert all(spec.query == specs[0].query for spec in specs)
        assert specs[0].support == specs[-1].support[:1]
        assert specs[1].support == specs[-1].support[:2]
        assert set(specs[-1].support).isdisjoint(specs[-1].query)


def test_manifest_payload_contains_cell_ids_not_affinities():
    data = synthetic_data()
    records = build_nested_manifest(data, "meta_test", (1, 5), 2, 1, 19)
    payload = manifest_payload(records, seed=19, support_sizes=(1, 5), query_size=2)

    assert payload["selection_uses_affinity"] is False
    assert "pK" not in str(payload)


def test_all_targets_with_one_query_after_max_k_are_included():
    data = synthetic_data()
    data.tasks["meta_test"]["t1"] = data.tasks["meta_test"]["t1"][:6]

    records = build_nested_manifest(data, "meta_test", (1, 5), 8, 1, 21)

    target = next(record for record in records if record.target == "t1")
    assert len(target.support_cell_ids) == 5
    assert len(target.query_cell_ids) == 1


def test_paired_component_bootstrap_averages_within_component_first():
    rows = []
    for seed in (1, 2, 3):
        for component, target, effect in (("a", "a1", 2.0), ("b", "b1", -1.0)):
            rows.extend([
                {"model_seed": seed, "component": component, "target": target,
                 "draw": 0, "k": 5, "arm": "full", "mse_pk": 3.0},
                {"model_seed": seed, "component": component, "target": target,
                 "draw": 0, "k": 5, "arm": "control", "mse_pk": 3.0 + effect},
            ])

    effects = paired_component_effects(rows, "full", "control", 5)
    result = paired_component_bootstrap(rows, "full", "control", 5, 99, 23)

    assert effects == {"a": 2.0, "b": -1.0}
    assert np.isclose(result["mean_mse_reduction_pk"], 0.5)
    assert result["components"] == 2
