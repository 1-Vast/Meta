"""Mechanical contracts for the thermodynamic relation census."""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.a2s import a2s_thermodynamic_relation_census as census


def test_single_cut_fragmentation_finds_a_common_series_core():
    methyl = census.fragment_smiles("Cc1ccccc1")
    ethyl = census.fragment_smiles("CCc1ccccc1")
    methyl_by_core = {fragment.core: fragment.substituent for fragment in methyl}
    ethyl_by_core = {fragment.core: fragment.substituent for fragment in ethyl}
    shared = set(methyl_by_core) & set(ethyl_by_core)
    assert shared
    assert any(methyl_by_core[core] != ethyl_by_core[core] for core in shared)


def test_transform_key_and_orientation_are_direction_stable():
    assert census.transform_key("C[*:1]", "CC[*:1]") == census.transform_key(
        "CC[*:1]", "C[*:1]"
    )
    assert census.orient_replacement("z", 3.0, "a", 1.0) == ("a", 1.0, "z", 3.0)


def test_edge_delta_uses_the_canonical_replacement_direction():
    measurements = pd.DataFrame(
        [
            {
                "role": "fit",
                "component": "c",
                "target": "t",
                "assays": "a",
                "conn": "Cc1ccccc1",
                "affinity": 5.0,
                "base": 4.5,
            },
            {
                "role": "fit",
                "component": "c",
                "target": "t",
                "assays": "a",
                "conn": "CCc1ccccc1",
                "affinity": 7.0,
                "base": 5.5,
            },
        ]
    )
    fragments = census.build_fragment_frame(measurements.conn)
    edges = census.build_mmp_edges(measurements, fragments)
    assert len(edges) == 1
    edge = edges.iloc[0]
    values = {
        "Cc1ccccc1": (5.0, 4.5),
        "CCc1ccccc1": (7.0, 5.5),
    }
    affinity_a, base_a = values[edge.conn_a]
    affinity_b, base_b = values[edge.conn_b]
    assert np.isclose(edge.delta_affinity, affinity_b - affinity_a)
    assert np.isclose(edge.delta_base, base_b - base_a)
    assert np.isclose(edge.delta_residual, (affinity_b - base_b) - (affinity_a - base_a))


def test_edge_builder_allows_a_molecule_measured_on_multiple_targets():
    rows = []
    for target, component in (("t1", "c1"), ("t2", "c2")):
        rows.extend(
            [
                {
                    "role": "fit",
                    "component": component,
                    "target": target,
                    "assays": f"a-{target}",
                    "conn": "Cc1ccccc1",
                    "affinity": 5.0,
                    "base": 4.5,
                },
                {
                    "role": "fit",
                    "component": component,
                    "target": target,
                    "assays": f"a-{target}",
                    "conn": "CCc1ccccc1",
                    "affinity": 7.0,
                    "base": 5.5,
                },
            ]
        )
    measurements = pd.DataFrame(rows)
    fragments = census.build_fragment_frame(measurements.conn)
    edges = census.build_mmp_edges(measurements, fragments)
    assert len(edges) == 2
    assert edges.target.nunique() == 2


def test_numpy_tanimoto_does_not_overflow_uint8_intersections():
    bits = np.zeros((2, 512), dtype=np.uint8)
    bits[:, :300] = 1
    counts = bits.sum(axis=1, dtype=np.int32)
    similarity = census.nearest_tanimoto_bits(
        bits,
        counts,
        np.asarray([0]),
        np.asarray([1]),
    )
    assert np.allclose(similarity, [1.0])


def test_decision_requires_power_at_k3_and_k5_and_positive_transfer():
    coverage = {
        "k3": {"low_similarity": {"components_with_at_least_8_known_queries": 47}},
        "k5": {"low_similarity": {"components_with_at_least_8_known_queries": 47}},
    }
    transfer = {"strata": {"low_similarity": {"proper_gain": {"lower95": 0.01}}}}
    assert "ADMITTED" in census.decide(coverage, transfer)["verdict"]
    coverage["k3"]["low_similarity"]["components_with_at_least_8_known_queries"] = 46
    assert census.decide(coverage, transfer)["verdict"].endswith("NOT_ADMITTED")
