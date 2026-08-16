"""Contracts of the A2S-TRACE Q1 stratum-resolved information gate."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from research.a2s import a2s_trace_stratum as stratum
from research.a2s.a2s_information_gate import metric_loss as reference_metric_loss


def make_group(n_rows: int = 60, seed: int = 0) -> pd.DataFrame:
    """A target rich enough that all three policies are feasible at k=5.

    A policy legitimately returns ``None`` on a target that cannot satisfy it;
    that path is covered by ``test_policies_decline_a_target_they_cannot_serve``.
    """

    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "source_row": np.arange(n_rows, dtype=np.int64),
            "scaffold": [f"S{value % 20}" for value in range(n_rows)],
            "conn": [f"C{value}" for value in range(n_rows)],
            "docs": [f"D{value % 12}" for value in range(n_rows)],
            "assays": [f"A{value % 11}" for value in range(n_rows)],
            "affinity": rng.normal(7.0, 1.0, n_rows),
        }
    )


@pytest.mark.parametrize("k", [1, 3, 5])
def test_every_policy_returns_disjoint_support_and_query(k):
    group = make_group()
    rng = np.random.default_rng(k)
    for name, draw in stratum.POLICY_FUNCTIONS.items():
        selected = draw(group, k, rng)
        assert selected is not None, name
        support, query = selected
        assert len(support) == k, name
        assert len(query) >= stratum.MIN_QUERY, name
        assert not set(support) & set(query), name


def test_policies_decline_a_target_they_cannot_serve():
    """Too few rows or too few distinct scaffolds must yield no episode at all."""

    thin = make_group(n_rows=8)
    assert stratum.draw_random_within_target(thin, 5, np.random.default_rng(0)) is None
    narrow = make_group(n_rows=40)
    narrow["scaffold"] = "S0"
    assert stratum.draw_scaffold_disjoint(narrow, 3, np.random.default_rng(0)) is None


def test_scaffold_disjoint_policy_excludes_support_scaffolds():
    group = make_group()
    rng = np.random.default_rng(7)
    support, query = stratum.draw_scaffold_disjoint(group, 3, rng)
    scaffolds = group.set_index("source_row").scaffold.to_dict()
    support_scaffolds = {scaffolds[row] for row in support}
    assert not {scaffolds[row] for row in query} & support_scaffolds


def test_provenance_disjoint_policy_is_deterministic():
    group = make_group()
    first = stratum.draw_provenance_disjoint(group, 3, np.random.default_rng(1))
    second = stratum.draw_provenance_disjoint(group, 3, np.random.default_rng(999))
    assert first == second


def test_random_policy_depends_on_the_declared_seed():
    group = make_group(n_rows=40)
    a = stratum.draw_random_within_target(group, 5, stratum.episode_rng(1729, "p", "t", 5, 0))
    b = stratum.draw_random_within_target(group, 5, stratum.episode_rng(1730, "p", "t", 5, 0))
    assert a != b


def test_metric_loss_matches_the_gate_implementation():
    rng = np.random.default_rng(0)
    for _ in range(50):
        size = int(rng.integers(2, 30))
        label = rng.normal(size=size).round(2)
        prediction = rng.normal(size=size).round(2)
        fast = stratum.metric_loss(label, prediction)
        slow = reference_metric_loss(label, prediction)
        for key, value in slow.items():
            other = fast[key]
            assert (np.isnan(value) and np.isnan(other)) or abs(value - other) < 1e-9


def test_level_estimator_cannot_change_a_ranking():
    rng = np.random.default_rng(3)
    label = rng.normal(size=16)
    base = rng.normal(size=16)
    shifted = base + 0.83
    assert stratum.metric_loss(label, base)["ci"] == stratum.metric_loss(label, shifted)["ci"]


def test_strata_partition_the_unit_interval():
    values = np.linspace(0.0, 1.0, 101)
    assigned = stratum.stratum_of(values)
    assert set(assigned) <= set(stratum.STRATUM_NAMES)
    assert len(assigned) == len(values)


def test_tanimoto_matrix_matches_the_definition():
    rng = np.random.default_rng(5)
    left = (rng.random((4, 32)) < 0.3).astype(np.float64)
    right = (rng.random((3, 32)) < 0.3).astype(np.float64)
    computed = stratum.tanimoto_matrix(left, right)
    for i in range(left.shape[0]):
        for j in range(right.shape[0]):
            intersection = float(np.sum(left[i] * right[j]))
            union = float(left[i].sum() + right[j].sum() - intersection)
            assert abs(computed[i, j] - intersection / max(union, 1.0)) < 1e-12


def test_paired_bootstrap_aggregates_targets_before_components():
    frame = pd.DataFrame(
        {
            "component": ["a", "a", "a", "b"],
            "target": ["t1", "t1", "t2", "t3"],
            "value": [0.0, 2.0, 3.0, 5.0],
        }
    )
    result = stratum.paired_bootstrap(frame, "value", draws=200)
    # component a -> mean(mean(0,2), 3) = 2.0 ; component b -> 5.0
    assert result["components"] == 2
    assert abs(result["mean"] - 3.5) < 1e-9


def test_admission_requires_both_gain_and_assignment():
    summary = {
        "policy": {
            "k5": {
                "strong": {
                    "contrasts": {
                        "krr_minus_base": {"ci": {"lower95": 0.02, "mean": 0.03, "components": 40}},
                        "krr_correct_minus_deranged": {"ci": {"lower95": 0.01}},
                    }
                },
                "shortcut": {
                    "contrasts": {
                        "krr_minus_base": {"ci": {"lower95": 0.02, "mean": 0.03, "components": 40}},
                        "krr_correct_minus_deranged": {"ci": {"lower95": -0.01}},
                    }
                },
            }
        }
    }
    verdict = stratum.admission(summary)
    by_stratum = {record["stratum"]: record for record in verdict["records"]}
    assert by_stratum["strong"]["admitted"] is True
    assert by_stratum["shortcut"]["admitted"] is False
