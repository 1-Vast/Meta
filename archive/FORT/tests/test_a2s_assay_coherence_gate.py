"""Contracts for the exact-assay passive information gate."""

from __future__ import annotations

import numpy as np

from research.a2s import a2s_assay_coherence_gate as gate


def test_similarity_bins_cover_the_unit_interval():
    values = np.asarray([0.0, 0.19, 0.20, 0.34, 0.35, 0.54, 0.55, 1.0])
    assert list(gate.similarity_bin(values)) == [
        "t00_20",
        "t00_20",
        "t20_35",
        "t20_35",
        "t35_55",
        "t35_55",
        "t55_100",
        "t55_100",
    ]


def test_k1_permutation_is_a_sign_control():
    episode = gate.Episode("e", "t", "c", "a", "probe", 0, (0, 1, 2, 3, 4), (5,))
    rows, values = gate.arm_rows(episode, episode, 1, "permuted", np.arange(8.0))
    assert list(rows) == [0]
    assert list(values) == [-0.0]


def test_permutation_preserves_support_values_but_changes_assignment():
    episode = gate.Episode("e", "t", "c", "a", "probe", 0, (0, 1, 2, 3, 4), (5,))
    residual = np.asarray([1.0, 2.0, 4.0, 8.0, 16.0, 0.0])
    rows, values = gate.arm_rows(episode, episode, 5, "permuted", residual)
    assert list(rows) == [0, 1, 2, 3, 4]
    assert sorted(values) == [1.0, 2.0, 4.0, 8.0, 16.0]
    assert not np.array_equal(values, residual[rows])


def test_wrong_target_uses_the_donor_rows_and_labels():
    episode = gate.Episode("e", "t", "c", "a", "probe", 0, (0, 1, 2, 3, 4), (10,))
    donor = gate.Episode("d", "u", "z", "b", "probe", 0, (5, 6, 7, 8, 9), (11,))
    residual = np.arange(12.0)
    rows, values = gate.arm_rows(episode, donor, 3, "wrong_target", residual)
    assert list(rows) == [5, 6, 7]
    assert list(values) == [5.0, 6.0, 7.0]


def test_eb_delta_is_zero_when_posterior_equals_the_prior():
    from research.a2s.a2s_mode_generalization import build_subspace

    rng = np.random.default_rng(3)
    heads = rng.normal(size=(100, 4))
    subspace = build_subspace(heads, sigma=0.5)
    design = rng.normal(size=(8, 4))
    support = np.arange(5)
    query = np.arange(5, 8)
    residual = design[support] @ subspace.mean_head
    delta = gate.eb_delta(design, subspace, support, residual, query)
    assert np.allclose(delta, 0.0, atol=1e-8)


def test_decision_requires_both_k3_and_k5_low_similarity_cells():
    def cell(passed: bool) -> dict:
        lower = 0.01 if passed else -0.01
        metrics = {"ci": {"lower95": lower, "mean": lower + 0.01}}
        return {
            "components": 40,
            "contrasts": {
                "eb_desc_minus_base": metrics,
                "eb_desc_correct_minus_permuted": metrics,
                "eb_desc_correct_minus_wrong": metrics,
                "eb_original_minus_base": metrics,
                "eb_original_correct_minus_permuted": metrics,
                "eb_original_correct_minus_wrong": metrics,
            },
        }

    summary = {"probe": {"k3": {"t00_20": cell(True)}, "k5": {"t00_20": cell(True)}}}
    assert "ADMITTED" in gate.decide(summary, {"pass": True})["verdict"]
    summary["probe"]["k3"]["t00_20"] = cell(False)
    assert gate.decide(summary, {"pass": True})["k3_low_similarity_pass"] is False
