"""Contracts of the A2S-MODE pre-implementation gates (A0-A4).

These check the measurement machinery, not the mechanism: that the compact basis
is label-free, that a head recovers a known linear response, that mode evidence
identifies a mode when one exists, that the k=1 rank channel is structurally
near-silent as designed, and that the gate verdicts require what they claim.
"""

from __future__ import annotations

import numpy as np
import pytest

from research.a2s import a2s_mode_gates as gates


def make_modes(n_modes: int = 3, dimension: int = 6, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    modes = rng.normal(0.0, 1.0, (n_modes, dimension))
    return np.vstack([np.zeros(dimension), modes])


def test_fit_head_recovers_a_linear_response():
    rng = np.random.default_rng(1)
    design = rng.normal(size=(400, 6))
    truth = rng.normal(size=6)
    residual = design @ truth + 1.25
    weight, intercept = gates.fit_head(design, residual, ridge=1e-6)
    assert np.allclose(weight, truth, atol=1e-4)
    assert abs(intercept - 1.25) < 1e-4


def test_fit_head_intercept_is_not_penalised():
    rng = np.random.default_rng(2)
    design = rng.normal(size=(200, 4))
    residual = np.full(200, 3.5)
    weight, intercept = gates.fit_head(design, residual, ridge=1e6)
    assert np.allclose(weight, 0.0, atol=1e-6)
    assert abs(intercept - 3.5) < 1e-6


def test_mode_evidence_identifies_the_true_mode_at_k5():
    rng = np.random.default_rng(3)
    modes = make_modes(n_modes=4, dimension=6, seed=3)
    correct = 0
    for trial in range(60):
        support = rng.normal(size=(5, 6))
        true_index = int(rng.integers(1, modes.shape[0]))
        residual = support @ modes[true_index] + 0.4 + rng.normal(0.0, 0.15, 5)
        evidence, _ = gates.mode_log_evidence(support, residual, modes, sigma=0.2, level_sd=1.0)
        correct += int(np.argmax(evidence) == true_index)
    assert correct / 60 > 0.8


def test_null_mode_wins_when_the_base_is_already_right():
    rng = np.random.default_rng(4)
    modes = make_modes(n_modes=4, dimension=6, seed=4)
    support = rng.normal(size=(5, 6))
    evidence, _ = gates.mode_log_evidence(
        support, np.zeros(5), modes, sigma=0.2, level_sd=1.0
    )
    assert int(np.argmax(evidence)) == 0


def test_k1_rank_channel_is_structurally_near_silent():
    """One label plus a shrunk level cannot separate the modes."""

    rng = np.random.default_rng(5)
    modes = make_modes(n_modes=4, dimension=6, seed=5)
    spreads = []
    for _ in range(50):
        support = rng.normal(size=(1, 6))
        residual = np.asarray([rng.normal(0.0, 1.5)])
        single, _ = gates.mode_log_evidence(support, residual, modes, sigma=1.0, level_sd=1.9)
        five_support = rng.normal(size=(5, 6))
        five_residual = five_support @ modes[2] + rng.normal(0.0, 1.0, 5)
        five, _ = gates.mode_log_evidence(five_support, five_residual, modes, sigma=1.0, level_sd=1.9)
        spreads.append((np.ptp(single), np.ptp(five)))
    single_spread = np.mean([value[0] for value in spreads])
    five_spread = np.mean([value[1] for value in spreads])
    assert single_spread < 0.25 * five_spread


def test_level_shrinkage_grows_with_the_support_budget():
    modes = make_modes(n_modes=2, dimension=3, seed=6)
    rng = np.random.default_rng(6)
    levels = []
    for k in (1, 3, 5, 40):
        support = rng.normal(size=(k, 3))
        residual = np.full(k, 2.0)
        _, level = gates.mode_log_evidence(support, residual, modes, sigma=1.0, level_sd=1.0)
        levels.append(level[0])
    assert levels == sorted(levels)
    assert levels[-1] > 0.9 * 2.0


def test_mode_matrix_puts_the_null_mode_first():
    dictionary = gates.Dictionary(
        modes=np.ones((3, 4)),
        global_head=np.zeros(4),
        sigma=1.0,
        level_sd=1.0,
        n_targets=10,
        assignment={},
    )
    matrix = gates.mode_matrix(dictionary)
    assert matrix.shape == (4, 4)
    assert np.allclose(matrix[0], 0.0)


def test_gate_a1_requires_beating_the_global_head_not_the_base():
    summary = {
        "heldout": {
            "k5": {
                "t55_100": {
                    "contrasts": {
                        # Large gain over the base, none over the global head:
                        # this is "a better global model", not a target state.
                        "a1_splithalf_minus_base": {"ci": {"lower95": 0.05, "mean": 0.06, "components": 40}},
                        "a1_splithalf_minus_global": {"ci": {"lower95": -0.001, "mean": 0.0, "components": 40}},
                        "a2_kshot_minus_global": {"ci": {"lower95": -0.001, "mean": 0.0, "components": 40}},
                    },
                    "absolute_ci": {},
                    "kshot_selection_accuracy": 0.9,
                }
            }
        }
    }
    verdict = gates.decide(summary, n_modes=4)
    assert verdict["gates"]["A1"]["pass"] is False
    assert verdict["verdict"] == "MODE_ROUTE_NOT_ADMITTED"


def test_gate_a2_requires_both_accuracy_and_a_gain():
    def summary_with(accuracy: float, lower: float) -> dict:
        return {
            "heldout": {
                "k5": {
                    "t55_100": {
                        "contrasts": {
                            "a1_splithalf_minus_global": {"ci": {"lower95": 0.02, "mean": 0.03, "components": 40}},
                            "a2_kshot_minus_global": {"ci": {"lower95": lower, "mean": lower, "components": 40}},
                        },
                        "absolute_ci": {},
                        "kshot_selection_accuracy": accuracy,
                    }
                }
            }
        }

    assert gates.decide(summary_with(0.9, 0.02), n_modes=4)["gates"]["A2"]["pass"] is True
    # Chance for a 4-mode dictionary plus the null mode is 1/5.
    assert gates.decide(summary_with(0.15, 0.02), n_modes=4)["gates"]["A2"]["pass"] is False
    assert gates.decide(summary_with(0.9, 0.001), n_modes=4)["gates"]["A2"]["pass"] is False


def test_gate_a3_only_counts_the_transport_null_strata():
    summary = {
        "heldout": {
            "k5": {
                "t55_100": {
                    "contrasts": {
                        "a1_splithalf_minus_global": {"ci": {"lower95": 0.02, "mean": 0.03, "components": 40}},
                        "a2_kshot_minus_global": {"ci": {"lower95": 0.02, "mean": 0.03, "components": 40}},
                    },
                    "absolute_ci": {},
                    "kshot_selection_accuracy": 0.9,
                },
                "t00_20": {
                    "contrasts": {
                        "a2_kshot_minus_global": {"ci": {"lower95": -0.01, "mean": 0.0, "components": 40}},
                    },
                    "absolute_ci": {},
                    "kshot_selection_accuracy": 0.9,
                },
            }
        }
    }
    verdict = gates.decide(summary, n_modes=4)
    assert verdict["gates"]["A2"]["pass"] is True
    assert verdict["gates"]["A3"]["pass"] is False


@pytest.mark.parametrize("n_modes", [2, 4, 6])
def test_chance_level_tracks_the_dictionary_size(n_modes):
    summary = {"heldout": {"k5": {"t55_100": {"contrasts": {}, "absolute_ci": {}}}}}
    records = gates.decide(summary, n_modes=n_modes)["gates"]["A2"]["records"]
    assert all(abs(record["chance"] - 1.0 / (n_modes + 1)) < 1e-12 for record in records)
