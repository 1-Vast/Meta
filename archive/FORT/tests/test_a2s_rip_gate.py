"""Contracts of the A2S-RIP Gate R0 measurement machinery.

These check that the posterior is a real posterior, that the oracle selection is
an oracle and the margin rule is not, that the AUC is a proper AUC, and that the
verdict logic reports the ceiling and the implementable rule as separate things.
"""

from __future__ import annotations

import numpy as np
import pytest

from research.a2s import a2s_rip_gate as rip
from research.a2s.a2s_mode_generalization import build_subspace


def make_subspace(dimension: int = 6, n_targets: int = 120, sigma: float = 0.5, seed: int = 0):
    rng = np.random.default_rng(seed)
    return build_subspace(rng.normal(0.0, 1.0, (n_targets, dimension)), sigma=sigma)


def test_posterior_covariance_is_symmetric_positive_semidefinite():
    subspace = make_subspace()
    rng = np.random.default_rng(1)
    design = rng.normal(size=(5, 6))
    posterior = rip.head_posterior(subspace, design, rng.normal(size=5))
    assert np.allclose(posterior.covariance, posterior.covariance.T, atol=1e-8)
    assert np.linalg.eigvalsh(posterior.covariance).min() > -1e-8


def test_posterior_uncertainty_shrinks_as_labels_accumulate():
    subspace = make_subspace()
    rng = np.random.default_rng(2)
    traces = []
    for k in (2, 5, 20, 100):
        design = rng.normal(size=(k, 6))
        residual = design @ subspace.mean_head + rng.normal(0.0, 0.5, k)
        traces.append(np.trace(rip.head_posterior(subspace, design, residual).covariance))
    assert traces == sorted(traces, reverse=True)


def test_posterior_mean_returns_the_prior_when_evidence_adds_nothing():
    subspace = make_subspace()
    rng = np.random.default_rng(3)
    design = rng.normal(size=(4, 6))
    residual = design @ subspace.mean_head
    posterior = rip.head_posterior(subspace, design, residual)
    assert np.allclose(posterior.mean, subspace.mean_head, atol=1e-8)


def test_marginal_gain_is_positive_exactly_for_helpful_edits():
    base = np.asarray([0.0, 1.0, 2.0, 3.0])
    labels = np.asarray([3.0, 2.0, 1.0, 0.0])  # base ordering is exactly inverted
    delta = np.asarray([10.0, 0.0, 0.0, -10.0])  # move the extremes the right way
    gain = rip.marginal_concordance_gain(base, labels, delta)
    assert gain[0] > 0 and gain[3] > 0
    assert gain[1] == 0 and gain[2] == 0


def test_marginal_gain_is_negative_for_harmful_edits():
    base = np.asarray([0.0, 1.0, 2.0])
    labels = np.asarray([0.0, 1.0, 2.0])  # base ordering is already correct
    gain = rip.marginal_concordance_gain(base, labels, np.asarray([10.0, 0.0, 0.0]))
    assert gain[0] < 0


@pytest.mark.parametrize("coverage", [0.0, 0.25, 0.5, 1.0])
def test_selection_respects_the_coverage_budget(coverage):
    rng = np.random.default_rng(4)
    statistics = {name: rng.normal(size=20) for name in rip.RULES}
    for rule in rip.RULES:
        chosen = rip.select(rule, coverage, statistics, rng)
        assert int(chosen.sum()) == int(round(coverage * 20))


def test_selection_takes_the_highest_scoring_compounds():
    statistics = {"margin": np.asarray([0.1, 0.9, 0.5, 0.7]), "random": np.zeros(4)}
    chosen = rip.select("margin", 0.5, statistics, np.random.default_rng(0))
    assert list(np.flatnonzero(chosen)) == [1, 3]


def test_rank_auc_matches_a_brute_force_count():
    rng = np.random.default_rng(5)
    scores = rng.normal(size=60)
    positive = rng.random(60) < 0.4
    pairs = [
        (1.0 if a > b else 0.5 if a == b else 0.0)
        for a in scores[positive]
        for b in scores[~positive]
    ]
    assert abs(rip.rank_auc(scores, positive) - float(np.mean(pairs))) < 1e-9


def test_rank_auc_is_undefined_without_both_classes():
    scores = np.arange(10, dtype=np.float64)
    assert np.isnan(rip.rank_auc(scores, np.ones(10, dtype=bool)))
    assert np.isnan(rip.rank_auc(scores, np.zeros(10, dtype=bool)))


def test_verdict_separates_the_ceiling_from_the_implementable_rule():
    """A real ceiling with a useless margin must not read as a working mechanism."""

    def summary_with(margin_gain: float) -> dict:
        points = {
            f"c{coverage:.1f}": {
                "vs_base": {"mean": 0.05, "lower95": 0.04, "upper95": 0.06, "components": 50},
                "vs_wholesale": {"mean": 0.06, "lower95": 0.05, "upper95": 0.07, "components": 50},
                "vs_magnitude_matched": {"mean": 0.05, "lower95": 0.04, "upper95": 0.06, "components": 50},
                "mean_harm_rate": 0.3,
                "mean_abs_delta": 0.5,
            }
            for coverage in (0.4,)
        }
        margin_points = {
            "c0.4": {
                "vs_base": {"mean": margin_gain, "lower95": margin_gain - 0.01, "upper95": margin_gain + 0.01, "components": 50},
                "vs_wholesale": {"mean": margin_gain, "lower95": margin_gain - 0.01, "upper95": margin_gain + 0.01, "components": 50},
                "vs_magnitude_matched": {"mean": margin_gain, "lower95": margin_gain - 0.01, "upper95": margin_gain + 0.01, "components": 50},
                "mean_harm_rate": 0.4,
                "mean_abs_delta": 0.2,
            }
        }
        auc = {"auc_margin": {"mean": 0.55, "lower95": 0.52, "upper95": 0.58, "components": 50}}
        return {
            "probe": {
                "k5": {"oracle": points, "margin": margin_points, "auc": auc},
                "k3": {"oracle": points, "margin": margin_points, "auc": auc},
            }
        }

    useless = rip.decide(summary_with(0.0), {})
    assert useless["gates"]["R0a"]["pass"] is True
    assert useless["gates"]["R0d"]["pass"] is False
    assert useless["gates"]["R0d"]["oracle_pass"] is True

    useful = rip.decide(summary_with(0.05), {})
    assert useful["gates"]["R0d"]["pass"] is True


def test_no_ceiling_is_reported_when_the_oracle_does_not_beat_wholesale():
    summary = {
        "probe": {
            "k5": {
                "oracle": {
                    "c0.4": {
                        "vs_base": {"mean": 0.0, "lower95": -0.01, "upper95": 0.01, "components": 50},
                        "vs_wholesale": {"mean": 0.0, "lower95": -0.01, "upper95": 0.01, "components": 50},
                        "vs_magnitude_matched": {"mean": 0.0, "lower95": -0.01, "upper95": 0.01, "components": 50},
                        "mean_harm_rate": 0.4,
                        "mean_abs_delta": 0.5,
                    }
                },
                "auc": {},
            }
        }
    }
    assert rip.decide(summary, {})["verdict"] == "NO_SELECTION_CEILING"
