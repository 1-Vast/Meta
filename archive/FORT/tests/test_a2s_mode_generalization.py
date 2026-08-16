"""Contracts of the A2S-MODE generalizability gates (G1-G4).

These check the measurement machinery: that within-target splits really are
scaffold-disjoint, that the empirical-Bayes head recovers a known response and
shrinks to the source prior when the evidence is empty, that a rank restriction
actually restricts, and that the verdict logic distinguishes "the object is
local" from "the object is not few-shot reachable".
"""

from __future__ import annotations

import numpy as np
import pytest

from research.a2s import a2s_mode_generalization as generalization


def make_subspace(dimension: int = 6, n_targets: int = 40, seed: int = 0):
    rng = np.random.default_rng(seed)
    heads = rng.normal(0.0, 1.0, (n_targets, dimension))
    return generalization.build_subspace(heads, sigma=0.5), heads


def test_subspace_directions_are_orthonormal_and_ordered():
    subspace, _ = make_subspace()
    gram = subspace.directions.T @ subspace.directions
    assert np.allclose(gram, np.eye(gram.shape[0]), atol=1e-8)
    assert list(subspace.spectrum) == sorted(subspace.spectrum, reverse=True)


def test_subspace_variances_match_the_source_dispersion():
    subspace, heads = make_subspace()
    projected = (heads - heads.mean(axis=0)) @ subspace.directions
    assert np.allclose(subspace.variances, projected.var(axis=0, ddof=1), atol=1e-8)


def test_empirical_bayes_head_recovers_a_response_given_enough_labels():
    subspace, _ = make_subspace(dimension=6, n_targets=200, seed=1)
    rng = np.random.default_rng(2)
    truth = subspace.mean_head + subspace.directions[:, 0] * 1.5
    design = rng.normal(size=(400, 6))
    residual = design @ truth + rng.normal(0.0, 0.05, 400)
    estimate = generalization.empirical_bayes_head(subspace, design, residual, rank=6)
    assert np.corrcoef(estimate, truth)[0, 1] > 0.95


def test_empirical_bayes_head_returns_the_prior_when_evidence_adds_nothing():
    """Residuals exactly explained by the source mean head leave it unchanged."""

    subspace, _ = make_subspace(dimension=6, n_targets=200, seed=3)
    rng = np.random.default_rng(4)
    design = rng.normal(size=(3, 6))
    residual = design @ subspace.mean_head
    estimate = generalization.empirical_bayes_head(subspace, design, residual, rank=6)
    assert np.allclose(estimate, subspace.mean_head, atol=1e-8)


def test_zero_residuals_pull_the_head_toward_no_correction():
    """`r == 0` is evidence that the base is already right, not absence of evidence.

    The posterior must therefore sit between the source prior and the zero head,
    strictly closer to zero than the prior is.
    """

    subspace, _ = make_subspace(dimension=6, n_targets=200, seed=3)
    rng = np.random.default_rng(4)
    design = rng.normal(size=(40, 6))
    estimate = generalization.empirical_bayes_head(subspace, design, np.zeros(40), rank=6)
    prior_prediction = np.linalg.norm(design @ subspace.mean_head)
    posterior_prediction = np.linalg.norm(design @ estimate)
    assert posterior_prediction < prior_prediction


def test_stronger_shrinkage_pulls_the_estimate_back_to_the_prior():
    subspace, _ = make_subspace(dimension=6, n_targets=200, seed=9)
    rng = np.random.default_rng(10)
    design = rng.normal(size=(20, 6))
    residual = design @ (subspace.mean_head + 3.0 * subspace.directions[:, 0]) + rng.normal(0.0, 0.1, 20)
    loose = generalization.empirical_bayes_head(subspace, design, residual, rank=6)
    tight = generalization.build_subspace(
        rng.normal(0.0, 1e-3, (200, 6)) + subspace.mean_head, sigma=subspace.sigma
    )
    tight.mean_head = subspace.mean_head
    shrunk = generalization.empirical_bayes_head(tight, design, residual, rank=6)
    assert np.linalg.norm(shrunk - subspace.mean_head) < np.linalg.norm(loose - subspace.mean_head)


def test_rank_restriction_confines_the_estimate_to_the_subspace():
    subspace, _ = make_subspace(dimension=6, n_targets=200, seed=5)
    rng = np.random.default_rng(6)
    design = rng.normal(size=(80, 6))
    residual = design @ (subspace.mean_head + rng.normal(size=6)) + rng.normal(0.0, 0.1, 80)
    for rank in (1, 2, 3):
        estimate = generalization.empirical_bayes_head(subspace, design, residual, rank=rank)
        offset = estimate - subspace.mean_head
        outside = subspace.directions[:, rank:].T @ offset
        assert np.allclose(outside, 0.0, atol=1e-8)


def test_similarity_bins_partition_the_unit_interval():
    values = np.asarray([0.0, 0.29, 0.30, 0.44, 0.45, 0.59, 0.60, 1.0])
    assigned = generalization.similarity_bin(values)
    assert list(assigned) == [
        "s00_30", "s00_30", "s30_45", "s30_45", "s45_60", "s45_60", "s60_100", "s60_100",
    ]


def test_verdict_separates_local_objects_from_unreachable_ones():
    def summary_with(g1_lower: float, g4_lower: float) -> dict:
        cell = {
            "vs_base": {
                "target_head": {"lower95": g1_lower, "mean": g1_lower, "components": 50},
                "eb_k5_rank26": {"lower95": g4_lower, "mean": g4_lower, "components": 50},
                "protein_zero_shot": {"lower95": -0.05, "mean": -0.01, "components": 50},
            }
        }
        for rank in generalization.RANK_SWEEP:
            cell["vs_base"][f"oracle_rank{rank}"] = {"lower95": 0.0, "mean": 0.01, "components": 50}
        return {"scaffold_disjoint": {"all": cell}}

    subspace, _ = make_subspace()
    assert (
        generalization.decide(summary_with(0.03, 0.02), subspace)["verdict"]
        == "GENERALIZABLE_AND_FEW_SHOT_REACHABLE"
    )
    assert (
        generalization.decide(summary_with(0.03, -0.01), subspace)["verdict"]
        == "GENERALIZABLE_BUT_NOT_YET_FEW_SHOT_REACHABLE"
    )
    assert (
        generalization.decide(summary_with(-0.01, -0.01), subspace)["verdict"]
        == "NOT_GENERALIZABLE_OBJECT_IS_LOCAL"
    )


def test_zero_shot_gate_requires_clearing_the_mde():
    subspace, _ = make_subspace()
    summary = {
        "scaffold_disjoint": {
            "all": {
                "vs_base": {
                    "target_head": {"lower95": 0.03, "mean": 0.04, "components": 50},
                    "protein_zero_shot": {"lower95": 0.004, "mean": 0.02, "components": 50},
                }
            }
        }
    }
    assert generalization.decide(summary, subspace)["gates"]["G3"]["pass"] is False


@pytest.mark.parametrize("rank", [1, 2, 3, 5])
def test_low_rank_estimates_are_nested_in_higher_rank_ones(rank):
    """A rank-r estimate must lie inside every rank-r' >= r subspace."""

    subspace, _ = make_subspace(dimension=6, n_targets=120, seed=7)
    rng = np.random.default_rng(8)
    design = rng.normal(size=(50, 6))
    residual = rng.normal(size=50)
    low = generalization.empirical_bayes_head(subspace, design, residual, rank=rank)
    offset = low - subspace.mean_head
    reconstructed = subspace.directions[:, :rank] @ (subspace.directions[:, :rank].T @ offset)
    assert np.allclose(offset, reconstructed, atol=1e-8)
