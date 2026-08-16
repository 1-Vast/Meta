"""Contracts of the Gate T0 measurement machinery.

T0 fits only closed-form ridge heads and trains no gradient model, so the only
way it can be wrong is if a statistic does not measure what its name claims.
These tests pin the pieces that carry the verdict: the document-partitioned
concordance, the level-free support fit, the document-offset oracle, and the
decision logic's refusal to read a measurement-context effect as a transferable
object.

Four of them exist because an external review found defects the first suite did
not cover: end-to-end determinism, provenance overlap between support and query,
artifact-hash integrity, and the same-document estimand whose absence let an
all-pair transfer gain be reported as real transfer.
"""

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch

from research.a2s import a2s_transfer_object_gate as t0
from research.a2s.a2s_information_gate import canonical


def test_within_group_concordance_ignores_a_pure_context_offset():
    """Two documents with different offsets must not create within-document skill."""

    rng = np.random.default_rng(0)
    chemistry = rng.normal(size=40)
    document = np.asarray(["A"] * 20 + ["B"] * 20)
    offset = np.where(document == "A", -2.0, 2.0)
    label = chemistry + offset
    metrics = t0.pair_concordance(label, offset, document)
    assert metrics["ci_within"] == pytest.approx(0.5, abs=1e-9)
    assert metrics["ci_between"] > 0.99


def test_between_group_concordance_rewards_the_offset_predictor():
    document = np.asarray(["A"] * 10 + ["B"] * 10)
    label = np.where(document == "A", 0.0, 1.0) + np.arange(20) * 1e-6
    metrics = t0.pair_concordance(label, np.where(document == "A", 0.0, 1.0), document)
    assert metrics["ci_between"] > 0.99
    assert metrics["pairs_within"] == 90


def test_concordance_is_undefined_when_a_partition_is_too_small():
    """Every row in its own document leaves no within-document pair to score."""

    label = np.arange(6, dtype=np.float64)
    group = np.asarray(list("abcdef"))
    metrics = t0.pair_concordance(label, label, group)
    assert metrics["pairs_within"] == 0
    assert np.isnan(metrics["ci_within"])
    assert metrics["ci_between"] == pytest.approx(1.0)


def test_support_fit_loss_is_minimised_by_the_generating_head():
    rng = np.random.default_rng(1)
    heads = rng.normal(size=(12, 6))
    design = rng.normal(size=(30, 6))
    truth = 7
    residual = design @ heads[truth]
    assert int(np.argmin(t0.support_fit_loss(design, residual, heads))) == truth


def test_support_fit_loss_ignores_an_additive_level_shift():
    rng = np.random.default_rng(2)
    heads = rng.normal(size=(9, 5))
    design = rng.normal(size=(20, 5))
    residual = design @ heads[3]
    shifted = t0.support_fit_loss(design, residual + 4.7, heads)
    assert np.allclose(shifted, t0.support_fit_loss(design, residual, heads))


def test_support_fit_loss_is_flat_for_a_single_observation():
    """One label cannot distinguish shapes once the level is removed."""

    rng = np.random.default_rng(3)
    heads = rng.normal(size=(8, 4))
    loss = t0.support_fit_loss(rng.normal(size=(1, 4)), np.asarray([1.3]), heads)
    assert np.allclose(loss, 0.0)


def test_document_offset_oracle_returns_group_means_and_zero_for_unseen():
    prediction = t0.document_offset_prediction(
        np.asarray(["A", "A", "B"]),
        np.asarray([1.0, 3.0, -5.0]),
        np.asarray(["A", "B", "C"]),
    )
    assert prediction == pytest.approx([2.0, -5.0, 0.0])


def test_head_scores_match_an_explicit_loop():
    rng = np.random.default_rng(4)
    design, heads = rng.normal(size=(7, 5)), rng.normal(size=(3, 5))
    scores = t0.head_scores(design, heads)
    for column in range(3):
        assert np.allclose(scores[:, column], design @ heads[column])


def _summary(**overrides) -> dict:
    def band(mean: float, width: float = 0.004) -> dict:
        return {"mean": mean, "lower95": mean - width, "upper95": mean + width, "components": 50}

    cells = {
        "own_minus_base": band(0.052),
        "own_minus_base_within_document": band(0.030),
        "own_minus_base_between_document": band(0.060),
        "docoffset_minus_base": band(0.010),
        "source_best_minus_base": band(0.040),
        "source_median_minus_base": band(0.000),
        "source_pooled_minus_base": band(0.005),
        "select_full_minus_base": band(0.020),
        "select_full_minus_base_within_document": band(0.018),
        "proposal_protein_shortlist_minus_random": band(0.008),
        "proposal_chemotype_shortlist_minus_random": band(0.008),
        "select_k5_minus_random": band(0.012),
        "select_k5_minus_base": band(0.015),
        "select_k5_minus_pooled": band(0.015),
        "select_k3_minus_random": band(0.008),
        "select_k1_minus_random": band(0.000),
    }
    cells.update(overrides)
    return {"scaffold_disjoint": cells}


PASSING_CONTROL = {"pass": True}
BUDGET = {"bits_available": {"k5": 0.1}}


def test_full_pass_reports_the_mechanism_as_admitted():
    decision = t0.decide(_summary(), PASSING_CONTROL, BUDGET)
    assert decision["verdict"] == "DISCRETE_TRANSFER_ADMITTED_PROCEED_TO_MECHANISM"
    assert all(gate["pass"] for gate in decision["gates"].values())


def test_a_context_only_headroom_is_not_reported_as_a_transferable_object():
    """The decisive failure mode: all the gain lives between documents."""

    decision = t0.decide(
        _summary(own_minus_base_within_document={
            "mean": 0.002, "lower95": -0.006, "upper95": 0.010, "components": 50}),
        PASSING_CONTROL,
        BUDGET,
    )
    assert decision["verdict"] == "NO_TRANSFERABLE_CHEMICAL_HEADROOM_OBJECT_IS_MEASUREMENT_CONTEXT"
    assert decision["gates"]["T0A"]["pass"] is False


def test_low_retention_stops_at_non_transferability_even_with_good_selection():
    decision = t0.decide(
        _summary(source_best_minus_base={
            "mean": 0.004, "lower95": -0.001, "upper95": 0.009, "components": 50}),
        PASSING_CONTROL,
        BUDGET,
    )
    assert decision["verdict"] == "ADAPTATION_OBJECT_IS_NOT_TRANSFERABLE_ACROSS_TARGETS"
    assert decision["gates"]["T0B"]["retention"] < t0.T0B_MIN_TRANSFER_RETENTION


def test_transferable_but_unselectable_is_its_own_verdict():
    decision = t0.decide(
        _summary(select_k5_minus_base={
            "mean": -0.010, "lower95": -0.029, "upper95": 0.009, "components": 50}),
        PASSING_CONTROL,
        BUDGET,
    )
    assert decision["verdict"] == "TRANSFERABLE_AT_FULL_SUPPORT_BUT_NOT_IDENTIFIABLE_AT_K5"
    assert decision["gates"]["T0D"]["pass"] is False


def test_beating_random_selection_alone_cannot_admit_the_mechanism():
    """A library whose median member is harmful makes 'better than random' free."""

    decision = t0.decide(
        _summary(
            select_k5_minus_random={
                "mean": 0.016, "lower95": 0.008, "upper95": 0.024, "components": 50},
            select_k5_minus_base={
                "mean": -0.010, "lower95": -0.029, "upper95": 0.009, "components": 50},
            select_k5_minus_pooled={
                "mean": -0.002, "lower95": -0.020, "upper95": 0.015, "components": 50},
        ),
        PASSING_CONTROL,
        BUDGET,
    )
    assert decision["gates"]["T0D"]["pass"] is False
    assert decision["verdict"] == "TRANSFERABLE_AT_FULL_SUPPORT_BUT_NOT_IDENTIFIABLE_AT_K5"


def test_a_failed_positive_control_overrides_every_other_verdict():
    decision = t0.decide(_summary(), {"pass": False}, BUDGET)
    assert decision["verdict"] == "T0_HARNESS_INVALID_NO_POWER"


def test_transfer_is_refused_when_it_does_not_survive_the_same_document_control():
    """The defect that produced a false positive: an all-pair gain, alone.

    T0A measures a per-document offset oracle that outscores the whole chemical
    head, so an all-pair transfer gain is confounded by construction.
    """

    decision = t0.decide(
        _summary(select_full_minus_base_within_document={
            "mean": -0.018, "lower95": -0.044, "upper95": 0.005, "components": 50}),
        PASSING_CONTROL,
        BUDGET,
    )
    assert decision["gates"]["T0B"]["pass"] is False
    assert decision["verdict"] != "DISCRETE_TRANSFER_ADMITTED_PROCEED_TO_MECHANISM"


def test_the_bit_account_ships_its_caveats_with_its_numbers():
    """It was quoted as a parameter-free closure; it is neither."""

    artifact = Path(t0.DEFAULT_OUTPUT)
    if not artifact.exists():
        pytest.skip("gate has not been run in this workspace")
    budget = json.loads(artifact.read_text(encoding="utf-8"))["information_budget"]
    assert budget["status"] == "HEURISTIC_ORDER_OF_MAGNITUDE_ONLY"
    assert len(budget["caveats"]) >= 4
    assert "bits_required_note" in budget


def test_provenance_audit_detects_document_reuse_between_support_and_query():
    frame = pd.DataFrame({
        "docs": ["D1", "D1", "D2", "D1", "D3"],
        "assays": ["A1", "A1", "A2", "A1", "A3"],
    })

    class FakeSubstrate:
        labeled = frame

    class FakeSplit:
        split = "scaffold_disjoint"
        train_rows = np.asarray([0, 1, 2])
        eval_rows = np.asarray([3, 4])

    report = t0.provenance_audit(FakeSubstrate(), [FakeSplit()])["scaffold_disjoint"]
    assert report["query_rows_reusing_a_support_document"] == pytest.approx(0.5)
    assert report["targets_sharing_any_document"] == 1
    assert report["targets_with_every_query_document_seen"] == 0


def test_provenance_audit_reports_a_clean_split_as_clean():
    frame = pd.DataFrame({"docs": ["D1", "D1", "D9"], "assays": ["A1", "A1", "A9"]})

    class FakeSubstrate:
        labeled = frame

    class FakeSplit:
        split = "scaffold_disjoint"
        train_rows = np.asarray([0, 1])
        eval_rows = np.asarray([2])

    report = t0.provenance_audit(FakeSubstrate(), [FakeSplit()])["scaffold_disjoint"]
    assert report["query_rows_reusing_a_support_document"] == 0.0
    assert report["targets_sharing_any_document"] == 0


def test_recorded_content_hash_matches_the_file_on_disk():
    """The runner used to hash, append an artifact block, then rewrite."""

    artifact = Path(t0.DEFAULT_OUTPUT)
    if not artifact.exists():
        pytest.skip("gate has not been run in this workspace")
    payload = json.loads(artifact.read_text(encoding="utf-8"))
    recorded = payload.pop("content_sha256")
    assert sha256(canonical(payload).encode()).hexdigest() == recorded


def test_the_ligand_basis_is_deterministic_across_calls():
    """``torch.svd_lowrank`` took no generator, so every gate was irreproducible."""

    pytest.importorskip("torch")
    from research.a2s.a2s_mode_gates import build_basis
    from research.a2s.a2s_trace import load_substrate
    from research.a2s.a2s_trace_stratum import DEFAULT_LOCK, DEFAULT_OOF

    if not Path(DEFAULT_LOCK).exists():
        pytest.skip("source lock is not present in this workspace")
    substrate, _ = load_substrate(Path(DEFAULT_LOCK), Path(DEFAULT_OOF))
    first, stats = build_basis(substrate)
    second, _ = build_basis(substrate)
    assert np.array_equal(first, second)
    assert "deterministic" in stats["decomposition"]
