"""Contracts of the NEA precondition gates D0, N0 and N1.

The load-bearing claim of D0 is that its split really is separated.  The
strongest available check is structural rather than statistical: on a
document-disjoint split the document-mean oracle predicts a constant, so its
measured gain must be exactly zero.  These tests pin that logic, the coverage
accounting, and the decision rules.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from research.a2s import a2s_nea_preconditions as nea
from research.a2s.a2s_transfer_object_gate import document_offset_prediction, pair_concordance


def test_document_oracle_is_powerless_on_a_document_disjoint_split():
    """The validity check D0 relies on, verified directly."""

    train_documents = np.asarray(["A", "A", "B", "B"])
    train_residual = np.asarray([1.0, 2.0, -3.0, -1.0])
    eval_documents = np.asarray(["X", "X", "Y", "Y"])
    prediction = document_offset_prediction(train_documents, train_residual, eval_documents)
    assert np.all(prediction == 0.0)


def test_document_oracle_is_powerful_when_documents_leak():
    train_documents = np.asarray(["A", "A", "B", "B"])
    train_residual = np.asarray([1.0, 1.0, -1.0, -1.0])
    eval_documents = np.asarray(["A", "B", "A", "B"])
    prediction = document_offset_prediction(train_documents, train_residual, eval_documents)
    assert prediction.tolist() == [1.0, -1.0, 1.0, -1.0]


def test_a_pure_context_offset_produces_no_within_document_ranking_skill():
    document = np.asarray(["A"] * 12 + ["B"] * 12)
    offset = np.where(document == "A", -3.0, 3.0)
    rng = np.random.default_rng(0)
    label = offset + rng.normal(size=24)
    metrics = pair_concordance(label, offset, document)
    assert metrics["ci_within"] == pytest.approx(0.5, abs=1e-9)
    assert metrics["ci_between"] > 0.99


def _summary(separated: dict) -> dict:
    def band(mean: float, width: float = 0.004) -> dict:
        return {"mean": mean, "lower95": mean - width, "upper95": mean + width, "components": 92}

    cells = {
        "own_minus_base": band(0.030),
        "own_minus_base_within_document": band(0.030),
        "docoffset_minus_base": {"mean": 0.0, "lower95": 0.0, "upper95": 0.0, "components": 92},
        "descriptives": {"mean_document_overlap": 0.0},
    }
    cells.update(separated)
    return {"separated": cells}


COVER_OK = {"docs": {"k5": {"fraction_with_a_context_of_2plus": 0.86}}}
COVER_BAD = {"docs": {"k5": {"fraction_with_a_context_of_2plus": 0.10}}}


def test_a_clean_split_with_headroom_admits_the_preconditions():
    decision = nea.decide(_summary({}), {}, COVER_OK)
    assert decision["verdict"] == "NEA_PRECONDITIONS_MET_PROCEED_TO_MECHANISM_DESIGN"
    assert decision["gates"]["D0"]["split_is_clean"] is True


def test_no_surviving_headroom_stops_the_programme():
    decision = nea.decide(
        _summary({"own_minus_base_within_document": {
            "mean": 0.012, "lower95": -0.004, "upper95": 0.028, "components": 92}}),
        {},
        COVER_OK,
    )
    assert decision["verdict"] == "NO_CHEMICAL_ADAPTATION_OBJECT_SURVIVES_SEPARATION_STOP_PROGRAMME"
    assert decision["gates"]["D0"]["pass"] is False


def test_a_leaking_document_oracle_invalidates_the_harness_before_any_verdict():
    """If the oracle still works, the split is not separated and nothing is readable."""

    decision = nea.decide(
        _summary({"docoffset_minus_base": {
            "mean": 0.061, "lower95": 0.039, "upper95": 0.082, "components": 92}}),
        {},
        COVER_OK,
    )
    assert decision["verdict"] == "D0_SPLIT_NOT_CLEAN_HARNESS_INVALID"


def test_document_overlap_alone_invalidates_the_split():
    summary = _summary({})
    summary["separated"]["descriptives"] = {"mean_document_overlap": 0.5}
    assert nea.decide(summary, {}, COVER_OK)["verdict"] == "D0_SPLIT_NOT_CLEAN_HARNESS_INVALID"


def test_missing_coverage_is_reported_as_a_deployment_failure_not_a_data_failure():
    decision = nea.decide(_summary({}), {}, COVER_BAD)
    assert decision["verdict"] == "OBJECT_SURVIVES_BUT_NEA_HAS_NO_K5_DEPLOYMENT_PATH"
    assert decision["gates"]["N1"]["pass"] is False


def test_separated_splits_are_document_and_assay_disjoint_by_construction():
    """The split builder's guarantee, checked on a synthetic target."""

    frame = pd.DataFrame({
        "target": ["T"] * 60,
        "component": ["c"] * 60,
        "role": ["fit"] * 60,
        "docs": ["D1"] * 20 + ["D2"] * 20 + ["D3"] * 20,
        "assays": ["A1"] * 20 + ["A2"] * 20 + ["A3"] * 20,
        "scaffold": [f"s{i}" for i in range(60)],
    })

    class FakeSubstrate:
        labeled = frame

    for split in nea.separated_splits(FakeSubstrate()):
        if split.regime != "separated":
            continue
        train, evaluate = split.train_rows, split.eval_rows
        assert not (set(frame.docs[train]) & set(frame.docs[evaluate]))
        assert not (set(frame.assays[train]) & set(frame.assays[evaluate]))
        assert not (set(frame.scaffold[train]) & set(frame.scaffold[evaluate]))
