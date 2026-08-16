"""Invariants for the A2S-BIR mechanism and the Gate D0-R roster.

These are protocol tests, not performance tests. They assert the leakage and
identifiability contracts that the scientific claim depends on.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch

from research.a2s import a2s_bir as bir
from research.a2s import a2s_d0r as d0r

ROSTER = Path(__file__).resolve().parents[1] / "dataset" / "formal_training" / "a2s_d0r_roster.v2"
pytestmark = pytest.mark.skipif(not ROSTER.exists(), reason="Gate D0-R roster not built")


# ------------------------------------------------------------------ roster
def _report() -> dict:
    return json.loads((ROSTER / "d0r_report.json").read_text(encoding="utf-8"))


def test_roster_meets_frozen_d0r_rules() -> None:
    r = _report()
    assert r["status"] == "PASS"
    assert r["counts"]["recipients"] >= d0r.RECIPIENT_FLOOR
    assert r["counts"]["independent_components"] >= d0r.COMPONENT_FLOOR
    assert r["counts"]["largest_component_share"] <= d0r.MAX_COMPONENT_SHARE
    assert r["query_depth"]["min"] >= d0r.QUERY_FLOOR


def test_roster_hard_overlaps_are_zero() -> None:
    overlap = _report()["overlap"]
    for axis in ("source_recipient_target_uid", "source_recipient_accession",
                 "source_recipient_document_uid", "source_recipient_parent",
                 "source_recipient_assay_id", "support_query_parent"):
        assert overlap[axis] == 0, f"{axis} leaked: {overlap[axis]}"


def test_roster_selection_is_label_blind() -> None:
    assert _report()["labels_used_for_roster_selection"] is False


def test_support_draws_are_nested_and_distinct() -> None:
    draws = pd.read_parquet(ROSTER / "support_draws.parquet")
    for (target, draw_id), frame in draws.groupby(["target_uid", "draw_id"]):
        by_k = {k: set(g.compound_parent_uid) for k, g in frame.groupby("k")}
        for small, large in ((1, 3), (3, 5)):
            if small in by_k and large in by_k:
                assert by_k[small] <= by_k[large], f"{target} draw {draw_id} not nested"
    # the five draws per recipient must be distinct at k=5
    top = draws[draws.k == draws.k.max()]
    for target, frame in top.groupby("target_uid"):
        keys = {tuple(sorted(g.compound_parent_uid)) for _, g in frame.groupby("draw_id")}
        assert len(keys) == frame.draw_id.nunique(), f"{target} repeats a support set"


def test_support_precedes_query_in_document_time() -> None:
    """The estimand requires strict temporal precedence, per recipient."""
    corpus_rows = pd.read_parquet(
        d0r.CORPUS / "canonical" / "pki_measurements_exact.parquet",
        columns=["target_uid", "compound_parent_uid", "document_year"])
    latest_year = corpus_rows.groupby(["target_uid", "compound_parent_uid"]).document_year.min()
    recipients = pd.read_parquet(ROSTER / "recipients.parquet").set_index("target_uid")
    draws = pd.read_parquet(ROSTER / "support_draws.parquet")
    query = pd.read_parquet(ROSTER / "query.parquet")
    for target, tau in recipients.tau.items():
        sup = draws[draws.target_uid == target].compound_parent_uid.unique()
        qry = query[query.target_uid == target].compound_parent_uid.unique()
        assert all(latest_year.loc[(target, p)] <= tau for p in sup)
        assert all(latest_year.loc[(target, p)] > tau for p in qry)


# ------------------------------------------------------------------ mechanism
def _hyper() -> bir.Hyper:
    return bir.Hyper(tau_b=0.28, tau_z=0.18, sigma=1.0, drift=0.0)


def test_anchor_is_shrunk_toward_the_global_mean() -> None:
    """tau_b << sigma must produce James-Stein shrinkage, not the raw mean."""
    h = _hyper()
    r = torch.tensor([2.0, 2.0, 2.0], device=bir.DEVICE)
    g = torch.zeros((3, 0), device=bir.DEVICE)
    mean, _ = bir.hierarchical_posterior(g, r, h)
    assert 0.0 < float(mean[0]) < 2.0
    # more support must shrink less
    r5 = torch.full((5,), 2.0, device=bir.DEVICE)
    mean5, _ = bir.hierarchical_posterior(torch.zeros((5, 0), device=bir.DEVICE), r5, h)
    assert float(mean5[0]) > float(mean[0])


def test_certificate_rejects_an_unidentifiable_code() -> None:
    """Collinear support cannot identify a code direction."""
    h = _hyper()
    g = torch.ones((3, 1), device=bir.DEVICE) * 0.01      # near-zero design
    r = torch.tensor([0.1, -0.1, 0.05], device=bir.DEVICE)
    ok, contraction = bir.certificate(g, r, h)
    assert not ok and contraction < bir.CERT_CONTRACTION


def test_certificate_accepts_a_well_conditioned_code() -> None:
    h = bir.Hyper(tau_b=0.28, tau_z=1.0, sigma=0.1, drift=0.0)
    g = torch.tensor([[-2.0], [0.0], [2.0]], device=bir.DEVICE)
    r = torch.tensor([-1.0, 0.0, 1.0], device=bir.DEVICE)
    ok, contraction = bir.certificate(g, r, h)
    assert ok and contraction >= bir.CERT_CONTRACTION


def test_local_residual_effective_dof_is_bounded_by_k() -> None:
    """The budget constraint: realised degrees of freedom never exceed k."""
    torch.manual_seed(0)
    model = bir.LocalResidual(64, rank=8).to(bir.DEVICE)
    for k in (1, 3, 5):
        g = torch.randn(k, 64, device=bir.DEVICE)
        assert 0.0 <= model.effective_dof(g) <= k + 1e-6


def test_local_residual_falls_back_to_anchor_without_coverage() -> None:
    """A query orthogonal to every support compound must receive the anchor."""
    torch.manual_seed(0)
    model = bir.LocalResidual(4, rank=2).to(bir.DEVICE)
    with torch.no_grad():
        model.gate_slope.fill_(50.0)
        model.gate_bias.fill_(0.9)
    g_s = torch.tensor([[1.0, 0, 0, 0], [0, 1.0, 0, 0]], device=bir.DEVICE)
    g_q = torch.tensor([[0, 0, 0, 1.0]], device=bir.DEVICE)
    anchor = torch.tensor(0.7, device=bir.DEVICE)
    r_s = torch.tensor([3.0, -3.0], device=bir.DEVICE)
    out = model(g_q, g_s, r_s, anchor)
    assert float(out[0]) == pytest.approx(float(anchor), abs=1e-3)


def test_code_budget_is_monotone_and_zero_at_k_equals_one() -> None:
    assert bir.CODE_BUDGET[1] == 0, "one label identifies only the anchor"
    assert bir.CODE_BUDGET[1] <= bir.CODE_BUDGET[3] <= bir.CODE_BUDGET[5]
    for k, m in bir.CODE_BUDGET.items():
        assert m < k, "code dimension must leave a degree of freedom for the anchor"


def test_episode_protocols_respect_their_contracts() -> None:
    rng = np.random.default_rng(0)
    frame = pd.DataFrame({
        "compound_parent_uid": [f"c{i}" for i in range(20)],
        "document_uid": [f"d{i // 5}" for i in range(20)],
        "document_year": [2000 + i // 5 for i in range(20)],
    })
    sup, qry = bir.document_ordered_episode(frame, 5, rng)
    years = dict(zip(frame.compound_parent_uid, frame.document_year))
    assert max(years[p] for p in sup) < min(years[p] for p in qry)
    assert not set(sup) & set(qry)
    sup_r, qry_r = bir.random_episode(frame, 5, rng)
    assert not set(sup_r) & set(qry_r)


def test_constant_predictor_scores_at_chance_not_zero() -> None:
    y = np.array([5.0, 6.0, 7.0, 8.0])
    const = np.full(4, 6.5)
    assert bir.spearman(y, const) == 0.0
    assert bir.pairwise_accuracy(y, const) == pytest.approx(0.5)
