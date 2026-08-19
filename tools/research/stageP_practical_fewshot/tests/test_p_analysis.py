"""P1 strata / P2 screening artifact tests (prereg 59a90ef2 + AD1)."""
import hashlib
import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ART = HERE.parent / "artifacts"


def _sha(name):
    return hashlib.sha256((ART / name).read_bytes()).hexdigest()


def test_strata_artifact_pinned():
    manifest = json.loads((ART / "P1_STRATA.json.manifest.json").read_text(encoding="utf-8"))
    assert _sha("P1_STRATA.json") == manifest["sha256"]
    data = json.loads((ART / "P1_STRATA.json").read_text(encoding="utf-8"))
    assert data["schema"] == "MetaSieve.StageP.P1Strata.v1"
    assert set(data["strata"]) == {"ligand_only", "tanimoto"}


def test_screening_artifact_pinned():
    manifest = json.loads((ART / "P2_SCREENING.json.manifest.json").read_text(encoding="utf-8"))
    assert _sha("P2_SCREENING.json") == manifest["sha256"]
    data = json.loads((ART / "P2_SCREENING.json").read_text(encoding="utf-8"))
    assert data["schema"] == "MetaSieve.StageP.P2Screening.v1"
    assert data["active_threshold_pki"] == 6.0
    k0 = data["screening"]["tanimoto"]["p_test:k0"]
    assert k0["degenerate_ranker"] is True
    k5 = data["screening"]["tanimoto"]["p_test:k5"]
    assert 0 <= k5["pr_auc"] <= 1
    assert k5["ef_5pct"] > 1.0


def test_high_similarity_band_easier_than_low():
    data = json.loads((ART / "P1_STRATA.json").read_text(encoding="utf-8"))
    for k in ("k5", "k10", "k20"):
        agg = data["strata"]["tanimoto"][f"p_test:{k}"]
        assert agg["high"]["mse"] < agg["low"]["mse"]
