"""P1 baseline artifact invariant tests (prereg 59a90ef2 + AD1)."""
import hashlib
import json
import math
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ART = STAGE / "artifacts"


@pytest.fixture(scope="module")
def bl():
    return json.loads((ART / "P1_BASELINES.json").read_text(encoding="utf-8"))


def test_baseline_artifact_schema_and_sha(bl):
    assert bl["schema"] == "MetaSieve.StageP.P1Baselines.v1"
    manifest = json.loads((ART / "P1_BASELINES.manifest.json").read_text(encoding="utf-8"))
    sha = hashlib.sha256((ART / "P1_BASELINES.json").read_bytes()).hexdigest()
    assert sha == manifest["sha256"]
    assert set(bl["arms"]) == {"ligand_only", "tanimoto"}


def test_k0_arms_are_global_mean(bl):
    lo = bl["arms"]["ligand_only"]["p_test:k0"]
    tn = bl["arms"]["tanimoto"]["p_test:k0"]
    assert lo is not None and tn is not None
    assert abs(lo["mse"] - tn["mse"]) < 1e-12
    g = bl["global_mean_pki"]
    assert g > 0
    assert lo["rmse"] > 1.0


def test_tanimoto_beats_ligand_only_at_k5(bl):
    lo = bl["arms"]["ligand_only"]["p_test:k5"]
    tn = bl["arms"]["tanimoto"]["p_test:k5"]
    assert tn["mse"] < lo["mse"]
    assert tn["ci"] > lo["ci"]


def test_k40_records_exist_and_k0_has_most(bl):
    n0 = bl["arms"]["tanimoto"]["p_test:k0"]["n_records"]
    n40 = bl["arms"]["tanimoto"]["p_test:k40"]["n_records"]
    assert n40 > 0
    assert n0 > n40


def test_metrics_within_sane_ranges(bl):
    for arm in ("ligand_only", "tanimoto"):
        for k in (0, 1, 2, 3, 5, 10, 20, 40):
            a = bl["arms"][arm].get(f"p_test:k{k}")
            if a is None:
                continue
            assert 0 < a["mse"] < 10
            assert 0 <= a["ci"] <= 1
            assert 0 <= a["train_seen_frac"] <= 1
