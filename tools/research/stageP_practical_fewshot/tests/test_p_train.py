"""P1 arm-3 (ordinary FT) invariant tests (CPU; dry artifacts allowed)."""
import hashlib
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))
ART = STAGE / "artifacts"

import torch  # noqa: E402
import p_train as T  # noqa: E402


def test_backbone_shapes_and_grad():
    m = T.PTrunk()
    xp = torch.randn(7, 640)
    xl = torch.randn(7, 2048)
    out = m(xp, xl)
    assert out["yhat"].shape == (7,)
    assert out["inter"].shape == (7,)
    loss = (out["yhat"] ** 2).mean()
    loss.backward()
    names = [n for n, p in m.named_parameters() if p.grad is not None and p.grad.abs().sum() > 0]
    assert "A.weight" in names and "B.weight" in names and "mu" in names
    assert m.A.weight.shape == (16, 64) and m.B.weight.shape == (16, 64)


def test_k0_path_has_no_adaptation():
    art = json.loads((ART / "P1_ARM3_ORDINARYFT.json").read_text(encoding="utf-8"))
    seed0 = art["seeds"]["1"]
    for rec in seed0["records"]:
        if rec["k"] == 0:
            assert rec["best_support_loss"] is None


def test_artifact_schema_and_manifest_sha():
    art = json.loads((ART / "P1_ARM3_ORDINARYFT.json").read_text(encoding="utf-8"))
    assert art["schema"] == "MetaSieve.StageP.P1Arm3OrdinaryFT.v1"
    assert art["backbone_spec_sha256"] == T.BACKBONE_SHA
    manifest = json.loads(
        (ART / "P1_ARM3_ORDINARYFT.json.manifest.json").read_text(encoding="utf-8"))
    sha = hashlib.sha256((ART / "P1_ARM3_ORDINARYFT.json").read_bytes()).hexdigest()
    assert sha == manifest["sha256"]
    recs = art["seeds"]["1"]["records"]
    assert any(r["split"] == "p_test" and r["k"] == 40 for r in recs)
    for r in recs:
        assert 0 < r["mse"] < 1e6
        assert 0 <= r["ci"] <= 1


def test_protocol_values_frozen():
    assert T.STEPS == 6000 and T.ADAPT_STEPS == 50 and T.ADAPT_LR == 1e-3
    assert T.MONITOR_KS == (0, 1, 2, 3, 5, 10)
    assert T.TRAIN_K == (5, 10, 20)
