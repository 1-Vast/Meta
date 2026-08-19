"""P2-arms structure tests (CPU; dry artifacts allowed)."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

import p2_arms as PA  # noqa: E402


def test_arm_metrics_structure():
    art = json.loads((PA.OUT / "P1_ARM3_ORDINARYFT.json").read_text(encoding="utf-8"))
    rec = art["seeds"]["1"]["records"][0]
    assert "yhat" in rec and "y" in rec and len(rec["yhat"]) == len(rec["y"])
    m = PA.arm_metrics(PA.OUT / "P1_ARM3_ORDINARYFT.json")
    assert "p_test:k5" in m and "p_test:k0" in m
    for key, row in m.items():
        assert row["n_records"] > 0
        assert 0 <= row["pr_auc"] <= 1
        assert row["ef1"] > 0 and np.isfinite(row["ef1"])
        assert 0 <= row["top1_hit"] <= 1


def test_pr_auc_and_ef_reference():
    yhat = [9.0, 8.0, 7.0, 6.0]
    y = [9.0, 8.0, 5.0, 5.0]
    assert PA.pr_auc(yhat, y) == 1.0
    # top-50% contains ALL actives (base rate 50%) -> enrichment 2.0
    assert PA.ef(yhat, y, 0.5) == 2.0
    assert PA.top1_hit(yhat, y) == 1.0
    # active ranked last among 2 -> average precision 1/2
    assert PA.pr_auc([6.0, 5.0], [5.0, 6.0]) == 0.5
