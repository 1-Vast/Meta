"""Stage Q2d-1d runner tests (forensic + regression).

The misalignment regression: the original loss_fn paired minibatch outputs
with the first n targets in train-set order. With rng-drawn minibatch
indices the gradients were noise and every arm collapsed to the zero
predictor (monitor loss ~= 2.0 = z variance). Fixed: loss_fn receives the
minibatch-aligned targets. These tests freeze the fix.
"""
import sys
import json
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
X0C = STAGE.parent / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(STAGE))

import torch  # noqa: E402
import truth_d as truth  # noqa: E402
import runner_d as R  # noqa: E402


@pytest.fixture(scope="module")
def data():
    fz = np.load(STAGE / "q2d1d_features.npz", allow_pickle=False)
    splits = json.loads((STAGE / "Q2D1D_SPLITS.json").read_text(encoding="utf-8"))
    for k in ("train_cells", "pc", "lc", "dc"):
        splits[k] = np.asarray(splits[k], dtype=np.int64)
    for k in ("cold_row", "cold_lig", "train_row", "train_lig"):
        splits[k] = np.asarray(splits[k], dtype=bool)
    return fz["P_t"].astype(np.float32), fz["L_t"].astype(np.float32), splits


def test_span_rank_and_truth_in_span(data):
    P_t, L_t, splits = data
    Xtr = P_t[splits["train_row"]].astype(np.float64)
    _, S, _ = np.linalg.svd(Xtr, full_matrices=False)
    r = int((S > 1e-6).sum())
    assert r < P_t.shape[1], "premise: train rows under-basis the features"
    t = truth.generate_truth("M1", 0, P_t, L_t, splits)
    V, r2 = truth._span_projection(P_t, splits)
    assert r2 == r
    resid = t["A"] - V @ (V.T @ t["A"])
    assert np.abs(resid).max() < 1e-9, "M1 truth map must lie in the train-row span"


def test_oracle_precheck_json_pass():
    pre = json.loads((STAGE / "Q2D1D_ORACLE_PRECHECK.json").read_text(encoding="utf-8"))
    assert pre["M1_identifiable_on_all_surfaces"] is True
    assert pre["STOP_before_training"] is False


def test_minibatch_target_alignment(data):
    """Regression: with aligned targets the learner fits the signal (monitor
    loss well below the zero-predictor loss ~2.0); with the old misaligned
    pairing it stayed at chance."""
    P_t, L_t, splits = data
    device = "cuda" if torch.cuda.is_available() else "cpu"
    t = truth.generate_truth("M1", 0, P_t, L_t, splits)
    P_or = (P_t.astype(np.float64) @ t["A"]).astype(np.float32)
    Lt_dev = torch.from_numpy(L_t).float().to(device)
    model, val, _ = R.train_level(P_or, "oracle_diagnostic", t, "A", 0, splits,
                                  device, 0, Lt_dev, max_steps=600)
    assert val < 1.5, "monitor loss {0:.3f} - learner failed to fit aligned targets".format(val)
    res = R.eval_arm(model, P_or, "oracle_diagnostic", t, "A", splits, device, Lt_dev)
    assert res["dc"]["dz"] > 0.75, "dc dz {0:.3f}".format(res["dc"]["dz"])
