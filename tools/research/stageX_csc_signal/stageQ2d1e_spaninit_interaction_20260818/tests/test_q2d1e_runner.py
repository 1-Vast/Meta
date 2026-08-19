"""Stage Q2d-1e tests: span-init + L2 regularizer freeze the Q2d-1d
null-space drift fix.

Q2d-1d attribution (frozen): trained correct arm fits train dz 0.99 while
33-43% of the learned protein map energy sits in the train-row null space
(4 unidentifiable directions), degrading cold surfaces. Q2d-1e initializes
A in the feature-only span basis and adds lambda=1e-3 L2 on the factor
maps. These tests pin both mechanisms.
"""
import sys
import json
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent.parent / "stageQ2d1d_spanrestricted_interaction_20260818"
X0C = STAGE_D.parent / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(STAGE_D))
sys.path.insert(0, str(HERE))

import torch  # noqa: E402
import truth_d as truth  # noqa: E402
import runner_e as R  # noqa: E402


@pytest.fixture(scope="module")
def data():
    fz = np.load(STAGE_D / "q2d1d_features.npz", allow_pickle=False)
    splits = json.loads((STAGE_D / "Q2D1D_SPLITS.json").read_text(encoding="utf-8"))
    for k in ("train_cells", "pc", "lc", "dc"):
        splits[k] = np.asarray(splits[k], dtype=np.int64)
    for k in ("cold_row", "cold_lig", "train_row", "train_lig"):
        splits[k] = np.asarray(splits[k], dtype=bool)
    return fz["P_t"].astype(np.float32), fz["L_t"].astype(np.float32), splits


def test_span_init_null_fraction_low(data):
    P_t, L_t, splits = data
    device = "cuda" if torch.cuda.is_available() else "cpu"
    t = truth.generate_truth("M1", 0, P_t, L_t, splits)
    Vsp, _ = truth._span_projection(P_t, splits)
    Lt_dev = torch.from_numpy(L_t).float().to(device)
    model, val, _ = R.train_level(P_t, "correct", t, "A", 0, splits, device, 0,
                                  Lt_dev, Vsp=Vsp, max_steps=600)
    A = model.A.weight.detach().cpu().numpy().T
    null_frac = 1 - float(np.linalg.norm(Vsp @ (Vsp.T @ A)) ** 2 / np.linalg.norm(A) ** 2)
    assert null_frac < 0.20, f"null fraction {null_frac:.3f} >= 0.20 (Q2d-1d was 0.33-0.43)"
    assert val < 1.5


def test_l2_penalty_enters_monitor_loss(data):
    P_t, L_t, splits = data
    device = "cuda" if torch.cuda.is_available() else "cpu"
    t = truth.generate_truth("M1", 0, P_t, L_t, splits)
    Vsp, _ = truth._span_projection(P_t, splits)
    Lt_dev = torch.from_numpy(L_t).float().to(device)
    model, val_pen, _ = R.train_level(P_t, "correct", t, "A", 0, splits, device, 0,
                                      Lt_dev, Vsp=Vsp, max_steps=600)
    old = R.L2_PEN
    R.L2_PEN = 0.0
    try:
        model0, val_nopen, _ = R.train_level(P_t, "correct", t, "A", 0, splits, device, 0,
                                             Lt_dev, Vsp=Vsp, max_steps=600)
    finally:
        R.L2_PEN = old
    nA = float(model.A.weight.square().sum() + model.B.weight.square().sum())
    nB = float(model0.A.weight.square().sum() + model0.B.weight.square().sum())
    assert nA < nB, "L2 penalty must shrink factor-map norms"


def test_oracle_precondition_reads_true():
    pre = json.loads((STAGE_D / "Q2D1D_ORACLE_PRECHECK.json").read_text(encoding="utf-8"))
    assert pre["M1_identifiable_on_all_surfaces"] is True
