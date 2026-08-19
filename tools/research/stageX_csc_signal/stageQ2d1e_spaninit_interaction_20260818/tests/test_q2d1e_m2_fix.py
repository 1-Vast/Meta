"""Q2d-1e AD1 regression tests: truth_e repairs (CPU-only).

Asserts M1/M2/M3 streams bit-identical to frozen truth_d; NC1 returns the
frozen description (I = 0, A/B None); NC2 returns A/B None; runner_e's
oracle arm uses the zero bound for NC mechanisms.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
STAGE_D = STAGE.parent / "stageQ2d1d_spanrestricted_interaction_20260818"
sys.path.insert(0, str(STAGE))
sys.path.insert(0, str(STAGE_D))

import truth_d  # noqa: E402
import truth_e  # noqa: E402
import runner_e as R  # noqa: E402


@pytest.fixture(scope="module")
def data():
    fz = np.load(STAGE_D / "q2d1d_features.npz", allow_pickle=False)
    splits = json.loads((STAGE_D / "Q2D1D_SPLITS.json").read_text(encoding="utf-8"))
    for k in ("train_cells", "pc", "lc", "dc"):
        splits[k] = np.asarray(splits[k], dtype=np.int64)
    for k in ("cold_row", "cold_lig", "train_row", "train_lig"):
        splits[k] = np.asarray(splits[k], dtype=bool)
    truth_d.PCA_VT = fz["PCA_VT"].astype(np.float64)
    return (fz["P_t"].astype(np.float32), fz["L_t"].astype(np.float32),
            splits)


@pytest.mark.parametrize("mech,seed", [("M1", 0), ("M1", 1), ("M2", 0),
                                       ("M2", 1), ("M3", 0), ("M3", 1)])
def test_m_streams_bit_identical(data, mech, seed):
    P_t, L_t, splits = data
    a = truth_d.generate_truth(mech, seed, P_t, L_t, splits)
    b = truth_e.generate_truth(mech, seed, P_t, L_t, splits)
    for key in a:
        va, vb = a[key], b[key]
        if va is None and vb is None:
            continue
        np.testing.assert_array_equal(np.asarray(va), np.asarray(vb),
                                      err_msg=f"{mech} {seed} {key}")


def test_nc1_frozen_description(data):
    P_t, L_t, splits = data
    t = truth_e.generate_truth("NC1", 0, P_t, L_t, splits)
    assert t["A"] is None and t["B"] is None
    assert np.all(np.asarray(t["I"]) == 0.0)
    assert t["mu"] == 0.5


def test_nc2_no_feature_map(data):
    P_t, L_t, splits = data
    t = truth_e.generate_truth("NC2", 0, P_t, L_t, splits)
    assert t["A"] is None and t["B"] is None
    assert np.isfinite(np.asarray(t["I"])).all()
    assert np.asarray(t["I"]).shape == (P_t.shape[0], L_t.shape[0])


def test_oracle_zero_bound_for_nc(data):
    P_t, L_t, splits = data
    t = truth_e.generate_truth("NC1", 0, P_t, L_t, splits)
    shuf = np.arange(len(P_t))
    ai = R.build_arm_inputs(P_t, t, shuf, shuf, np.zeros_like(P_t))
    assert np.all(ai["oracle_diagnostic"] == 0.0)
    t2 = truth_e.generate_truth("M1", 0, P_t, L_t, splits)
    ai2 = R.build_arm_inputs(P_t, t2, shuf, shuf, np.zeros_like(P_t))
    assert not np.all(ai2["oracle_diagnostic"] == 0.0)
