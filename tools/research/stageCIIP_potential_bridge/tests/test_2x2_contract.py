"""Contract tests for the 2x2 matched-subset data and frozen gates."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
PREREG2X2_SHA = "ee844b2b29f0009cf97c1bd18b8a92f68dcd2dc8ea1268731c740c2224f47a8b"


@pytest.fixture(scope="module")
def d2():
    art = json.loads((STAGE / "DATA2X2.json").read_text(encoding="utf-8"))
    z = np.load(STAGE / "DATA2X2.npz", allow_pickle=False)
    d1 = json.loads((STAGE / "DATA1A.json").read_text(encoding="utf-8"))
    return art, z, d1


def test_prereg_sha(d2):
    import hashlib
    got = hashlib.sha256((STAGE / "PREREGISTRATION_STAGE1_2X2.md").read_bytes()).hexdigest()
    assert got == PREREG2X2_SHA


def test_matched_49_subset(d2):
    art, z, d1 = d2
    assert art["schema"] == "MetaSieve.StageCIIP1A.2x2.Data.v1"
    covered = list(art["covered_pair_indices"])
    assert len(covered) == 49
    assert z["esm_wt"].shape == (65, 640) and z["esm_var"].shape == (65, 640)
    # missing rows are exactly zero; covered rows are all nonzero
    assert not np.any(np.all(z["esm_wt"][covered] == 0, axis=1))
    rest = [i for i in range(65) if i not in covered]
    assert np.all(z["esm_wt"][rest] == 0)
    # all four cells share the same subset: contract stored once
    sp = np.asarray(d1["split"]["pair_split"])
    cnt = {k: int((sp[covered] == k).sum()) for k in (0, 1, 2)}
    assert cnt == {0: 32, 1: 8, 2: 9}
    assert art["coverage_bias"]["covered"]["split_counts"] == {"train": 32, "val": 8, "test": 9}


def test_oracle_semantics(d2):
    art, _, _ = d2
    o = art["oracle_local_esm"]
    assert "ORACLE" in o["note"] or "oracle" in o["note"]
    assert o["window_radius"] == 6 and o["dim"] == 640
    assert "NOT a deployable" in o["note"]
    assert "oracle-covered subset diagnostic" in art["coverage_bias"]["verdict"]


def test_covered_test_pairs(d2):
    art, _, d1 = d2
    sp = np.asarray(d1["split"]["pair_split"])
    expect = sorted(i for i in art["covered_pair_indices"] if sp[i] == 2)
    assert art["covered_test_pairs"] == expect
    assert art["n_covered_test_parents"] == 6
    assert len(expect) == 9


def test_missing_bias_recorded(d2):
    art, _, _ = d2
    b = art["coverage_bias"]
    assert b["missing"]["n_pairs"] == 16
    assert b["missing"]["n_parents"] == 4
    assert all(x["reason"] == "pos > ESM_MAX_LEN" for x in b["missing_detail"])
    assert b["covered"]["var_true_median"] < b["missing"]["var_true_median"]
    assert b["covered"]["informative_frac_median"] < b["missing"]["informative_frac_median"]


def test_endpoint_language(d2):
    art, _, _ = d2
    assert "percent inhibition" in art["endpoint"]
    for bad in ("pK", "Ki", "Kd", "affinity"):
        assert bad not in art["endpoint"]
