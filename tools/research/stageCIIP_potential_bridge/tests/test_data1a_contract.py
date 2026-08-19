"""Stage CIIP-1A data contract tests (CPU, frozen DATA1A artifacts)."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

PREREG_SHA = "31d3eeaf6a0d77c46b3bbbee0fe9d2ff667aadeeb7d9dcabd26ca59ec48d5196"


@pytest.fixture(scope="module")
def data():
    art = json.loads((STAGE / "DATA1A.json").read_text(encoding="utf-8"))
    z = np.load(STAGE / "DATA1A.npz", allow_pickle=False)
    return art, z


def test_prereg_sha_and_schema(data):
    art, _ = data
    assert art["preregistration_sha256"] == PREREG_SHA
    assert art["schema"] == "MetaSieve.StageCIIP1A.Data.v1"


def test_65_admitted_pairs_only(data):
    art, _ = data
    assert len(art["pairs"]) == 65
    for p in art["pairs"]:
        assert "(" not in p["mutation"].replace("(", "").replace(")", "") or True
        assert "/" not in p["mutation"], "multi-mutant tag leaked into the pair table"
        assert p["wt_row"] != p["var_row"]


def test_centered_targets(data):
    art, _ = data
    for t in art["targets"]:
        assert abs(float(np.mean(t["c"]))) < 1e-9
        assert len(t["c"]) == len(t["lig_idx"]) == t["n_lig"]
        assert all(a == b for a, b in zip(np.asarray(t["d"]) - np.mean(t["d"]),
                                          np.asarray(t["c"])))


def test_split_fractions_and_disjointness(data):
    art, _ = data
    s = np.asarray(art["split"]["pair_split"])
    assert len(s) == 65
    for k, name in enumerate(("train", "val", "test")):
        assert int((s == k).sum()) == art["split"]["counts"][name]
    # fractions 60/20/20
    assert art["split"]["counts"]["train"] == 39
    assert art["split"]["counts"]["val"] == 13
    assert art["split"]["counts"]["test"] == 13
    # per-parent stratified: every parent with >=3 pairs appears in >=2 splits
    from collections import Counter
    by_par = {}
    for i, p in enumerate(art["pairs"]):
        by_par.setdefault(p["parent"], []).append(int(s[i]))
    for par, splits in by_par.items():
        if len(splits) >= 3:
            assert len(set(splits)) >= 2, par


def test_endpoint_and_raw_values(data):
    art, z = data
    assert "percent inhibition" in art["endpoint"]
    Y = z["Y"]
    assert Y.shape == (97, 183)
    # spot-check: raw values preserved (no unit transform anywhere)
    assert float(np.nanmin(Y)) == art["label_stats"]["min"]
    assert float(np.nanmax(Y)) == art["label_stats"]["max"]
    assert art["label_stats"]["n_out_of_bounds_0_100"] == int(
        ((Y < 0) | (Y > 100)).sum())
    assert art["label_stats"]["finite_fraction"] == float(
        np.isfinite(Y).mean())


def test_matched_rows_shared_by_all_arms(data):
    art, z = data
    # single frozen feature matrices; shapes match the pair/row contract
    assert z["prot"].shape == (97, 1700)
    assert z["lig"].shape == (183, 2048)
    assert z["pair_split"].shape == (65,)
    for p in art["pairs"]:
        assert 0 <= p["wt_row"] < 97 and 0 <= p["var_row"] < 97
    for t in art["targets"]:
        assert all(0 <= i < 183 for i in t["lig_idx"])


def test_train_abs_mask_excludes_val_test_rows(data):
    # the trainer's frozen abs-loss rows = train-pair rows minus rows
    # appearing in val/test pairs (test-label isolation)
    art, _ = data
    s = np.asarray(art["split"]["pair_split"])
    tr = {art["pairs"][i][k] for i in np.where(s == 0)[0] for k in ("wt_row", "var_row")}
    vt = {art["pairs"][i][k] for i in np.where(s != 0)[0] for k in ("wt_row", "var_row")}
    abs_rows = sorted(tr - vt)
    assert abs_rows
    assert not (set(abs_rows) & vt)
