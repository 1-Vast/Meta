"""Contract tests for the read-only collapse audit."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


@pytest.fixture(scope="module")
def audit():
    return json.loads((STAGE / "STAGE1_COLLAPSE_AUDIT.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def frozen():
    art = json.loads((STAGE / "DATA1A.json").read_text(encoding="utf-8"))
    res = json.loads((STAGE / "RESULT_SCREENING.json").read_text(encoding="utf-8"))
    z = np.load(STAGE / "DATA1A.npz", allow_pickle=False)
    return art, res, z


def test_frozen_sha_pins_unchanged(audit, frozen):
    import hashlib
    for name, path in (("RESULT_SCREENING.json", STAGE / "RESULT_SCREENING.json"),
                       ("DATA1A.json", STAGE / "DATA1A.json")):
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        assert audit["frozen_inputs_sha256"][name] == got


def test_structural_identity(audit, frozen):
    art, _, z = frozen
    s = audit["structural_collapse_identity"]
    assert s["constant_prediction_set_equals_zero_input_set"] is True
    assert len(s["test_pairs_with_zero_input_difference"]) == 10
    # verify the zero-difference claim independently
    te = [i for i, x in enumerate(art["split"]["pair_split"]) if x == 2]
    zero = {i for i in te if np.linalg.norm(
        z["prot"][art["pairs"][i]["var_row"]] - z["prot"][art["pairs"][i]["wt_row"]]) == 0}
    assert zero == set(s["test_pairs_with_zero_input_difference"])


def test_38_zero_diff_pairs(audit, frozen):
    art, _, z = frozen
    n = sum(1 for p in art["pairs"] if np.linalg.norm(
        z["prot"][p["var_row"]] - z["prot"][p["wt_row"]]) == 0)
    assert n == 38
    assert audit["variance_sources"]["klifs_diff"]["n_zero_difference_pairs"] == 38
    assert audit["variance_sources"]["klifs_diff"]["zero_difference_by_split"] == {
        "train": 20, "val": 8, "test": 10}


def test_gradient_competition_values(audit):
    g = audit["gradient_audit"]
    assert g["R_g"] > 100
    assert g["C_g"] < 0
    assert g["per_branch"]["g_ctr"]["b_P"] == 0.0
    assert g["per_branch"]["g_ctr"]["b_L"] == 0.0


def test_q1_evidence_recorded(audit):
    q1 = audit["variance_sources"]["q1_evidence"]
    assert q1["pair_centered_local_esm"]["bootstrap_ci_lo"] > 0
    assert q1["klifs_pocket"]["selectivity"] < 0


def test_finite_pairs_counting(audit):
    c = audit["collapse"]["counts"]
    for arm in ("unified_local", "free_pairwise"):
        assert c[arm]["finite_pairs/total_pairs"] == "3/13"
    assert c["ligand_only"]["finite_pairs/total_pairs"] == "0/13"
    # NaN spearman never treated as zero: constant pairs carry nulls
    for r in audit["collapse"]["per_arm"]["ligand_only"]:
        assert r["spearman"] is None
