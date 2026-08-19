"""Contract tests for the control-arm stage (prereg 39d02166...)."""
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

import controls_train as C  # noqa: E402

PREREG_SHA = "39d02166f69acf235a34d351b649a4cdbf3b828491a0994901bf2378777463f7"


@pytest.fixture(scope="module")
def env():
    return C.load()


def test_prereg_sha_matches_frozen():
    import hashlib
    got = hashlib.sha256((STAGE / "PREREGISTRATION_STAGE1_CONTROLS.md").read_bytes()).hexdigest()
    assert got == PREREG_SHA


def test_random_window_validity(env):
    d1, z1, d2, z2, esm = env
    wt, var, meta = C.build_windows("random_local_window", d1, d2, z2, esm)
    covered = d2["covered_pair_indices"]
    for i in covered:
        q = meta["winpos"][int(i)]
        true_pos = d1["pairs"][i]["pos"]
        assert abs(q - true_pos) > C.R
        # windows differ from the zero sentinel and are matched (same pos for wt/var)
    assert not np.any(np.all(wt[covered] == 0, axis=1))
    # deterministic under the same keyed stream
    wt2, var2, meta2 = C.build_windows("random_local_window", d1, d2, z2, esm)
    assert np.array_equal(wt, wt2) and meta == meta2


def test_family_shuffle_stays_within_parent(env):
    d1, z1, d2, z2, esm = env
    wt, var, meta = C.build_windows("family_preserving_shuffle", d1, d2, z2, esm)
    covered = d2["covered_pair_indices"]
    parent_of_var = {}
    for i in covered:
        parent_of_var[d1["pairs"][i]["parent"]] = parent_of_var.get(
            d1["pairs"][i]["parent"], []) + [i]
    for i in covered:
        # every reassigned window must equal SOME same-parent variant window
        same_parent_pairs = parent_of_var[d1["pairs"][i]["parent"]]
        assert any(np.array_equal(var[i], z2["esm_var"][j])
                   for j in same_parent_pairs), (i, d1["pairs"][i]["parent"])
    # at least one pair changed its window relative to the correct arm
    assert not np.array_equal(var, z2["esm_var"])


def test_random_protein_changes_windows(env):
    d1, z1, d2, z2, esm = env
    wt, var, meta = C.build_windows("random_protein", d1, d2, z2, esm)
    assert not np.array_equal(var, z2["esm_var"])


def test_ligand_only_zero_input(env):
    d1, z1, d2, z2, esm = env
    wt, var, meta = C.build_windows("ligand_only", d1, d2, z2, esm)
    assert np.all(wt == 0) and np.all(var == 0)


def test_ligand_invariant_zero_predictions(env):
    d1, z1, d2, z2, esm = env
    info, model, wt, var, meta = C.train_arm("ligand_invariant_shift",
                                             d1, z1, d2, z2, esm, "cpu")
    assert info["trained"] is False
    rows, agg = C.per_pair_metrics(None, "ligand_invariant_shift", d1, z1, wt, var, "cpu")
    for r in rows:
        assert r["var_pred"] == 0.0 and not r["nonconstant"]
        assert r["spearman"] is None  # undefined, never zero
    assert agg["n_nonconstant"] == 0


def test_arms_share_matched_subset(env):
    d1, z1, d2, z2, esm = env
    covered = d2["covered_pair_indices"]
    sp = np.asarray(d1["split"]["pair_split"])
    assert {k: int((sp[covered] == k).sum()) for k in (0, 1, 2)} == {0: 32, 1: 8, 2: 9}
    for arm in C.ARMS:
        wt, var, meta = C.build_windows(arm, d1, d2, z2, esm)
        assert wt.shape == (65, 640) and var.shape == (65, 640)


def test_effect_uses_observed_point(env):
    d1, z1, d2, z2, esm = env
    def R(parent, r2):
        return {"parent": parent, "r2": r2}
    rows_a = [R("A", 0.5), R("A", 0.5), R("B", 0.1)]
    rows_b = [R("A", 0.0), R("B", 0.0), R("C", 0.0)]
    e = C.effect(rows_a, rows_b, "test")
    assert abs(e["observed_pair_mean_effect"] - 0.3666666667) < 1e-6
    assert "bootstrap_ci" in e and e["bootstrap_ci"]["draws"] == 2000
    assert "leave_one_parent_out_sign_stable" in e


def test_annotation_audit_uses_verified_positions(env):
    d1, z1, d2, z2, esm = env
    a = C.annotation_shortcut_audit(d1, d2, z2, esm)
    assert a["n_pairs"] == 49
    for r in a["per_pair"]:
        assert r["correct_delta_norm"] > 0
        assert r["random_delta_norm"] > 0
