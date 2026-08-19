"""Implementation-amendment regression tests for the 2x2 trainer."""
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

import train2x2 as M  # noqa: E402


def test_abs_sampling_differs_across_minibatches():
    pool = np.arange(2000)
    s0 = M.rng_key("esm_joint", "abs", 0, 0).permutation(len(pool))[:512]
    s1 = M.rng_key("esm_joint", "abs", 0, 1).permutation(len(pool))[:512]
    assert not np.array_equal(s0, s1)
    assert len(set(s0) & set(s1)) < 512  # not the identical batch


def test_abs_sampling_deterministic_for_same_key():
    pool = np.arange(2000)
    a = M.rng_key("esm_joint", "abs", 3, 7).permutation(len(pool))[:512]
    b = M.rng_key("esm_joint", "abs", 3, 7).permutation(len(pool))[:512]
    assert np.array_equal(a, b)


def _fixture_rows():
    def R(parent, r2):
        return {"parent": parent, "r2": r2}
    return {
        "esm_joint": [R("A", 0.5), R("A", 0.5), R("B", 0.1)],
        "klifs_joint": [R("A", 0.0), R("B", 0.0), R("C", 0.0)],
        "esm_centered": [R("A", 0.6), R("A", 0.6), R("B", 0.2)],
        "klifs_centered": [R("A", 0.1), R("B", 0.0), R("C", 0.0)],
    }


def test_point_estimate_is_observed_pair_mean():
    fx = _fixture_rows()
    out = M.effect_boot(fx, "test")
    e = out["rep_main_joint"]
    # observed pair-mean: (0.5+0.5+0.1)/3 - 0 = 0.3666...
    assert abs(e["observed_pair_mean_effect"] - 0.3666666667) < 1e-6
    # parent-mean per cell over ITS covered parents: esm_joint {A:0.5,B:0.1}
    # -> 0.3; klifs_joint {A:0,B:0,C:0} -> 0.0; difference = 0.3
    assert abs(e["observed_parent_mean_effect"] - 0.3) < 1e-6
    # the two estimates genuinely differ in this unbalanced fixture
    assert abs(e["observed_pair_mean_effect"] - e["observed_parent_mean_effect"]) > 0.05
    assert "bootstrap_ci" in e and e["bootstrap_ci"]["draws"] == 2000


def _balanced_rows():
    def R(parent, r2):
        return {"parent": parent, "r2": r2}
    return {
        "esm_joint": [R("A", 0.5), R("A", 0.5), R("B", 0.1), R("C", 0.2)],
        "klifs_joint": [R("A", 0.0), R("B", 0.0), R("C", 0.0)],
        "esm_centered": [R("A", 0.6), R("A", 0.6), R("B", 0.2), R("C", 0.3)],
        "klifs_centered": [R("A", 0.1), R("B", 0.0), R("C", 0.05)],
    }


def test_lopo_extends_to_all_five_effects():
    fx = _balanced_rows()
    out = M.effect_boot(fx, "test")
    for k in ("rep_main_joint", "rep_main_centered", "obj_main_klifs",
              "obj_main_esm", "interaction"):
        assert k in out, k
        assert "leave_one_parent_out_sign_stable" in out[k], k
        assert isinstance(out[k]["leave_one_parent_out_sign_stable"], bool)


def test_effect_boot_deterministic():
    fx = _balanced_rows()
    a = M.effect_boot(fx, "test")
    b = M.effect_boot(fx, "test")
    assert a == b


def test_status_rule_uses_observed_point_and_ci():
    fx = _balanced_rows()
    out = M.effect_boot(fx, "test")
    for k, e in out.items():
        lo = e["bootstrap_ci"]["lo2.5"]
        pt = e["observed_pair_mean_effect"]
        expect = ("established" if lo > 0 and abs(pt) >= 0.05
                  else "absent" if abs(pt) < 0.02 else "ambiguous")
        assert e["status"] == expect, k


def test_interaction_is_diff_of_obj_mains():
    fx = _balanced_rows()
    out = M.effect_boot(fx, "test")
    e = out["interaction"]
    expect = (out["obj_main_esm"]["observed_pair_mean_effect"]
              - out["obj_main_klifs"]["observed_pair_mean_effect"])
    assert abs(e["observed_pair_mean_effect"] - expect) < 1e-9
