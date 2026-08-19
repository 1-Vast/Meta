"""P1 arms 6/7 invariant tests (CPU; dry artifacts allowed).

Arm 6 (FS-CAP-style): k=0 zero correction + trunk equality, support
permutation invariance, query equivariance, query-label isolation,
input dim 2049, param counts. Arm 7 (ActFound-style): identity-zero +
antisymmetry bitwise, k=0 constant, anchor permutation invariance,
query-label isolation, dry train + eval structure for both.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

import torch  # noqa: E402
import p_train as PT  # noqa: E402
import p_fscap as PF  # noqa: E402
import p_actfound as PA  # noqa: E402


@pytest.fixture(scope="module")
def env():
    bank = json.loads((PT.OUT / "P_BANK.json").read_text(encoding="utf-8"))
    split_art = json.loads((PT.OUT / "P_SPLIT.json").read_text(encoding="utf-8"))
    pki, lid_of, lig = PT.load_labels_and_ligands()
    pfeat = PT.load_protein_features()
    fps_by_lid = PT.load_ecfp_cache(lig)
    return dict(bank=bank, split_art=split_art, pki=pki, lid_of=lid_of,
                lig=lig, pfeat=pfeat, fps_by_lid=fps_by_lid)


# ---------------- arm 6 ----------------

def _rand_l(n_sup=4, n_q=3):
    g = torch.Generator().manual_seed(21)
    xl_s = torch.randn(n_sup, 2048, generator=g)
    y_s = torch.randn(n_sup, generator=g)
    xl_q = torch.randn(n_q, 2048, generator=g)
    return xl_s, y_s, xl_q


def test_fscap_k0_zero_correction_and_trunk_equality():
    torch.manual_seed(0)
    model = PF.FSCAP()
    xl_s, y_s, xl_q = _rand_l()
    with torch.no_grad():
        out = model(xl_q, ctx=None)
        trunk = model.trunk(xl_q)
        ctx0 = model.context(xl_s[:0], y_s[:0])
        assert torch.all(ctx0 == 0.0)
        assert torch.equal(out["yhat"], trunk)
        assert torch.all(model.off_head(ctx0).squeeze(-1) == 0.0)


def test_fscap_support_permutation_invariance_and_query_equivariance():
    torch.manual_seed(1)
    model = PF.FSCAP()
    xl_s, y_s, xl_q = _rand_l(n_sup=5, n_q=3)
    with torch.no_grad():
        c1 = model.context(xl_s, y_s)
        perm = torch.tensor([3, 0, 4, 1, 2])
        c2 = model.context(xl_s[perm], y_s[perm])
        torch.testing.assert_close(c1, c2, atol=1e-6, rtol=1e-6)
        o1 = model(xl_q, ctx=c1)["yhat"]
        qperm = torch.tensor([2, 0, 1])
        o2 = model(xl_q[qperm], ctx=c1)["yhat"]
        torch.testing.assert_close(o1[qperm], o2, atol=1e-6, rtol=1e-6)


def test_fscap_query_label_isolation(env):
    torch.manual_seed(2)
    model = PF.FSCAP()
    rec = env["bank"]["records"][0]
    sup, q = rec["support_cell_ids"], rec["query_cell_ids"]
    pki_clean = dict(env["pki"])
    pki_corrupt = dict(env["pki"])
    for c in q:
        pki_corrupt[c] = 999.0
    y1 = PF.predict(model, sup, q, pki_clean, env["lid_of"], env["fps_by_lid"], "cpu")
    y2 = PF.predict(model, sup, q, pki_corrupt, env["lid_of"], env["fps_by_lid"], "cpu")
    np.testing.assert_array_equal(y1, y2)


def test_fscap_input_dim_and_param_counts():
    assert PF.IN_DIM == 2048 + 1 == 2049
    n_trunk = sum(p.numel() for p in PT.PTrunk().parameters())
    n_arm = sum(p.numel() for p in PF.FSCAP().parameters())
    assert n_trunk == 174339
    assert n_arm - n_trunk == 96446


def test_fscap_dry_train_and_eval(env):
    model, mon = PF.train_seed(1, "cpu", env["bank"], env["split_art"], env["pki"],
                               env["lid_of"], env["fps_by_lid"], None, dry=True)
    assert np.isfinite(mon)
    out = PF.eval_seed(model, 1, "cpu", env["bank"], env["pki"], env["lid_of"],
                       env["fps_by_lid"])
    for r in out["records"]:
        assert 0 < r["mse"] < 1e6
        assert r["best_support_loss"] is None


# ---------------- arm 7 ----------------

def _rand_pair(n=3):
    g = torch.Generator().manual_seed(31)
    xp_a = torch.randn(n, 640, generator=g)
    xl_a = torch.randn(n, 2048, generator=g)
    xp_b = torch.randn(n, 640, generator=g)
    xl_b = torch.randn(n, 2048, generator=g)
    return xp_a, xl_a, xp_b, xl_b


def test_pairnet_identity_zero_and_antisymmetry_bitwise():
    torch.manual_seed(0)
    model = PA.PairNet()
    xp_a, xl_a, xp_b, xl_b = _rand_pair()
    with torch.no_grad():
        same = model.d(xp_a, xl_a, xp_a, xl_a)
        assert torch.all(same == 0.0)
        fwd = model.d(xp_a, xl_a, xp_b, xl_b)
        bwd = model.d(xp_b, xl_b, xp_a, xl_a)
        assert torch.equal(fwd, -bwd)


def test_actfound_k0_constant_and_anchor_invariance(env):
    torch.manual_seed(1)
    model = PA.PairNet()
    rec = env["bank"]["records"][0]
    q = rec["query_cell_ids"]
    k0 = 3.14
    y0 = PA.predict(model, [], q, env["pki"], env["lid_of"], env["fps_by_lid"],
                    env["pfeat"], env["split_art"], "cpu", k0)
    np.testing.assert_array_equal(y0, np.full(len(q), k0, dtype=np.float32))
    sup = rec["support_cell_ids"]
    y1 = PA.predict(model, sup, q, env["pki"], env["lid_of"], env["fps_by_lid"],
                    env["pfeat"], env["split_art"], "cpu", k0)
    y2 = PA.predict(model, list(reversed(sup)), q, env["pki"], env["lid_of"],
                    env["fps_by_lid"], env["pfeat"], env["split_art"], "cpu", k0)
    np.testing.assert_allclose(y1, y2, atol=1e-6)


def test_actfound_query_label_isolation(env):
    torch.manual_seed(2)
    model = PA.PairNet()
    rec = env["bank"]["records"][0]
    sup, q = rec["support_cell_ids"], rec["query_cell_ids"]
    pki_clean = dict(env["pki"])
    pki_corrupt = dict(env["pki"])
    for c in q:
        pki_corrupt[c] = 999.0
    y1 = PA.predict(model, sup, q, pki_clean, env["lid_of"], env["fps_by_lid"],
                    env["pfeat"], env["split_art"], "cpu", 3.14)
    y2 = PA.predict(model, sup, q, pki_corrupt, env["lid_of"], env["fps_by_lid"],
                    env["pfeat"], env["split_art"], "cpu", 3.14)
    np.testing.assert_array_equal(y1, y2)


def test_actfound_dry_train_and_eval(env):
    model, mon = PA.train_seed(1, "cpu", env["bank"], env["split_art"], env["pki"],
                               env["lid_of"], env["fps_by_lid"], env["pfeat"], dry=True)
    assert np.isfinite(mon)
    out = PA.eval_seed(model, 1, "cpu", env["bank"], env["split_art"], env["pki"],
                       env["lid_of"], env["fps_by_lid"], env["pfeat"])
    for r in out["records"]:
        assert 0 < r["mse"] < 1e6
    k0recs = [r for r in out["records"] if r["k"] == 0]
    assert k0recs and all(np.isclose(r["mse"],
                                     np.mean((np.full(8, out["k0_mean"]) -
                                              np.array([env["pki"][c] for c in
                                                        env["bank"]["records"][
                                                            out["records"].index(r)][
                                                            "query_cell_ids"]])) ** 2),
                                     atol=1e-5)
                           for r in k0recs[:1])
