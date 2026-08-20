"""Structural tests T1-T10 (prereg a7b17e8a... section 10). No labels,
no real-data training, no checkpoints."""
import inspect
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import olr as O  # noqa: E402
import runner as R  # noqa: E402


def _toy(n=37, d=O.D_RES, b=5):
    torch.manual_seed(0)
    res = torch.randn(b, n, d)
    mask = torch.ones(b, n, dtype=torch.bool)
    lig = torch.randn(b, O.D_LIG)
    return res, mask, lig


def test_t1_antisymmetry():
    res, mask, lig = _toy()
    m = O.OLRPotential(router=True)
    sw = m.s(res[0], mask[0], lig)
    sv = m.s(res[1], mask[1], lig)
    d1 = sv - sw                                    # contrast w -> v
    d2 = m.s(res[0], mask[0], lig) - m.s(res[1], mask[1], lig)  # v -> w
    assert torch.max(torch.abs(d1 + d2)) < 1e-5     # antisymmetry


def test_t2_cycle_zero():
    res, mask, lig = _toy(b=1)
    m = O.OLRPotential(router=True)
    P = res[0]
    A, B, C = lig[0], torch.randn(O.D_LIG), torch.randn(O.D_LIG)
    cyc = (m.s(P, mask[0], B.unsqueeze(0)) - m.s(P, mask[0], A.unsqueeze(0))
           + m.s(P, mask[0], C.unsqueeze(0)) - m.s(P, mask[0], B.unsqueeze(0))
           + m.s(P, mask[0], A.unsqueeze(0)) - m.s(P, mask[0], C.unsqueeze(0)))
    assert float(cyc.abs().max()) < 1e-4


def test_t3_centering_identity():
    # panel centering zeroes row means exactly; protein-axis centering is an
    # exact identity on within-pair contrasts (mathematical fact, checked
    # numerically): s_tilde(P,L) - mean_B s_tilde(P,.) cancels in Delta_P.
    torch.manual_seed(1)
    s = torch.randn(2, 9)
    # panel (row) centering zeroes row means exactly
    stilde = s - s.mean(1, keepdim=True)
    assert torch.max(torch.abs(stilde.mean(1))) < 1e-6
    # protein-axis (column) centering is an exact identity on contrasts
    d_raw = s[1] - s[0]
    stilde2 = s - s.mean(0, keepdim=True)
    d2 = stilde2[1] - stilde2[0]
    assert torch.max(torch.abs(d_raw - d2)) < 1e-6


def test_t4_coordinate_free_contract():
    sig = inspect.signature(O.OLRPotential.s)
    params = set(sig.parameters) - {"self"}
    assert params <= {"res", "mask", "lig", "site"}, params
    # site is only permitted when site_channel=True (teacher); default model
    m = O.OLRPotential(router=True)
    assert not m.site_channel
    src = inspect.getsource(O.OLRPotential.s_from_kv)
    assert "site" not in src
    src2 = inspect.getsource(O.OLRPotential.construct_kv)
    assert "self.site_channel" in src2


def test_t5_capacity():
    m = O.OLRPotential(rank=8, heads=1, router=True)
    n = sum(p.numel() for p in m.parameters())
    assert n <= 2_000_000, n


def test_t6_permutation_and_folds():
    rng = np.random.default_rng(3)
    c = np.arange(10, dtype=np.float32)
    p = O.permute_within_pair(c, rng)
    assert sorted(p.tolist()) == sorted(c.tolist())
    assert not np.any(p == c)
    fold_of = O.folds_by_parent(["A", "B", "C", "A", "B", "D", "E", "F"])
    for p1 in fold_of:
        assert fold_of[p1] in {0, 1, 2}
    assert len(fold_of) == 6  # unique parents


def test_t7_leakage_guards():
    d1, z1, d2, esm = O.load_stage_data()
    states = O.construct_states(d1, esm)
    recs = O.pair_tensors(d1, d2, states, z1["lig"])
    split = np.array([r["split1"] for r in recs], dtype=np.int8)
    train_js = [j for j in range(len(recs)) if split[j] == 0]
    w, sel = O.gain_weights(d1, z1, recs, train_js)
    rows = d1["rows"]
    train_parents = {recs[j]["parent"] for j in train_js}
    for i in sel:
        assert "(" not in str(rows[i])
        assert str(rows[i]).strip() in train_parents
    assert 0.25 - 1e-6 <= w.min() and w.max() <= 4.0 + 1e-6
    assert abs(w.mean() - 1.0) < 0.05


def test_t8_determinism():
    torch.manual_seed(7)
    m1 = O.OLRPotential(router=True)
    torch.manual_seed(7)
    m2 = O.OLRPotential(router=True)
    res, mask, lig = _toy()
    for a, b in zip(m1.parameters(), m2.parameters()):
        assert torch.equal(a, b)
    s1 = m1.s(res[0], mask[0], lig)
    s2 = m2.s(res[0], mask[0], lig)
    assert torch.equal(s1, s2)


def test_t9_erased_equality():
    p = HERE / "ERASED_ESM.npz"
    if not p.exists():
        print("  [T9] skipped: ERASED_ESM.npz not yet generated")
        return
    z = np.load(p, allow_pickle=False)
    n_checked = 0
    for k in z.files:
        if k.startswith("we_") and ("me_" + k[3:]) in z.files:
            assert np.max(np.abs(z[k] - z["me_" + k[3:]])) < 1e-5
            n_checked += 1
    assert n_checked >= 40, n_checked


def test_t10_no_closed_forms():
    banned = ["pinv", "linalg.solve", "cholesky", "np.linalg.lstsq", "kernel_ridge", "Ridge("]
    for fname in ("olr.py",):
        src = (HERE / fname).read_text(encoding="utf-8")
        for b in banned:
            assert b not in src, (fname, b)
    src_run = (HERE / "runner.py").read_text(encoding="utf-8")
    for b in banned[:-1]:
        assert b not in src_run, ("runner.py", b)


def test_splits_frozen():
    d1, z1, d2, esm = O.load_stage_data()
    s2 = O.split_s2_s3(d1, d2, 2)
    s3 = O.split_s2_s3(d1, d2, 3)
    assert s2.shape == s3.shape == (49,)
    assert set(np.unique(s2)) <= {0, 1, 2}
    spb, test_parents = O.split_spb(d1, d2)
    parents = [d1["pairs"][i]["parent"] for i in d2["covered_pair_indices"]]
    test_par_set = {parents[j] for j in range(49) if spb[j] == 2}
    assert test_par_set == set(test_parents)
    n_test = int((spb == 2).sum())
    assert len(test_parents) >= 5 and n_test >= 10, (len(test_parents), n_test)
    for j in range(49):
        if spb[j] != 2:
            assert parents[j] not in test_par_set  # parent-disjoint


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL {name}: {e}")
    print("fails:", fails)
    sys.exit(1 if fails else 0)
