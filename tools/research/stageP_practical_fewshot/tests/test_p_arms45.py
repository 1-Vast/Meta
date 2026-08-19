"""P1 arms 4/5 invariant tests (CPU; dry artifacts allowed).

1. collect_tasks reproduces arm 3's exact first-256-cell minibatch.
2. MAML dry train + eval structure.
3. CNP AD2 invariants: k=0 exact zero context correction + trunk equality,
   support permutation invariance, query permutation equivariance,
   query-label isolation, input dim, parameter delta.
4. CNP dry train + eval structure.
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
import p_maml as PM  # noqa: E402
import p_cnp as PC  # noqa: E402


@pytest.fixture(scope="module")
def env():
    bank = json.loads((PT.OUT / "P_BANK.json").read_text(encoding="utf-8"))
    split_art = json.loads((PT.OUT / "P_SPLIT.json").read_text(encoding="utf-8"))
    pki, lid_of, lig = PT.load_labels_and_ligands()
    pfeat = PT.load_protein_features()
    fps_by_lid = PT.load_ecfp_cache(lig)
    return dict(bank=bank, split_art=split_art, pki=pki, lid_of=lid_of,
                lig=lig, pfeat=pfeat, fps_by_lid=fps_by_lid)


def _task_struct(split_art, lid_of):
    train_ids = [cid for cid, rec in split_art["cell_split"].items()
                 if rec["split"] == "p_train"]
    ligs_of_target = {}
    cells_of_target = {}
    for c in train_ids:
        t = split_art["cell_split"][c]["target_id"]
        ligs_of_target.setdefault(t, set()).add(lid_of[c])
        cells_of_target.setdefault(t, []).append(c)
    tasks = sorted(t for t in ligs_of_target if len(ligs_of_target[t]) >= PT.MIN_TRAIN_LIGS)
    first_cell = {t: {} for t in tasks}
    for t in tasks:
        for c in cells_of_target[t]:
            first_cell[t].setdefault(lid_of[c], c)
    return tasks, ligs_of_target, first_cell


def test_collect_tasks_matches_arm3_cells(env):
    tasks, ligs_of_target, first_cell = _task_struct(env["split_art"], env["lid_of"])
    for step in (0, 1, 7):
        rng = PT.stable_rng("stageP", "porder", 1, "step", step)
        cells_batch = []
        while len(cells_batch) < PT.BATCH_CELLS:
            t = tasks[int(rng.integers(len(tasks)))]
            k = PT.TRAIN_K[int(rng.integers(len(PT.TRAIN_K)))]
            ligs = sorted(ligs_of_target[t])
            perm = rng.permutation(len(ligs))
            ordered = [ligs[i] for i in perm[:k + PT.QUERY]]
            cells_batch.extend([first_cell[t][l] for l in ordered])
        cells_batch = cells_batch[:PT.BATCH_CELLS]
        rng2 = PT.stable_rng("stageP", "porder", 1, "step", step)
        kept = PM.collect_tasks(rng2, tasks, ligs_of_target, first_cell)
        flat = [c for sup, q in kept for c in sup + q]
        assert len(flat) == PT.BATCH_CELLS == len(cells_batch)
        assert flat == cells_batch


def test_maml_dry_train_and_eval(env):
    model, mon = PM.train_seed_maml(1, "cpu", env["bank"], env["split_art"],
                                    env["pki"], env["lid_of"], env["fps_by_lid"],
                                    env["pfeat"], dry=True)
    assert np.isfinite(mon)
    out = PT.eval_seed(model, 1, "cpu", env["bank"], env["split_art"], env["pki"],
                       env["lid_of"], env["fps_by_lid"], env["pfeat"], dry=True)
    recs = out["records"]
    assert any(r["k"] == 0 and r["best_support_loss"] is None for r in recs)
    for r in recs:
        assert 0 < r["mse"] < 1e6 and 0 <= r["ci"] <= 1


def _rand_cnp_inputs(n_sup=4, n_q=3):
    g = torch.Generator().manual_seed(11)
    xp_s = torch.randn(n_sup, 640, generator=g)
    xl_s = torch.randn(n_sup, 2048, generator=g)
    y_s = torch.randn(n_sup, generator=g)
    xp_q = torch.randn(n_q, 640, generator=g)
    xl_q = torch.randn(n_q, 2048, generator=g)
    return xp_s, xl_s, y_s, xp_q, xl_q


def test_cnp_k0_context_correction_is_exactly_zero():
    torch.manual_seed(0)
    model = PC.CNP()
    xp_s, xl_s, y_s, xp_q, xl_q = _rand_cnp_inputs()
    with torch.no_grad():
        out = model(xp_q, xl_q, ctx=None)
        trunk_out = model.trunk(xp_q, xl_q)["yhat"]
        ctx_zero = model.context(xp_s[:0], xl_s[:0], y_s[:0])
        assert torch.all(ctx_zero == 0.0)
        assert torch.equal(out["yhat"], trunk_out)  # bitwise
        off = model.off_head(ctx_zero).squeeze(-1)
        assert torch.all(off == 0.0)


def test_cnp_support_permutation_invariance():
    torch.manual_seed(1)
    model = PC.CNP()
    xp_s, xl_s, y_s, xp_q, xl_q = _rand_cnp_inputs(n_sup=5, n_q=2)
    with torch.no_grad():
        c1 = model.context(xp_s, xl_s, y_s)
        perm = torch.tensor([4, 0, 2, 1, 3])
        c2 = model.context(xp_s[perm], xl_s[perm], y_s[perm])
        torch.testing.assert_close(c1, c2, atol=1e-6, rtol=1e-6)
        o1 = model(xp_q, xl_q, ctx=c1)["yhat"]
        o2 = model(xp_q, xl_q, ctx=c2)["yhat"]
        torch.testing.assert_close(o1, o2, atol=1e-6, rtol=1e-6)


def test_cnp_query_permutation_equivariance():
    torch.manual_seed(2)
    model = PC.CNP()
    xp_s, xl_s, y_s, xp_q, xl_q = _rand_cnp_inputs(n_sup=4, n_q=3)
    with torch.no_grad():
        ctx = model.context(xp_s, xl_s, y_s)
        o1 = model(xp_q, xl_q, ctx=ctx)["yhat"]
        perm = torch.tensor([2, 0, 1])
        o2 = model(xp_q[perm], xl_q[perm], ctx=ctx)["yhat"]
        torch.testing.assert_close(o1[perm], o2, atol=1e-6, rtol=1e-6)


def test_cnp_query_label_isolation(env):
    torch.manual_seed(3)
    model = PC.CNP()
    rec = env["bank"]["records"][0]
    sup = rec["support_cell_ids"]
    q = rec["query_cell_ids"]
    pki_clean = dict(env["pki"])
    pki_corrupt = dict(env["pki"])
    for c in q:
        pki_corrupt[c] = 999.0  # query labels corrupted
    y_clean = PC.cnp_predict(model, sup, q, pki_clean, env["lid_of"],
                             env["fps_by_lid"], env["pfeat"], env["split_art"], "cpu")
    y_corrupt = PC.cnp_predict(model, sup, q, pki_corrupt, env["lid_of"],
                               env["fps_by_lid"], env["pfeat"], env["split_art"], "cpu")
    np.testing.assert_array_equal(y_clean, y_corrupt)


def test_cnp_input_dim_matches_code():
    assert PC.IN_DIM == 640 + 2048 + 1 == 2689
    model = PC.CNP()
    assert model.phi[0].in_features == PC.IN_DIM


def test_cnp_param_delta_recorded():
    assert PC.param_delta() == 180544


def test_cnp_dry_train_and_eval(env):
    model, mon = PC.train_seed_cnp(1, "cpu", env["bank"], env["split_art"],
                                   env["pki"], env["lid_of"], env["fps_by_lid"],
                                   env["pfeat"], dry=True)
    assert np.isfinite(mon)
    out = PC.eval_seed_cnp(model, 1, "cpu", env["bank"], env["split_art"],
                           env["pki"], env["lid_of"], env["fps_by_lid"], env["pfeat"])
    for r in out["records"]:
        assert 0 < r["mse"] < 1e6
        assert r["best_support_loss"] is None  # context adaptation: no grad steps
