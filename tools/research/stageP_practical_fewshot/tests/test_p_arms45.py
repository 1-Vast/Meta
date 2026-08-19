"""P1 arms 4/5 invariant tests (CPU; dry artifacts allowed).

1. collect_tasks reproduces arm 3's exact first-256-cell minibatch for a
   given (seed, step) — the shared frozen sampler.
2. MAML dry train + eval structure.
3. CNP k=0 prior path (frozen bypass assertion).
4. CNP dry train + eval structure + parameter delta vs PTrunk.
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
        # arm-3 accumulation (inline copy of p_train.train_seed)
        cells_batch = []
        while len(cells_batch) < PT.BATCH_CELLS:
            t = tasks[int(rng.integers(len(tasks)))]
            k = PT.TRAIN_K[int(rng.integers(len(PT.TRAIN_K)))]
            ligs = sorted(ligs_of_target[t])
            perm = rng.permutation(len(ligs))
            ordered = [ligs[i] for i in perm[:k + PT.QUERY]]
            cells_batch.extend([first_cell[t][l] for l in ordered])
        cells_batch = cells_batch[:PT.BATCH_CELLS]
        # collect_tasks (fresh rng, same key)
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


def test_cnp_k0_prior_path():
    torch.manual_seed(0)
    model = PC.CNP()
    xp = torch.randn(3, 640)
    xl = torch.randn(3, 2048)
    out = model(xp, xl, ctx=None)
    with torch.no_grad():
        trunk_out = model.trunk(xp, xl)["yhat"]
        manual = trunk_out + model.off_head(model.prior).squeeze(-1)
    assert torch.allclose(out["yhat"], manual, atol=1e-6)


def test_cnp_param_delta_positive():
    cnp = PC.CNP()
    trunk = PT.PTrunk()
    n_cnp = sum(p.numel() for p in cnp.parameters())
    n_trunk = sum(p.numel() for p in trunk.parameters())
    assert n_cnp == n_trunk + sum(p.numel() for p in
                                  list(cnp.enc.parameters()) +
                                  list(cnp.mu_head.parameters()) +
                                  list(cnp.lv_head.parameters()) +
                                  list(cnp.off_head.parameters()) +
                                  [cnp.prior])


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
