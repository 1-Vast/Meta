"""Stage P1 arm 4: first-order MAML (implementation addendum sha
84fea478...). Shared episodic sampler (identical minibatch cells as arm 3
for the first 256 cells), 5 inner SGD steps lr 1e-2 on support, outer
AdamW 3e-4 on query MSE; first-order only. Eval = frozen arm-3 protocol
(50 support steps, best SUPPORT loss). Artifact P1_ARM4_MAML.json.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

import p_train as PT

OUT = PT.OUT
SCHEMA = "MetaSieve.StageP.P1Arm4MAML.v1"
IMPL_SHA = "84fea478382fec1bf07be5a080b222046b18884d121792357e4ebb3418453ccf"
INNER_STEPS = 5
INNER_LR = 1e-2
OUTER_LR = 3e-4
OUTER_WD = 1e-4


def collect_tasks(rng, tasks, ligs_of_target, first_cell):
    """Sample tasks exactly like arm 3's cell accumulation: one task at a
    time until >= BATCH_CELLS cells, keep the first BATCH_CELLS cells."""
    task_list = []
    total = 0
    while total < PT.BATCH_CELLS:
        t = tasks[int(rng.integers(len(tasks)))]
        k = PT.TRAIN_K[int(rng.integers(len(PT.TRAIN_K)))]
        ligs = sorted(ligs_of_target[t])
        perm = rng.permutation(len(ligs))
        ordered = [ligs[i] for i in perm[:k + PT.QUERY]]
        sup = [first_cell[t][l] for l in ordered[:k]]
        q = [first_cell[t][l] for l in ordered[k:]]
        task_list.append((sup, q))
        total += k + PT.QUERY
    kept = []
    n = 0
    for sup, q in task_list:
        need = PT.BATCH_CELLS - n
        if need <= 0:
            break
        if need >= len(sup) + len(q):
            kept.append((sup, q))
            n += len(sup) + len(q)
        elif need > len(sup):
            kept.append((sup, q[:need - len(sup)]))
            n = need
        else:
            kept.append((sup[:need], []))
            n = need
    return kept


def task_fomaml_grad(model, xp_s, xl_s, y_s, xp_q, xl_q, y_q,
                     inner_steps=INNER_STEPS, inner_lr=INNER_LR):
    """Per-task first-order MAML gradient (frozen semantics): INNER_STEPS
    SGD steps on support, then the query-MSE gradient at the ADAPTED
    parameters. The adapted model's gradients are cleared before the query
    backward so the outer gradient is pure d(query)/dw' (first-order MAML),
    not query + last support-step. Returns (qloss, grads aligned with
    model.parameters())."""
    m = copy.deepcopy(model)
    inner = torch.optim.SGD(m.parameters(), lr=inner_lr)
    for _ in range(inner_steps):
        inner.zero_grad()
        loss = ((m(xp_s, xl_s)["yhat"] - y_s) ** 2).mean()
        loss.backward()
        inner.step()
    m.zero_grad()  # FIX (regression test test_p_maml_grad.py): drop the
    # last support-step gradient before the query backward accumulates.
    qloss = ((m(xp_q, xl_q)["yhat"] - y_q) ** 2).mean()
    qloss.backward()
    grads = [p.grad for p in m.parameters()]
    return float(qloss.detach()), grads


def train_seed_maml(seed, device, bank, split_art, pki, lid_of, fps_by_lid, pfeat, dry):
    torch.manual_seed(seed)
    np.random.seed(seed)
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
    model = PT.PTrunk().to(device)
    outer = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                              lr=OUTER_LR, weight_decay=OUTER_WD)
    best_mon = None
    best_state = None
    n_steps = 2 if dry else PT.STEPS
    for step in range(n_steps):
        rng = PT.stable_rng("stageP", "porder", seed, "step", step)
        kept = collect_tasks(rng, tasks, ligs_of_target, first_cell)
        outer.zero_grad()
        qloss_sum = 0.0
        qcount = 0
        for sup, q in kept:
            if not q:
                continue
            xp_s, xl_s, y_s = PT.build_cell_features(sup, pki, lid_of, fps_by_lid,
                                                     pfeat, split_art)
            xp_q, xl_q, y_q = PT.build_cell_features(q, pki, lid_of, fps_by_lid,
                                                     pfeat, split_art)
            qloss, grads = task_fomaml_grad(
                model,
                torch.from_numpy(xp_s).to(device), torch.from_numpy(xl_s).to(device),
                torch.from_numpy(y_s).to(device),
                torch.from_numpy(xp_q).to(device), torch.from_numpy(xl_q).to(device),
                torch.from_numpy(y_q).to(device))
            if not np.isfinite(qloss):
                # inner-loop divergence guard (implementation-level; the
                # frozen inner/outer hyperparameters are untouched): the
                # task contributes no gradient this step.
                continue
            for p_model, g in zip(model.parameters(), grads):
                if p_model.grad is None:
                    p_model.grad = torch.zeros_like(p_model)
                p_model.grad.add_(g)
            qloss_sum += qloss
            qcount += 1
        if qcount > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            outer.step()
        if step % PT.MONITOR_EVERY == PT.MONITOR_EVERY - 1 and not dry:
            mon = PT.monitor_value(model, bank, split_art, pki, lid_of, fps_by_lid,
                                   pfeat, device)
            print(f"maml seed {seed} step {step + 1} qloss {qloss_sum / max(qcount, 1):.4f} "
                  f"monitor {mon:.4f}", flush=True)
            if best_mon is None or mon < best_mon - 1e-9:
                best_mon = mon
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if dry or best_state is None:
        if best_state is None:
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        best_mon = qloss_sum / max(qcount, 1) if best_mon is None else best_mon
    model.load_state_dict(best_state)
    return model, best_mon


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--dry", action="store_true")
    args = ap.parse_args()
    device = args.device
    bank = json.loads((OUT / "P_BANK.json").read_text(encoding="utf-8"))
    split_art = json.loads((OUT / "P_SPLIT.json").read_text(encoding="utf-8"))
    pki, lid_of, lig = PT.load_labels_and_ligands()
    pfeat = PT.load_protein_features()
    fps_by_lid = PT.load_ecfp_cache(lig)
    seeds_out = {}
    for seed in args.seeds:
        model, mon = train_seed_maml(seed, device, bank, split_art, pki, lid_of,
                                     fps_by_lid, pfeat, args.dry)
        print(f"maml seed {seed} final monitor {mon}", flush=True)
        seeds_out[str(seed)] = PT.eval_seed(model, seed, device, bank, split_art, pki,
                                            lid_of, fps_by_lid, pfeat, args.dry)
        seeds_out[str(seed)]["schema"] = SCHEMA
    artifact = {
        "schema": SCHEMA,
        "impl_addendum_sha256": IMPL_SHA,
        "backbone_spec_sha256": PT.BACKBONE_SHA,
        "bank_sha256": PT.sha256_file(OUT / "P_BANK.json"),
        "split_sha256": PT.sha256_file(OUT / "P_SPLIT.json"),
        "dry": bool(args.dry),
        "protocol": {"outer_steps": (2 if args.dry else PT.STEPS),
                     "outer_lr": OUTER_LR, "outer_wd": OUTER_WD,
                     "inner_steps": INNER_STEPS, "inner_lr": INNER_LR,
                     "batch_cells": PT.BATCH_CELLS,
                     "adapt_steps": (2 if args.dry else PT.ADAPT_STEPS),
                     "adapt_lr": PT.ADAPT_LR,
                     "monitor_every": PT.MONITOR_EVERY,
                     "monitor_adapt_steps": PT.MONITOR_ADAPT_STEPS,
                     "train_k": list(PT.TRAIN_K),
                     "min_train_ligs": PT.MIN_TRAIN_LIGS},
        "seeds": seeds_out,
    }
    text = json.dumps(artifact, indent=1, sort_keys=True)
    path = OUT / "P1_ARM4_MAML.json"
    path.write_text(text, encoding="utf-8")
    art_sha = PT.sha256_file(path)
    (OUT / "P1_ARM4_MAML.json.manifest.json").write_text(json.dumps({
        "schema": SCHEMA + ".Manifest", "file": "P1_ARM4_MAML.json",
        "sha256": art_sha}, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    for split_name in ("p_val", "p_test"):
        print(f"--- arm4 MAML ({split_name}) ---")
        for k in PT.K_LIST:
            rows = [r for s in seeds_out.values() for r in s["records"]
                    if r["split"] == split_name and r["k"] == k]
            if not rows:
                continue
            print(f"  k={k:<2} n={len(rows)} MSE={np.mean([r['mse'] for r in rows]):.4f} "
                  f"cMSE={np.mean([r['centered_mse'] for r in rows]):.4f} "
                  f"CI={np.nanmean([r['ci'] for r in rows]):.4f} "
                  f"rho={np.nanmean([r['spearman'] for r in rows]):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
