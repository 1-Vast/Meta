"""Stage P1 arm 5: deterministic Deep-Sets CNP (addendum AD2 sha
2fde8d6a...). Per-support-item encoder phi (2889 -> 64 -> 64 -> 64),
mean-pool context r (empty support -> r = 0, no trainable prior), decoder
yhat = PTrunk(xp, xl) + off_head(r) with off_head bias=False so the k=0
context correction is exactly 0 and k=0 equals the shared trunk output.
Train: per-task query MSE (no KL, no sampling). Eval: context encoding
only, support-only, query labels never enter. Artifact P1_ARM5_CNP.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
import torch
import torch.nn as nn

import p_train as PT

OUT = PT.OUT
SCHEMA = "MetaSieve.StageP.P1Arm5CNP.v2"
IMPL_SHA = "84fea478382fec1bf07be5a080b222046b18884d121792357e4ebb3418453ccf"
AD2_SHA = "f8909eded8d3d11ec8e0207fdec9a30873ba383023cd57838c9ec4d462069e74"
IN_DIM = 640 + 2048 + 1  # protein | ECFP | y (documented in AD2)


class CNP(nn.Module):
    """Deterministic Deep-Sets CNP (AD2); parameter delta vs arm 3 =
    180,544."""

    def __init__(self):
        super().__init__()
        self.trunk = PT.PTrunk()
        self.phi = nn.Sequential(nn.Linear(IN_DIM, 64), nn.ReLU(),
                                 nn.Linear(64, 64), nn.ReLU(),
                                 nn.Linear(64, 64))
        self.off_head = nn.Linear(64, 1, bias=False)  # k=0 correction == 0
        for m in [self.phi[0], self.phi[2], self.phi[4]]:
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)
        nn.init.xavier_uniform_(self.off_head.weight)

    def context(self, xp, xl, y):
        """r = mean_i phi([xp_i | xl_i | y_i]); empty support -> fixed 0."""
        if xp.shape[0] == 0:
            return torch.zeros(1, 64)
        inp = torch.cat([xp, xl, y.unsqueeze(-1)], dim=-1)
        return self.phi(inp).mean(0, keepdim=True)

    def forward(self, xp, xl, ctx=None):
        out = self.trunk(xp, xl)
        if ctx is None:
            ctx = torch.zeros(1, 64)
        return {"yhat": out["yhat"] + self.off_head(ctx).squeeze(-1)}


def param_delta():
    cnp = CNP()
    trunk = PT.PTrunk()
    n_cnp = sum(p.numel() for p in cnp.parameters())
    n_trunk = sum(p.numel() for p in trunk.parameters())
    return n_cnp - n_trunk


def collect_tasks(rng, tasks, ligs_of_target, first_cell):
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


def train_seed_cnp(seed, device, bank, split_art, pki, lid_of, fps_by_lid, pfeat, dry):
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
    model = CNP().to(device)
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=PT.LR, weight_decay=PT.WD)
    best_mon = None
    best_state = None
    n_steps = 2 if dry else PT.STEPS
    for step in range(n_steps):
        rng = PT.stable_rng("stageP", "porder", seed, "step", step)
        kept = collect_tasks(rng, tasks, ligs_of_target, first_cell)
        opt.zero_grad()
        tot = 0.0
        cnt = 0
        for sup, q in kept:
            if not q:
                continue
            xp_s, xl_s, y_s = PT.build_cell_features(sup, pki, lid_of, fps_by_lid,
                                                     pfeat, split_art)
            ctx = model.context(torch.from_numpy(xp_s).to(device),
                                torch.from_numpy(xl_s).to(device),
                                torch.from_numpy(y_s).to(device))
            xp_q, xl_q, y_q = PT.build_cell_features(q, pki, lid_of, fps_by_lid,
                                                     pfeat, split_art)
            out = model(torch.from_numpy(xp_q).to(device),
                        torch.from_numpy(xl_q).to(device), ctx=ctx)
            mse = ((out["yhat"] - torch.from_numpy(y_q).to(device)) ** 2).mean()
            if not torch.isfinite(mse):
                continue
            mse.backward()
            tot += float(mse.detach())
            cnt += 1
        if cnt > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            opt.step()
        if step % PT.MONITOR_EVERY == PT.MONITOR_EVERY - 1 and not dry:
            mon = monitor_cnp(model, bank, split_art, pki, lid_of, fps_by_lid,
                              pfeat, device)
            print(f"cnp seed {seed} step {step + 1} mse {tot / max(cnt, 1):.4f} "
                  f"monitor {mon:.4f}", flush=True)
            if best_mon is None or mon < best_mon - 1e-9:
                best_mon = mon
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if dry or best_state is None:
        if best_state is None:
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        best_mon = tot / max(cnt, 1) if best_mon is None else best_mon
    model.load_state_dict(best_state)
    return model, best_mon


def cnp_predict(model, sup_ids, q_ids, pki, lid_of, fps_by_lid, pfeat,
                split_art, device):
    """Eval adaptation = context encoding (AD2): support-only by
    construction; query labels never enter."""
    model.eval()
    with torch.no_grad():
        if sup_ids:
            xp_s, xl_s, y_s = PT.build_cell_features(sup_ids, pki, lid_of,
                                                     fps_by_lid, pfeat, split_art)
            ctx = model.context(torch.from_numpy(xp_s).to(device),
                                torch.from_numpy(xl_s).to(device),
                                torch.from_numpy(y_s).to(device))
        else:
            ctx = None
        xp, xl, _ = PT.build_cell_features(q_ids, pki, lid_of, fps_by_lid,
                                           pfeat, split_art)
        out = model(torch.from_numpy(xp).to(device),
                    torch.from_numpy(xl).to(device), ctx=ctx)
    return out["yhat"].cpu().numpy()


def monitor_cnp(model, bank, split_art, pki, lid_of, fps_by_lid, pfeat, device):
    recs = [r for r in bank["records"] if r["split"] == "p_val" and r["draw"] == 0
            and r["k"] in PT.MONITOR_KS]
    mses = []
    for rec in recs:
        yh = cnp_predict(model, rec["support_cell_ids"], rec["query_cell_ids"],
                         pki, lid_of, fps_by_lid, pfeat, split_art, device)
        q_y = np.asarray([pki[c] for c in rec["query_cell_ids"]], dtype=np.float32)
        mses.append(float(np.mean((yh - q_y) ** 2)))
    return float(np.mean(mses))


def eval_seed_cnp(model, seed, device, bank, split_art, pki, lid_of, fps_by_lid, pfeat):
    out = {"schema": SCHEMA, "seed": seed, "records": []}
    for rec in bank["records"]:
        yh = cnp_predict(model, rec["support_cell_ids"], rec["query_cell_ids"],
                         pki, lid_of, fps_by_lid, pfeat, split_art, device)
        q_y = np.asarray([pki[c] for c in rec["query_cell_ids"]], dtype=np.float32)
        out["records"].append({
            "split": rec["split"], "k": rec["k"], "draw": rec["draw"],
            "target_id": rec["target_id"], "cluster": rec["cluster"],
            "mse": float(np.mean((yh - q_y) ** 2)),
            "centered_mse": float(np.mean(((yh - q_y) - np.mean(yh - q_y)) ** 2)),
            "ci": PT.ci(yh, q_y), "spearman": PT.spearman(yh, q_y),
            "pearson": PT.pearson(yh, q_y), "best_support_loss": None,
            "yhat": yh.tolist(), "y": q_y.tolist(),
        })
    return out


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
        model, mon = train_seed_cnp(seed, device, bank, split_art, pki, lid_of,
                                    fps_by_lid, pfeat, args.dry)
        print(f"cnp seed {seed} final monitor {mon}", flush=True)
        seeds_out[str(seed)] = eval_seed_cnp(model, seed, device, bank, split_art,
                                             pki, lid_of, fps_by_lid, pfeat)
    artifact = {
        "schema": SCHEMA,
        "impl_addendum_sha256": IMPL_SHA,
        "ad2_sha256": AD2_SHA,
        "adaptation_mechanism": "context-encoding (no gradient steps)",
        "note": ("different adaptation mechanism vs ordinary FT/MAML "
                 "(50 support gradient steps); same data protocol"),
        "param_delta_vs_arm3": param_delta(),
        "backbone_spec_sha256": PT.BACKBONE_SHA,
        "bank_sha256": PT.sha256_file(OUT / "P_BANK.json"),
        "split_sha256": PT.sha256_file(OUT / "P_SPLIT.json"),
        "dry": bool(args.dry),
        "protocol": {"steps": (2 if args.dry else PT.STEPS), "lr": PT.LR, "wd": PT.WD,
                     "batch_cells": PT.BATCH_CELLS, "in_dim": IN_DIM,
                     "eval_adaptation": "context-encoding-only",
                     "monitor_every": PT.MONITOR_EVERY,
                     "monitor_adapt_steps": PT.MONITOR_ADAPT_STEPS,
                     "train_k": list(PT.TRAIN_K),
                     "min_train_ligs": PT.MIN_TRAIN_LIGS},
        "seeds": seeds_out,
    }
    text = json.dumps(artifact, indent=1, sort_keys=True)
    path = OUT / "P1_ARM5_CNP.json"
    path.write_text(text, encoding="utf-8")
    art_sha = PT.sha256_file(path)
    (OUT / "P1_ARM5_CNP.json.manifest.json").write_text(json.dumps({
        "schema": SCHEMA + ".Manifest", "file": "P1_ARM5_CNP.json",
        "sha256": art_sha}, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    for split_name in ("p_val", "p_test"):
        print(f"--- arm5 CNP ({split_name}) ---")
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
