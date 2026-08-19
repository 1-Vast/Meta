"""Stage P1 arm 6: FS-CAP-style ligand-only support encoder (AD3 sha
602f8e98...). Ligand-only trunk (mu + l_head(l_enc(xl))), Deep-Sets
support encoder phi_l(2049 -> 64 -> 64 -> 64), mean-pool context,
off_head(64 -> 1, no bias) so k=0 correction is exactly 0. Train: query
MSE per task, AdamW 3e-4. Eval: context encoding only. Artifact
P1_ARM6_FSCAP.json.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import torch
import torch.nn as nn

import p_train as PT
import p_cnp as PC

OUT = PT.OUT
SCHEMA = "MetaSieve.StageP.P1Arm6FSCAP.v1"
AD3_SHA = "5c573132e8e9c82fef78837970ce0e9282d52786ec050e1855d42bf1933dba91"
IN_DIM = 2048 + 1  # ECFP | y


class FS_Trunk(nn.Module):
    """Ligand-only trunk (PTrunk minus protein path and interaction)."""

    def __init__(self):
        super().__init__()
        self.l_enc = nn.Linear(2048, 64)
        self.l_head = nn.Linear(64, 1, bias=False)
        self.mu = nn.Parameter(torch.zeros(1))
        nn.init.xavier_uniform_(self.l_enc.weight)
        nn.init.zeros_(self.l_enc.bias)
        nn.init.xavier_uniform_(self.l_head.weight)

    def forward(self, xl):
        return self.mu + self.l_head(torch.relu(self.l_enc(xl))).squeeze(-1)


class FSCAP(nn.Module):
    def __init__(self):
        super().__init__()
        self.trunk = FS_Trunk()
        self.phi_l = nn.Sequential(nn.Linear(IN_DIM, 64), nn.ReLU(),
                                   nn.Linear(64, 64), nn.ReLU(),
                                   nn.Linear(64, 64))
        self.off_head = nn.Linear(64, 1, bias=False)
        for m in [self.phi_l[0], self.phi_l[2], self.phi_l[4]]:
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)
        nn.init.xavier_uniform_(self.off_head.weight)

    def context(self, xl, y):
        if xl.shape[0] == 0:
            return torch.zeros(1, 64)
        inp = torch.cat([xl, y.unsqueeze(-1)], dim=-1)
        return self.phi_l(inp).mean(0, keepdim=True)

    def forward(self, xl, ctx=None):
        if ctx is None:
            ctx = torch.zeros(1, 64)
        return {"yhat": self.trunk(xl) + self.off_head(ctx).squeeze(-1)}


def build_cell_features_l(xl_ids, pki, lid_of, fps_by_lid):
    xl, y = [], []
    for c in xl_ids:
        xl.append(fps_by_lid[lid_of[c]])
        y.append(pki[c])
    return (np.stack(xl).astype(np.float32),
            np.asarray(y, dtype=np.float32))


def train_seed(seed, device, bank, split_art, pki, lid_of, fps_by_lid, pfeat, dry):
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
    model = FSCAP().to(device)
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=PT.LR, weight_decay=PT.WD)
    best_mon = None
    best_state = None
    n_steps = 2 if dry else PT.STEPS
    for step in range(n_steps):
        rng = PT.stable_rng("stageP", "porder", seed, "step", step)
        kept = PC.collect_tasks(rng, tasks, ligs_of_target, first_cell)
        opt.zero_grad()
        tot = 0.0
        cnt = 0
        for sup, q in kept:
            if not q:
                continue
            xl_s, y_s = build_cell_features_l(sup, pki, lid_of, fps_by_lid)
            ctx = model.context(torch.from_numpy(xl_s).to(device),
                                torch.from_numpy(y_s).to(device))
            xl_q, y_q = build_cell_features_l(q, pki, lid_of, fps_by_lid)
            out = model(torch.from_numpy(xl_q).to(device), ctx=ctx)
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
            mon = monitor(model, bank, pki, lid_of, fps_by_lid, device)
            print(f"fscap seed {seed} step {step + 1} mse {tot / max(cnt, 1):.4f} "
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


def predict(model, sup_ids, q_ids, pki, lid_of, fps_by_lid, device):
    model.eval()
    with torch.no_grad():
        if sup_ids:
            xl_s, y_s = build_cell_features_l(sup_ids, pki, lid_of, fps_by_lid)
            ctx = model.context(torch.from_numpy(xl_s).to(device),
                                torch.from_numpy(y_s).to(device))
        else:
            ctx = None
        xl_q, _ = build_cell_features_l(q_ids, pki, lid_of, fps_by_lid)
        out = model(torch.from_numpy(xl_q).to(device), ctx=ctx)
    return out["yhat"].cpu().numpy()


def monitor(model, bank, pki, lid_of, fps_by_lid, device):
    recs = [r for r in bank["records"] if r["split"] == "p_val" and r["draw"] == 0
            and r["k"] in PT.MONITOR_KS]
    mses = []
    for rec in recs:
        yh = predict(model, rec["support_cell_ids"], rec["query_cell_ids"],
                     pki, lid_of, fps_by_lid, device)
        q_y = np.asarray([pki[c] for c in rec["query_cell_ids"]], dtype=np.float32)
        mses.append(float(np.mean((yh - q_y) ** 2)))
    return float(np.mean(mses))


def eval_seed(model, seed, device, bank, pki, lid_of, fps_by_lid):
    out = {"schema": SCHEMA, "seed": seed, "records": []}
    for rec in bank["records"]:
        yh = predict(model, rec["support_cell_ids"], rec["query_cell_ids"],
                     pki, lid_of, fps_by_lid, device)
        q_y = np.asarray([pki[c] for c in rec["query_cell_ids"]], dtype=np.float32)
        out["records"].append({
            "split": rec["split"], "k": rec["k"], "draw": rec["draw"],
            "target_id": rec["target_id"], "cluster": rec["cluster"],
            "mse": float(np.mean((yh - q_y) ** 2)),
            "centered_mse": float(np.mean(((yh - q_y) - np.mean(yh - q_y)) ** 2)),
            "ci": PT.ci(yh, q_y), "spearman": PT.spearman(yh, q_y),
            "pearson": PT.pearson(yh, q_y), "best_support_loss": None,
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
    fps_by_lid = PT.load_ecfp_cache(lig)
    seeds_out = {}
    for seed in args.seeds:
        model, mon = train_seed(seed, device, bank, split_art, pki, lid_of,
                                fps_by_lid, None, args.dry)
        print(f"fscap seed {seed} final monitor {mon}", flush=True)
        seeds_out[str(seed)] = eval_seed(model, seed, device, bank, pki, lid_of,
                                         fps_by_lid)
    n_trunk = sum(p.numel() for p in PT.PTrunk().parameters())
    n_arm = sum(p.numel() for p in FSCAP().parameters())
    artifact = {
        "schema": SCHEMA,
        "ad3_sha256": AD3_SHA,
        "adaptation_mechanism": "context-encoding (no gradient steps)",
        "note": ("ligand-only arm: no protein features anywhere; different "
                 "adaptation mechanism vs ordinary FT/MAML (50 support "
                 "gradient steps); same data protocol"),
        "param_count_arm3": n_trunk, "param_count_arm6": n_arm,
        "param_delta_vs_arm3": n_arm - n_trunk,
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
    path = OUT / "P1_ARM6_FSCAP.json"
    path.write_text(text, encoding="utf-8")
    art_sha = PT.sha256_file(path)
    (OUT / "P1_ARM6_FSCAP.json.manifest.json").write_text(json.dumps({
        "schema": SCHEMA + ".Manifest", "file": "P1_ARM6_FSCAP.json",
        "sha256": art_sha}, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    for split_name in ("p_val", "p_test"):
        print(f"--- arm6 FSCAP ({split_name}) ---")
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
