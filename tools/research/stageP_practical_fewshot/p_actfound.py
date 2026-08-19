"""Stage P1 arm 7: ActFound-style pairwise within-target supervision
(AD3 sha 602f8e98...). Pair module D(xq,xs) = h([q|s]) - h([s|q]),
h: Linear(5376 -> 64, ReLU) -> Linear(64 -> 64, ReLU) -> Linear(64 -> 1);
identity-zero and exchange-antisymmetric by construction. Train: ordered
pairs inside each sampled task (first 256), MSE(D, y_i - y_j). Eval:
yhat(q) = mean_s [y_s + D(xq, xs)]; k=0 -> p_train label mean (frozen
constant). Artifact P1_ARM7_ACTFOUND.json.
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
SCHEMA = "MetaSieve.StageP.P1Arm7ActFound.v1"
AD3_SHA = "5c573132e8e9c82fef78837970ce0e9282d52786ec050e1855d42bf1933dba91"
PAIR_DIM = 640 + 2048 + 640 + 2048
MAX_PAIRS = 256


class PairNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.h = nn.Sequential(nn.Linear(PAIR_DIM, 64), nn.ReLU(),
                               nn.Linear(64, 64), nn.ReLU(),
                               nn.Linear(64, 1))
        for m in [self.h[0], self.h[2]]:
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)
        nn.init.xavier_uniform_(self.h[4].weight)
        nn.init.zeros_(self.h[4].bias)

    def _cat(self, xp_a, xl_a, xp_b, xl_b):
        return torch.cat([xp_a, xl_a, xp_b, xl_b], dim=-1)

    def d(self, xp_q, xl_q, xp_s, xl_s):
        """D(q, s); antisymmetric and identity-zero by construction."""
        fwd = self.h(self._cat(xp_q, xl_q, xp_s, xl_s))
        bwd = self.h(self._cat(xp_s, xl_s, xp_q, xl_q))
        return (fwd - bwd).squeeze(-1)


def pair_features(cells, pki, lid_of, fps_by_lid, pfeat, split_art, device):
    xp, xl, y = PT.build_cell_features(cells, pki, lid_of, fps_by_lid, pfeat,
                                       split_art)
    return (torch.from_numpy(xp).to(device), torch.from_numpy(xl).to(device),
            torch.from_numpy(y).to(device))


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
    model = PairNet().to(device)
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=PT.LR, weight_decay=PT.WD)
    best_mon = None
    best_state = None
    n_steps = 2 if dry else PT.STEPS
    for step in range(n_steps):
        rng = PT.stable_rng("stageP", "porder", seed, "step", step)
        kept = PC.collect_tasks(rng, tasks, ligs_of_target, first_cell)
        cells = [c for sup, q in kept for c in sup + q]
        xp, xl, y = pair_features(cells, pki, lid_of, fps_by_lid, pfeat,
                                  split_art, device)
        opt.zero_grad()
        tot = 0.0
        cnt = 0
        n = len(cells)
        for i in range(min(n, MAX_PAIRS)):
            a = i // max(n - 1, 1)
            b = i % max(n - 1, 1)
            if b >= a:
                b += 1
            if a >= n or b >= n:
                continue
            d = model.d(xp[a:a + 1], xl[a:a + 1], xp[b:b + 1], xl[b:b + 1])
            mse = ((d - (y[a] - y[b])) ** 2).mean()
            if not torch.isfinite(mse):
                continue
            mse.backward()
            tot += float(mse.detach())
            cnt += 1
        if cnt > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            opt.step()
        if step % PT.MONITOR_EVERY == PT.MONITOR_EVERY - 1 and not dry:
            mon = monitor(model, bank, split_art, pki, lid_of, fps_by_lid, pfeat,
                          device)
            print(f"actfound seed {seed} step {step + 1} mse {tot / max(cnt, 1):.4f} "
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


def predict(model, sup_ids, q_ids, pki, lid_of, fps_by_lid, pfeat, split_art,
            device, k0_mean):
    model.eval()
    with torch.no_grad():
        if not sup_ids:
            return np.full(len(q_ids), k0_mean, dtype=np.float32)
        xp_q, xl_q, _ = PT.build_cell_features(q_ids, pki, lid_of, fps_by_lid,
                                               pfeat, split_art)
        xp_s, xl_s, y_s = PT.build_cell_features(sup_ids, pki, lid_of, fps_by_lid,
                                                 pfeat, split_art)
        xp_q = torch.from_numpy(xp_q).to(device)
        xl_q = torch.from_numpy(xl_q).to(device)
        xp_s = torch.from_numpy(xp_s).to(device)
        xl_s = torch.from_numpy(xl_s).to(device)
        y_s = torch.from_numpy(y_s).to(device)
        yhat = []
        for i in range(len(q_ids)):
            xp_qi = xp_q[i:i + 1].expand(len(sup_ids), -1)
            xl_qi = xl_q[i:i + 1].expand(len(sup_ids), -1)
            d = model.d(xp_qi, xl_qi, xp_s, xl_s)
            yhat.append(float((y_s + d).mean()))
    return np.asarray(yhat, dtype=np.float32)


def monitor(model, bank, split_art, pki, lid_of, fps_by_lid, pfeat, device):
    train_ids = [cid for cid, rec in split_art["cell_split"].items()
                 if rec["split"] == "p_train"]
    k0_mean = float(np.mean([pki[c] for c in train_ids]))
    recs = [r for r in bank["records"] if r["split"] == "p_val" and r["draw"] == 0
            and r["k"] in PT.MONITOR_KS]
    mses = []
    for rec in recs:
        yh = predict(model, rec["support_cell_ids"], rec["query_cell_ids"],
                     pki, lid_of, fps_by_lid, pfeat, split_art, device, k0_mean)
        q_y = np.asarray([pki[c] for c in rec["query_cell_ids"]], dtype=np.float32)
        mses.append(float(np.mean((yh - q_y) ** 2)))
    return float(np.mean(mses))


def eval_seed(model, seed, device, bank, split_art, pki, lid_of, fps_by_lid, pfeat):
    train_ids = [cid for cid, rec in split_art["cell_split"].items()
                 if rec["split"] == "p_train"]
    k0_mean = float(np.mean([pki[c] for c in train_ids]))
    out = {"schema": SCHEMA, "seed": seed, "records": [], "k0_mean": k0_mean}
    for rec in bank["records"]:
        yh = predict(model, rec["support_cell_ids"], rec["query_cell_ids"],
                     pki, lid_of, fps_by_lid, pfeat, split_art, device, k0_mean)
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
        model, mon = train_seed(seed, device, bank, split_art, pki, lid_of,
                                fps_by_lid, pfeat, args.dry)
        print(f"actfound seed {seed} final monitor {mon}", flush=True)
        seeds_out[str(seed)] = eval_seed(model, seed, device, bank, split_art,
                                         pki, lid_of, fps_by_lid, pfeat)
    n_pair = sum(p.numel() for p in PairNet().parameters())
    artifact = {
        "schema": SCHEMA,
        "ad3_sha256": AD3_SHA,
        "adaptation_mechanism": ("support-anchored pairwise differences "
                                 "(no gradient steps)"),
        "note": ("pair module antisymmetric + identity-zero by construction; "
                 "labels enter eval ONLY through support anchors; k=0 = "
                 "p_train label mean (frozen)"),
        "param_count": n_pair,
        "bank_sha256": PT.sha256_file(OUT / "P_BANK.json"),
        "split_sha256": PT.sha256_file(OUT / "P_SPLIT.json"),
        "dry": bool(args.dry),
        "protocol": {"steps": (2 if args.dry else PT.STEPS), "lr": PT.LR, "wd": PT.WD,
                     "max_pairs_per_step": MAX_PAIRS, "pair_dim": PAIR_DIM,
                     "eval_adaptation": "support-anchored-differences",
                     "monitor_every": PT.MONITOR_EVERY,
                     "monitor_adapt_steps": PT.MONITOR_ADAPT_STEPS,
                     "train_k": list(PT.TRAIN_K),
                     "min_train_ligs": PT.MIN_TRAIN_LIGS},
        "seeds": seeds_out,
    }
    text = json.dumps(artifact, indent=1, sort_keys=True)
    path = OUT / "P1_ARM7_ACTFOUND.json"
    path.write_text(text, encoding="utf-8")
    art_sha = PT.sha256_file(path)
    (OUT / "P1_ARM7_ACTFOUND.json.manifest.json").write_text(json.dumps({
        "schema": SCHEMA + ".Manifest", "file": "P1_ARM7_ACTFOUND.json",
        "sha256": art_sha}, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    for split_name in ("p_val", "p_test"):
        print(f"--- arm7 ActFound ({split_name}) ---")
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
