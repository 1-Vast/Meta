"""Stage P1 arm 3: ordinary support fine-tuning (frozen backbone spec
2d733684...). Shared episodic sampler keyed per (seed, step); monitor on
p_val draw-0 records (frozen checkpoint rule); final eval on p_val+p_test
bank records with the frozen test-time protocol (50 support steps, best
SUPPORT loss; query labels never enter adaptation or selection).

Usage: python p_train.py --seeds 1 2 3 [--device cuda] [--dry]
Artifacts: P1_ARM3_ORDINARYFT.json + manifest (SHA-pinned).
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

ROOT = Path(__file__).resolve().parents[3]
CORPUS = ROOT / "dataset" / "processed" / "meta_fewshot" / "bindingdb_ki_main_v0"
PBANK = ROOT / "dataset" / "processed" / "meta_fewshot" / "bindingdb_ki_main_v0_protein_bank"
HERE = Path(__file__).resolve().parent
OUT = HERE / "artifacts"
SCHEMA = "MetaSieve.StageP.P1Arm3OrdinaryFT.v1"
BACKBONE_SHA = "2d73368490f77e6e590f7631c7313ff42bacafb5ad4e9ec79f75cca20297b799"
K_LIST = (0, 1, 2, 3, 5, 10, 20, 40)
TRAIN_K = (5, 10, 20)
QUERY = 8
BATCH_CELLS = 256
STEPS = 6000
LR = 3e-4
WD = 1e-4
ADAPT_STEPS = 50
ADAPT_LR = 1e-3
MONITOR_EVERY = 600
MONITOR_ADAPT_STEPS = 10
MONITOR_KS = (0, 1, 2, 3, 5, 10)
MIN_TRAIN_LIGS = 28


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def stable_rng(*parts):
    raw = "|".join(str(x) for x in parts).encode("utf-8")
    return np.random.default_rng(int(hashlib.sha256(raw).hexdigest()[:16], 16))


class PTrunk(nn.Module):
    """Frozen P1 backbone (spec P1_BACKBONE_SPEC.md)."""

    def __init__(self, d_p=640, d_l=2048, hid=64, rank=16):
        super().__init__()
        self.p_enc = nn.Linear(d_p, hid)
        self.l_enc = nn.Linear(d_l, hid)
        self.p_head = nn.Linear(hid, 1, bias=False)
        self.l_head = nn.Linear(hid, 1, bias=False)
        self.mu = nn.Parameter(torch.zeros(1))
        self.A = nn.Linear(hid, rank, bias=False)
        self.B = nn.Linear(hid, rank, bias=False)
        self.inter_scale = nn.Parameter(torch.ones(1))
        self.inter_bias = nn.Parameter(torch.zeros(1))
        for name, m in self.named_modules():
            if isinstance(m, nn.Linear) and name not in ("p_head", "l_head"):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(self, xp, xl):
        ep = torch.relu(self.p_enc(xp))
        el = torch.relu(self.l_enc(xl))
        inter = self.inter_scale * ((self.A(ep) * self.B(el)).sum(-1)
                                    + self.inter_bias)
        yhat = (self.mu + self.p_head(ep).squeeze(-1)
                + self.l_head(el).squeeze(-1) + inter)
        return {"yhat": yhat, "inter": inter}


def load_protein_features():
    feat = {}
    for shard in sorted(PBANK.glob("shard_*.npz")):
        z = np.load(shard, allow_pickle=False)
        for i, key in enumerate(z["keys"]):
            feat[str(key)] = z["pooled"][i].astype(np.float32)
    return feat


def load_labels_and_ligands():
    import gzip
    pki = {}
    lid_of = {}
    with gzip.open(CORPUS / "cells.jsonl.gz", "rt", encoding="utf-8") as f:
        for line in f:
            c = json.loads(line)
            pki[c["cell_id"]] = float(c["pK"])
            lid_of[c["cell_id"]] = c["ligand_id"]
    lig = {}
    with open(CORPUS / "ligands.jsonl", encoding="utf-8") as f:
        for line in f:
            d = json.loads(line)
            lig[d["drug_key"]] = d["smiles"]
    return pki, lid_of, lig


def ecfp_matrix(lig_ids, lig):
    fps = []
    for lid in lig_ids:
        mol = Chem.MolFromSmiles(lig[lid])
        if mol is None:
            fps.append(np.zeros(2048, dtype=np.float16))
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            v = np.zeros(2048, dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, v)
            fps.append(v.astype(np.float16))
    return np.stack(fps)


def load_ecfp_cache(lig):
    cache = OUT / "ecfp4_2048_cache.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        return {str(k): v.astype(np.float32) for k, v in zip(z["keys"], z["fps"])}
    lig_ids = sorted(lig)
    fps = ecfp_matrix(lig_ids, lig)
    np.savez(cache, keys=np.asarray(lig_ids), fps=fps)
    return {lid: v.astype(np.float32) for lid, v in zip(lig_ids, fps)}


def build_cell_features(cell_ids, pki, lid_of, fps_by_lid, pfeat, split_art):
    xp, xl, y = [], [], []
    for c in cell_ids:
        tid = split_art["cell_split"][c]["target_id"]
        xp.append(pfeat[tid])
        y.append(pki[c])
        xl.append(fps_by_lid[lid_of[c]])
    return (np.stack(xp).astype(np.float32), np.stack(xl).astype(np.float32),
            np.asarray(y, dtype=np.float32))


def spearman(a, b):
    from scipy.stats import spearmanr
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(spearmanr(a, b).correlation)


def pearson(a, b):
    from scipy.stats import pearsonr
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(pearsonr(a, b)[0])


def ci(yhat, y):
    if len(y) < 2:
        return float("nan")
    conc = disc = 0
    for i in range(len(y)):
        for j in range(i + 1, len(y)):
            if y[i] == y[j]:
                continue
            if (yhat[i] > yhat[j]) == (y[i] > y[j]):
                conc += 1
            else:
                disc += 1
    if conc + disc == 0:
        return float("nan")
    return conc / (conc + disc)


def adapt_and_predict(model, sup_ids, q_ids, pki, lid_of, fps_by_lid, pfeat,
                      split_art, device, steps=ADAPT_STEPS, lr=ADAPT_LR):
    """Frozen test-time protocol: adapt on support, best SUPPORT loss,
    predict query. Query labels never used."""
    m = copy.deepcopy(model)
    m.to(device)
    if len(sup_ids) == 0:
        m.eval()
        with torch.no_grad():
            xp, xl, _ = build_cell_features(q_ids, pki, lid_of, fps_by_lid, pfeat, split_art)
            out = m(torch.from_numpy(xp).to(device), torch.from_numpy(xl).to(device))
        return out["yhat"].cpu().numpy(), None
    xp_s, xl_s, y_s = build_cell_features(sup_ids, pki, lid_of, fps_by_lid, pfeat, split_art)
    xp_s = torch.from_numpy(xp_s).to(device)
    xl_s = torch.from_numpy(xl_s).to(device)
    y_s = torch.from_numpy(y_s).to(device)
    opt = torch.optim.AdamW([p for p in m.parameters() if p.requires_grad],
                            lr=lr, weight_decay=WD)
    m.train()
    best = None
    best_state = None
    for _ in range(steps):
        opt.zero_grad()
        loss = ((m(xp_s, xl_s)["yhat"] - y_s) ** 2).mean()
        loss.backward()
        opt.step()
        val = float(loss.detach())
        if best is None or val < best - 1e-9:
            best = val
            best_state = {k: v.detach().clone() for k, v in m.state_dict().items()}
    if best_state is not None:
        m.load_state_dict(best_state)
    m.eval()
    with torch.no_grad():
        xp, xl, _ = build_cell_features(q_ids, pki, lid_of, fps_by_lid, pfeat, split_art)
        out = m(torch.from_numpy(xp).to(device), torch.from_numpy(xl).to(device))
    return out["yhat"].cpu().numpy(), best


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
    tasks = sorted(t for t in ligs_of_target if len(ligs_of_target[t]) >= MIN_TRAIN_LIGS)
    first_cell = {t: {} for t in tasks}
    for t in tasks:
        for c in cells_of_target[t]:
            first_cell[t].setdefault(lid_of[c], c)
    model = PTrunk().to(device)
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=LR, weight_decay=WD)
    best_mon = None
    best_state = None
    n_steps = 2 if dry else STEPS
    for step in range(n_steps):
        rng = stable_rng("stageP", "porder", seed, "step", step)
        cells_batch = []
        while len(cells_batch) < BATCH_CELLS:
            t = tasks[int(rng.integers(len(tasks)))]
            k = TRAIN_K[int(rng.integers(len(TRAIN_K)))]
            ligs = sorted(ligs_of_target[t])
            perm = rng.permutation(len(ligs))
            ordered = [ligs[i] for i in perm[:k + QUERY]]
            cells_batch.extend([first_cell[t][l] for l in ordered])
        cells_batch = cells_batch[:BATCH_CELLS]
        xp, xl, y = build_cell_features(cells_batch, pki, lid_of, fps_by_lid, pfeat, split_art)
        opt.zero_grad()
        out = model(torch.from_numpy(xp).to(device), torch.from_numpy(xl).to(device))
        loss = ((out["yhat"] - torch.from_numpy(y).to(device)) ** 2).mean()
        loss.backward()
        opt.step()
        if step % MONITOR_EVERY == MONITOR_EVERY - 1 and not dry:
            mon = monitor_value(model, bank, split_art, pki, lid_of, fps_by_lid,
                                pfeat, device)
            print(f"seed {seed} step {step + 1} loss {float(loss):.4f} monitor {mon:.4f}",
                  flush=True)
            if best_mon is None or mon < best_mon - 1e-9:
                best_mon = mon
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if dry:
        mon = float(loss.detach())
        if best_state is None:
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            best_mon = mon
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, best_mon


def monitor_value(model, bank, split_art, pki, lid_of, fps_by_lid, pfeat, device):
    recs = [r for r in bank["records"] if r["split"] == "p_val" and r["draw"] == 0
            and r["k"] in MONITOR_KS]
    mses = []
    for rec in recs:
        yh, _ = adapt_and_predict(model, rec["support_cell_ids"],
                                  rec["query_cell_ids"], pki, lid_of, fps_by_lid,
                                  pfeat, split_art, device, steps=MONITOR_ADAPT_STEPS)
        q_y = np.asarray([pki[c] for c in rec["query_cell_ids"]], dtype=np.float32)
        mses.append(float(np.mean((yh - q_y) ** 2)))
    return float(np.mean(mses))


def eval_seed(model, seed, device, bank, split_art, pki, lid_of, fps_by_lid, pfeat,
              dry):
    out = {"schema": SCHEMA, "seed": seed, "records": []}
    for rec in bank["records"]:
        yh, best_sup = adapt_and_predict(
            model, rec["support_cell_ids"], rec["query_cell_ids"], pki, lid_of,
            fps_by_lid, pfeat, split_art, device,
            steps=(2 if dry else ADAPT_STEPS))
        q_y = np.asarray([pki[c] for c in rec["query_cell_ids"]], dtype=np.float32)
        out["records"].append({
            "split": rec["split"], "k": rec["k"], "draw": rec["draw"],
            "target_id": rec["target_id"], "cluster": rec["cluster"],
            "mse": float(np.mean((yh - q_y) ** 2)),
            "centered_mse": float(np.mean(((yh - q_y) - np.mean(yh - q_y)) ** 2)),
            "ci": ci(yh, q_y), "spearman": spearman(yh, q_y),
            "pearson": pearson(yh, q_y), "best_support_loss": best_sup,
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
    pki, lid_of, lig = load_labels_and_ligands()
    pfeat = load_protein_features()
    fps_by_lid = load_ecfp_cache(lig)
    seeds_out = {}
    for seed in args.seeds:
        model, mon = train_seed(seed, device, bank, split_art, pki, lid_of,
                                fps_by_lid, pfeat, args.dry)
        print(f"seed {seed} final monitor {mon}", flush=True)
        seeds_out[str(seed)] = eval_seed(model, seed, device, bank, split_art, pki,
                                         lid_of, fps_by_lid, pfeat, args.dry)
    artifact = {
        "schema": SCHEMA,
        "backbone_spec_sha256": BACKBONE_SHA,
        "bank_sha256": sha256_file(OUT / "P_BANK.json"),
        "split_sha256": sha256_file(OUT / "P_SPLIT.json"),
        "dry": bool(args.dry),
        "protocol": {"steps": (2 if args.dry else STEPS), "lr": LR, "wd": WD,
                     "batch_cells": BATCH_CELLS,
                     "adapt_steps": (2 if args.dry else ADAPT_STEPS),
                     "adapt_lr": ADAPT_LR,
                     "monitor_every": MONITOR_EVERY,
                     "monitor_adapt_steps": MONITOR_ADAPT_STEPS,
                     "train_k": list(TRAIN_K), "min_train_ligs": MIN_TRAIN_LIGS},
        "seeds": seeds_out,
    }
    text = json.dumps(artifact, indent=1, sort_keys=True)
    path = OUT / "P1_ARM3_ORDINARYFT.json"
    path.write_text(text, encoding="utf-8")
    art_sha = sha256_file(path)
    (OUT / "P1_ARM3_ORDINARYFT.json.manifest.json").write_text(json.dumps({
        "schema": SCHEMA + ".Manifest", "file": "P1_ARM3_ORDINARYFT.json",
        "sha256": art_sha}, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    for split_name in ("p_val", "p_test"):
        print(f"--- arm3 ordinary FT ({split_name}) ---")
        for k in K_LIST:
            rows = [r for s in seeds_out.values() for r in s["records"]
                    if r["split"] == split_name and r["k"] == k]
            if not rows:
                print(f"  k={k:<2} no records")
                continue
            print(f"  k={k:<2} n={len(rows)} MSE={np.mean([r['mse'] for r in rows]):.4f} "
                  f"cMSE={np.mean([r['centered_mse'] for r in rows]):.4f} "
                  f"CI={np.nanmean([r['ci'] for r in rows]):.4f} "
                  f"rho={np.nanmean([r['spearman'] for r in rows]):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
