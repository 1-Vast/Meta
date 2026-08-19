"""P1 cheap CPU baselines on the frozen episode bank (prereg 59a90ef2).

Arms implemented here (allowed baselines only, never final methods):
- ligand_only: per-target support-mean pKi; k=0 -> p_train global mean.
- tanimoto: ECFP4/2048 Tanimoto top-3 neighbour weighted mean of support
  pKi; k=0 -> p_train global mean.
Also reports exact-ligand-recall fractions (train-seen / support-seen)
required by the P1 promotion gate.

Metrics per record: MSE, RMSE, centered MSE, CI (within-query pairs),
Spearman, Pearson. Aggregation: record-mean per (split, k) and pooled
Spearman/Pearson/MSE over all query cells. SHA-pinned artifact out.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

ROOT = Path(__file__).resolve().parents[3]
CORPUS = ROOT / "dataset" / "processed" / "meta_fewshot" / "bindingdb_ki_main_v0"
HERE = Path(__file__).resolve().parent
OUT = HERE / "artifacts"
SCHEMA = "MetaSieve.StageP.P1Baselines.v1"

K_LIST = (0, 1, 2, 3, 5, 10, 20, 40)


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_labels():
    pki = {}
    train_cell_ids = []
    with gzip.open(CORPUS / "cells.jsonl.gz", "rt", encoding="utf-8") as f:
        for line in f:
            c = json.loads(line)
            pki[c["cell_id"]] = c["pK"]
    return pki


def load_ligands():
    lig = {}
    with open(CORPUS / "ligands.jsonl", encoding="utf-8") as f:
        for line in f:
            d = json.loads(line)
            lig[d["drug_key"]] = d["smiles"]
    return lig


def ecfp(smiles: str, radius: int = 2, nbits: int = 2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    return fp


def spearman(a, b):
    from scipy.stats import spearmanr, pearsonr
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


def main() -> int:
    bank = json.loads((OUT / "P_BANK.json").read_text(encoding="utf-8"))
    split_art = json.loads((OUT / "P_SPLIT.json").read_text(encoding="utf-8"))
    pki = load_labels()
    lig = load_ligands()
    # global p_train mean
    train_ids = [cid for cid, rec in split_art["cell_split"].items()
                 if rec["split"] == "p_train"]
    gmean = float(np.mean([pki[c] for c in train_ids]))
    train_ligs = {split_art["cell_split"][c]["ligand_id"] for c in train_ids}
    fp_cache = {}
    def get_fp(lid):
        if lid not in fp_cache:
            fp_cache[lid] = ecfp(lig[lid])
        return fp_cache[lid]

    out = {"schema": SCHEMA,
           "bank_sha256": sha256_file(OUT / "P_BANK.json"),
           "split_sha256": sha256_file(OUT / "P_SPLIT.json"),
           "arms": {"ligand_only": {}, "tanimoto": {}},
           "global_mean_pki": gmean,
           "n_train_cells": len(train_ids)}
    per_arm_records = {"ligand_only": {}, "tanimoto": {}}
    for arm in ("ligand_only", "tanimoto"):
        per_arm_records[arm] = {str(k): [] for k in K_LIST}
        per_arm_records[arm]["split"] = {}
    splits_of = {}
    for rec in bank["records"]:
        k = rec["k"]
        sup_y = np.array([pki[c] for c in rec["support_cell_ids"]], dtype=np.float64)
        q_ids = rec["query_cell_ids"]
        q_y = np.array([pki[c] for c in q_ids], dtype=np.float64)
        q_ligs = [split_art["cell_split"][c]["ligand_id"] for c in q_ids]
        if k == 0:
            yh_lo = np.full(len(q_ids), gmean)
            yh_t = np.full(len(q_ids), gmean)
        else:
            yh_lo = np.full(len(q_ids), sup_y.mean())
            yh_t = []
            sup_ligs = [split_art["cell_split"][c]["ligand_id"]
                        for c in rec["support_cell_ids"]]
            for ql in q_ligs:
                qf = get_fp(ql)
                sims = []
                for sl in sup_ligs:
                    sf = get_fp(sl)
                    if qf is None or sf is None:
                        sims.append(0.0)
                    else:
                        sims.append(DataStructs.TanimotoSimilarity(qf, sf))
                sims = np.asarray(sims)
                top = np.argsort(-sims)[:3]
                w = sims[top]
                if w.sum() <= 0:
                    yh_t.append(float(sup_y.mean()))
                else:
                    yh_t.append(float((w * sup_y[top]).sum() / w.sum()))
            yh_t = np.asarray(yh_t)
        for arm, yh in (("ligand_only", yh_lo), ("tanimoto", yh_t)):
            mse = float(np.mean((yh - q_y) ** 2))
            cmse = float(np.mean(((yh - q_y) - (np.mean(yh - q_y))) ** 2))
            rec_metrics = {
                "mse": mse,
                "rmse": float(math.sqrt(mse)),
                "centered_mse": cmse,
                "ci": ci(yh, q_y),
                "spearman": spearman(yh, q_y),
                "pearson": pearson(yh, q_y),
                "train_seen_frac": float(np.mean([l in train_ligs for l in q_ligs])),
                "support_seen_frac": float(np.mean(
                    [l in set(split_art["cell_split"][c]["ligand_id"]
                              for c in rec["support_cell_ids"]) for l in q_ligs])),
                "split": rec["split"],
            }
            per_arm_records[arm][str(k)].append(rec_metrics)
    for arm in ("ligand_only", "tanimoto"):
        for split_name in ("p_val", "p_test"):
            for k in K_LIST:
                rows = [r for r in per_arm_records[arm][str(k)]
                        if r["split"] == split_name]
                if not rows:
                    out["arms"][arm][f"{split_name}:k{k}"] = None
                    continue
                agg = {key: float(np.nanmean([r[key] for r in rows]))
                       for key in ("mse", "rmse", "centered_mse", "ci",
                                   "spearman", "pearson")}
                agg["n_records"] = len(rows)
                agg["train_seen_frac"] = float(np.mean(
                    [r["train_seen_frac"] for r in rows]))
                agg["support_seen_frac"] = float(np.mean(
                    [r["support_seen_frac"] for r in rows]))
                out["arms"][arm][f"{split_name}:k{k}"] = agg
    text = json.dumps(out, indent=1, sort_keys=True)
    path = OUT / "P1_BASELINES.json"
    path.write_text(text, encoding="utf-8")
    art_sha = sha256_file(path)
    (OUT / "P1_BASELINES.manifest.json").write_text(json.dumps({
        "schema": SCHEMA + ".Manifest",
        "file": "P1_BASELINES.json",
        "sha256": art_sha,
    }, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    for arm in ("ligand_only", "tanimoto"):
        print(f"--- {arm} (p_test) ---")
        for k in K_LIST:
            key = f"p_test:k{k}"
            a = out["arms"][arm].get(key)
            if a is None:
                print(f"  k={k:<2} no records")
                continue
            print(f"  k={k:<2} n={a['n_records']:<4} MSE={a['mse']:.4f} "
                  f"RMSE={a['rmse']:.4f} cMSE={a['centered_mse']:.4f} "
                  f"CI={a['ci']:.4f} rho={a['spearman']:.4f} "
                  f"trainSeen={a['train_seen_frac']:.3f} "
                  f"supSeen={a['support_seen_frac']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
