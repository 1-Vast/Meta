"""XP1-A: does a reproducible protein-by-ligand interaction term exist?

Registered in PREREG_XP1.md sections 2, 7, 8.  No transfer arm is run here.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panels import (additive_fit, load_klaeger, load_klifs, load_metz,  # noqa: E402
                    load_pdsp, map_kinases, verify_releases)
from lowrank import cv_rank_curve  # noqa: E402

REPORT = r"D:\MetaSieve\report\crossed_panel_identification"
os.makedirs(REPORT, exist_ok=True)
OUT = {}


def var_decomposition(Y, M, label):
    mu, a, b, fit = additive_fit(Y, M)
    y = Y[M]
    g = (Y - fit)[M]
    tot = float(((y - y.mean()) ** 2).sum())
    # sequential sums of squares on observed cells
    ss_a = float((((a[:, None] + np.zeros_like(Y))[M]) ** 2).sum())
    ss_b = float((((b[None, :] + np.zeros_like(Y))[M]) ** 2).sum())
    ss_g = float((g ** 2).sum())
    res = {
        "cells": int(M.sum()),
        "n_ligand": int(Y.shape[0]),
        "n_protein": int(Y.shape[1]),
        "sd_total": float(y.std(ddof=1)),
        "sd_alpha": float(a.std(ddof=1)),
        "sd_beta": float(b.std(ddof=1)),
        "sd_gamma_plus_noise": float(g.std(ddof=1)),
        "share_alpha": ss_a / tot,
        "share_beta": ss_b / tot,
        "share_gamma_plus_noise": ss_g / tot,
    }
    print(f"\n[{label}] two-way decomposition on {res['cells']} observed cells")
    print(f"   sd(y)={res['sd_total']:.3f}  sd(alpha)={res['sd_alpha']:.3f} "
          f"sd(beta)={res['sd_beta']:.3f}  sd(gamma+eps)={res['sd_gamma_plus_noise']:.3f}")
    print(f"   variance share  ligand={res['share_alpha']:.3f}  "
          f"protein={res['share_beta']:.3f}  interaction+noise={res['share_gamma_plus_noise']:.3f}")
    return res


def interaction_geometry(Y, M, min_overlap=25):
    """Kinase x kinase Pearson correlation of interaction residual columns."""
    mu, a, b, fit = additive_fit(Y, M)
    G = np.where(M, Y - fit, np.nan)
    p = Y.shape[1]
    C = np.full((p, p), np.nan)
    for j in range(p):
        gj = G[:, j]
        for k in range(j + 1, p):
            gk = G[:, k]
            ok = np.isfinite(gj) & np.isfinite(gk)
            if ok.sum() >= min_overlap:
                x, y = gj[ok], gk[ok]
                if x.std() > 1e-9 and y.std() > 1e-9:
                    C[j, k] = C[k, j] = float(np.corrcoef(x, y)[0, 1])
    return C


def compare_geometry(C1, C2, labels, rng, n_perm=2000):
    p = C1.shape[0]
    iu = np.triu_indices(p, 1)
    x, y = C1[iu], C2[iu]
    ok = np.isfinite(x) & np.isfinite(y)
    r = float(np.corrcoef(x[ok], y[ok])[0, 1])
    rho = float(stats.spearmanr(x[ok], y[ok]).statistic)
    null = []
    for _ in range(n_perm):
        q = rng.permutation(p)
        C2p = C2[np.ix_(q, q)]
        yp = C2p[iu]
        o2 = np.isfinite(x) & np.isfinite(yp)
        if o2.sum() > 10:
            null.append(abs(float(np.corrcoef(x[o2], yp[o2])[0, 1])))
    null = np.array(null)
    pval = float((null >= abs(r)).mean())
    return {"pairs": int(ok.sum()), "pearson": r, "spearman": rho,
            "perm_null_mean": float(null.mean()), "perm_null_p95": float(np.percentile(null, 95)),
            "perm_p": pval, "labels": labels}


# ==========================================================================
print("=" * 78)
print("XP1-A  release verification")
OUT["release_verification"] = verify_releases()
for k, v in OUT["release_verification"].items():
    print(f"  {k:22s} sha256 match = {v['match']}")

rng = np.random.default_rng(20260808)

# --------------------------------------------------------------- A1 / A2
print("\n" + "=" * 78)
print("XP1-A1  variance decomposition + XP1-A2 in-matrix rank curve")
for dens, tag in [(0.60, "BLK-METZ-60"), (0.70, "BLK-METZ-70")]:
    Y, M, cid, kin = load_metz(dens)
    OUT[f"vardecomp_{tag}"] = var_decomposition(Y, M, tag)
    OUT[f"vardecomp_{tag}"]["kinases"] = list(map(str, kin))

Y, M, cid, kin = load_metz(0.60)
ranks = [1, 2, 3, 5, 8, 12, 20, 30]
print("\n  random-cell 5-fold CV, additive vs additive+rank r")
curve = cv_rank_curve(Y, M, ranks, lam=2.0, folds=5, seed=0)
print(f"   additive RMSE = {curve['additive']['rmse']:.4f} log units")
rows = []
for r in ranks:
    c = curve[r]
    d = np.array(curve["additive"]["per_fold_mse"]) - np.array(c["per_fold_mse"])
    lo, hi = (d.mean() - 2.776 * d.std(ddof=1) / np.sqrt(len(d)),
              d.mean() + 2.776 * d.std(ddof=1) / np.sqrt(len(d)))
    print(f"   rank {r:3d}: RMSE={c['rmse']:.4f}  R2_gamma={c['r2_gamma']:+.4f}   "
          f"per-fold MSE gain {d.mean():+.5f} [{lo:+.5f},{hi:+.5f}]")
    rows.append({"rank": r, "rmse": c["rmse"], "r2_gamma": c["r2_gamma"],
                 "mse_gain": float(d.mean()), "ci95": [float(lo), float(hi)]})
OUT["rank_curve_BLK-METZ-60"] = {"additive_rmse": curve["additive"]["rmse"], "rows": rows}

# --------------------------------------------------------------- A3
print("\n" + "=" * 78)
print("XP1-A3  protein-side interaction geometry: disjoint compound halves (Metz)")
n = Y.shape[0]
perm = rng.permutation(n)
h1, h2 = perm[: n // 2], perm[n // 2:]
C1 = interaction_geometry(Y[h1], M[h1])
C2 = interaction_geometry(Y[h2], M[h2])
res = compare_geometry(C1, C2, ("metz_half1", "metz_half2"), rng)
print(f"   kinase pairs compared: {res['pairs']}")
print(f"   Pearson r = {res['pearson']:.4f}   Spearman = {res['spearman']:.4f}")
print(f"   label-permutation null |r|: mean {res['perm_null_mean']:.4f}, "
      f"p95 {res['perm_null_p95']:.4f}, p = {res['perm_p']:.4g}")
OUT["geometry_split_half_metz"] = res

# main-effect contrast: same test on the *raw* columns (contains beta)
def raw_geometry(Y, M, min_overlap=25):
    p = Y.shape[1]
    C = np.full((p, p), np.nan)
    Yn = np.where(M, Y, np.nan)
    for j in range(p):
        for k in range(j + 1, p):
            ok = np.isfinite(Yn[:, j]) & np.isfinite(Yn[:, k])
            if ok.sum() >= min_overlap:
                x, y = Yn[ok, j], Yn[ok, k]
                if x.std() > 1e-9 and y.std() > 1e-9:
                    C[j, k] = C[k, j] = float(np.corrcoef(x, y)[0, 1])
    return C


R1 = raw_geometry(Y[h1], M[h1])
R2 = raw_geometry(Y[h2], M[h2])
res_raw = compare_geometry(R1, R2, ("metz_half1_raw", "metz_half2_raw"), rng)
print(f"   [raw, main effects retained] Pearson r = {res_raw['pearson']:.4f}  "
      f"perm p = {res_raw['perm_p']:.4g}")
OUT["geometry_split_half_metz_raw"] = res_raw

# --------------------------------------------------------------- A4
print("\n" + "=" * 78)
print("XP1-A4  cross-platform replication: Metz (pKi) vs Klaeger (kinobeads)")
W, hit, drug, kinK = load_klaeger()
keep_k = hit.sum(axis=0) >= 5
keep_d = hit.sum(axis=1) >= 5
Wb = hit[np.ix_(keep_d, keep_k)].astype(float)
kinK2 = kinK[keep_k]
shared = [s for s in kin if s in set(kinK2)]
print(f"   shared kinases usable in both panels: {len(shared)}")
if len(shared) >= 20:
    im = {s: i for i, s in enumerate(kin)}
    ik = {s: i for i, s in enumerate(kinK2)}
    jm = [im[s] for s in shared]
    jk = [ik[s] for s in shared]
    Cm = interaction_geometry(Y[:, jm], M[:, jm], min_overlap=25)
    Mk = np.ones_like(Wb, bool)
    Ck = interaction_geometry(Wb[:, jk], Mk[:, jk], min_overlap=25)
    res_x = compare_geometry(Cm, Ck, ("metz_gamma", "klaeger_gamma"), rng)
    print(f"   kinase pairs compared: {res_x['pairs']}")
    print(f"   Pearson r = {res_x['pearson']:.4f}   Spearman = {res_x['spearman']:.4f}")
    print(f"   label-permutation null |r|: mean {res_x['perm_null_mean']:.4f}, "
          f"p95 {res_x['perm_null_p95']:.4f}, p = {res_x['perm_p']:.4g}")
    OUT["geometry_cross_platform"] = res_x
    OUT["geometry_cross_platform"]["shared_kinases"] = shared

# --------------------------------------------------------------- A5
print("\n" + "=" * 78)
print("XP1-A5  measurement-noise ceiling from PDSP independent replicates")
p = load_pdsp()
g = p.groupby(["target", "ligand"])["pKi"]
sizes = g.size()
rep = sizes[sizes >= 2]
print(f"   human uncensored rows {len(p)}, cells {len(sizes)}, replicated cells {len(rep)}")
sub = p.set_index(["target", "ligand"]).loc[rep.index].reset_index()
halves = []
rs = np.random.default_rng(7)
for (t, l), d in sub.groupby(["target", "ligand"]):
    v = d["pKi"].to_numpy()
    rs.shuffle(v)
    k = len(v) // 2
    if k >= 1 and len(v) - k >= 1:
        halves.append((t, l, v[:k].mean(), v[k:].mean()))
H = pd.DataFrame(halves, columns=["target", "ligand", "h1", "h2"])
d = H["h1"] - H["h2"]
sd_pair = float(d.std(ddof=1))
sigma_rep = sd_pair / np.sqrt(2.0)
print(f"   replicate half-split pairs: {len(H)}")
print(f"   sd(h1-h2) = {sd_pair:.4f} log units -> per-report sigma ~ {sigma_rep:.4f}")
print(f"   r(h1,h2) = {np.corrcoef(H['h1'], H['h2'])[0,1]:.4f}")
OUT["pdsp_replicate_noise"] = {
    "replicated_cells": int(len(H)),
    "sd_half_difference": sd_pair,
    "sigma_per_report": float(sigma_rep),
    "r_half_half": float(np.corrcoef(H["h1"], H["h2"])[0, 1]),
}

# KLIFS coverage of the Metz block
kl = load_klifs()
hitk, missk = map_kinases(kin, kl)
print(f"\n   KLIFS pocket coverage of BLK-METZ-60: {len(hitk)}/{len(kin)}; missing={missk}")
OUT["klifs_coverage_metz60"] = {"mapped": len(hitk), "total": int(len(kin)), "missing": missk}

with open(os.path.join(REPORT, "xp1a_existence.json"), "w") as f:
    json.dump(OUT, f, indent=2, default=float)
print("\nwrote", os.path.join(REPORT, "xp1a_existence.json"))
