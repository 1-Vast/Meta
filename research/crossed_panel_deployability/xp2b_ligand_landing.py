"""XP2-B — can the ligand loading u(L) be predicted from chemistry?

Scaffold-component folds.  For each fold the interaction basis V and the main
effects are fitted on TRAINING ligands only; the target parameters of a held-out
ligand are its own least-squares (alpha, u) against that frozen basis; the
predictor chi sees chemistry only.  No ligand ID, no lookup, no test-scaffold
information anywhere in fitting or model selection.

Reported gauge-invariantly wherever the quantity is basis-dependent.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from xp2_core import (cluster_bootstrap, interaction_basis, ridge_cv,  # noqa: E402
                      solve_ligand_params)
from xp2_panel import build, ligand_features  # noqa: E402

OUT = r"D:\MetaSieve\report\crossed_panel_deployability"
os.makedirs(OUT, exist_ok=True)
RANKS = (1, 2, 3, 5)
PRIMARY_RANK = 3
LAM_LR = 2.0
LAM_LIG = 1.0
N_FOLDS = 5


def run(rank=PRIMARY_RANK, seed=0):
    d, man = build()
    Y, M = d["Y"], d["M"]
    smiles = d["smiles"]
    comp = d["scaffold_component"]
    FEAT = ligand_features(smiles)
    comps = sorted(set(comp))
    rng = np.random.default_rng(seed)
    order = list(comps)
    rng.shuffle(order)
    fold_of = {c: i % N_FOLDS for i, c in enumerate(order)}

    arms = list(FEAT) + ["L-MEAN"]
    # per scaffold-component accumulators
    acc = {a: {} for a in arms}          # gauge-invariant interaction reconstruction
    acc_alpha = {a: {} for a in arms}    # ligand main effect prediction
    tgt = {}                             # target interaction energy per component
    coord_r2 = {a: [] for a in arms}
    sign_rate = {a: [] for a in arms}
    spear = {a: [] for a in arms}
    lam_used = {a: [] for a in arms}

    for f in range(N_FOLDS):
        te_c = [c for c in comps if fold_of[c] == f]
        te = np.isin(comp, te_c)
        tr = ~te
        m = interaction_basis(Y[tr], M[tr], rank, lam=LAM_LR, seed=100 + f)
        mu, beta, V = m["mu"], m["beta"], m["V"]

        # target parameters for every ligand against the frozen training basis
        A_star = np.zeros(len(Y))
        U_star = np.zeros((len(Y), rank))
        for i in range(len(Y)):
            if M[i].sum() >= rank + 2:
                A_star[i], U_star[i] = solve_ligand_params(Y[i], M[i], mu, beta, V, LAM_LIG)

        usable = M.sum(1) >= rank + 2
        tr_u = tr & usable
        te_u = te & usable
        if te_u.sum() < 3:
            continue
        inner = comp[tr_u]

        for a in arms:
            if a == "L-MEAN":
                Uhat = np.repeat(U_star[tr_u].mean(0)[None, :], te_u.sum(), axis=0)
                Ahat = np.repeat(A_star[tr_u].mean(), te_u.sum())
                lam = None
            else:
                X = FEAT[a]
                Uhat, lam = ridge_cv(X[tr_u], U_star[tr_u], X[te_u], inner)
                Ahat, _ = ridge_cv(X[tr_u], A_star[tr_u][:, None], X[te_u], inner)
                Ahat = Ahat[:, 0]
            lam_used[a].append(lam)

            # coordinate-wise R2 in the fold-local gauge (basis-dependent, descriptive)
            ss = ((U_star[te_u] - Uhat) ** 2).sum()
            st = ((U_star[te_u] - U_star[tr_u].mean(0)) ** 2).sum()
            coord_r2[a].append(float(1 - ss / max(st, 1e-12)))

            # gauge-invariant: reconstruct the interaction energy on measured cells
            te_idx = np.where(te_u)[0]
            for t, i in enumerate(te_idx):
                j = np.where(M[i])[0]
                true_g = U_star[i] @ V[j].T
                pred_g = Uhat[t] @ V[j].T
                c = comp[i]
                s, n = acc[a].get(c, (0.0, 0))
                acc[a][c] = (s + float(((true_g - pred_g) ** 2).sum()), n + len(j))
                if a == arms[0]:
                    s2, n2 = tgt.get(c, (0.0, 0))
                    tgt[c] = (s2 + float((true_g ** 2).sum()), n2 + len(j))
                if true_g.std() > 1e-9 and pred_g.std() > 1e-9:
                    spear[a].append(float(stats.spearmanr(true_g, pred_g).statistic))
                    sign_rate[a].append(float(np.mean(np.sign(true_g) == np.sign(pred_g))))
                sA, nA = acc_alpha[a].get(c, (0.0, 0))
                acc_alpha[a][c] = (sA + float((A_star[i] - Ahat[t]) ** 2), nA + 1)

    res = {"rank": rank, "n_folds": N_FOLDS, "panel": man["panel"],
           "n_scaffold_components": man["n_scaffold_components"], "arms": {}}
    tgt_ms = sum(v[0] for v in tgt.values()) / sum(v[1] for v in tgt.values())
    res["target_interaction_mean_square"] = float(tgt_ms)
    alpha_var = None
    for a in arms:
        bs = cluster_bootstrap(acc[a], seed=11)
        r2 = 1.0 - bs["mse"] / tgt_ms
        # bootstrap the gauge-invariant R2 directly
        u = sorted(acc[a])
        sa = np.array([acc[a][k][0] for k in u])
        sb = np.array([tgt[k][0] for k in u])
        rr = np.random.default_rng(12)
        idx = rr.integers(0, len(u), (2000, len(u)))
        bsr = 1.0 - sa[idx].sum(1) / np.maximum(sb[idx].sum(1), 1e-12)
        ab = cluster_bootstrap(acc_alpha[a], seed=13)
        if alpha_var is None:
            alpha_var = ab["mse"]
        res["arms"][a] = {
            "interaction_reconstruction_r2_gauge_invariant": {
                "point": float(r2),
                "ci95": [float(np.percentile(bsr, 2.5)), float(np.percentile(bsr, 97.5))]},
            "loading_coord_r2_fold_local_gauge_mean": float(np.mean(coord_r2[a])),
            "loading_coord_r2_per_fold": [round(x, 4) for x in coord_r2[a]],
            "sign_agreement": float(np.mean(sign_rate[a])) if sign_rate[a] else None,
            "spearman_within_ligand": float(np.mean(spear[a])) if spear[a] else None,
            "alpha_prediction_mse": ab["mse"],
            "lambda_selected": [None if x is None else float(x) for x in lam_used[a]],
            "cells": bs["n"],
        }
    for a in arms:
        res["arms"][a]["alpha_prediction_r2_vs_mean_arm"] = float(
            1.0 - res["arms"][a]["alpha_prediction_mse"]
            / res["arms"]["L-MEAN"]["alpha_prediction_mse"])
    return res


if __name__ == "__main__":
    out = {"stage": "XP2-B", "prereg": "PREREG_XP2.md", "ranks": {}}
    for r in RANKS:
        print(f"\n{'='*78}\nrank d = {r}")
        res = run(rank=r)
        out["ranks"][str(r)] = res
        print(f"{'arm':16s} {'interaction recon R2 (gauge-inv)':>34s} "
              f"{'coordR2':>9s} {'sign':>6s} {'rho':>6s} {'alphaR2':>8s}")
        for a, v in res["arms"].items():
            g = v["interaction_reconstruction_r2_gauge_invariant"]
            print(f"{a:16s} {g['point']:+.4f} [{g['ci95'][0]:+.4f},{g['ci95'][1]:+.4f}]"
                  f"{'':>7s} {v['loading_coord_r2_fold_local_gauge_mean']:+9.4f} "
                  f"{(v['sign_agreement'] or 0):6.3f} {(v['spearman_within_ligand'] or 0):6.3f} "
                  f"{v['alpha_prediction_r2_vs_mean_arm']:+8.4f}")
    p = os.path.join(OUT, "LIGAND_LANDING_AUDIT.json")
    json.dump(out, open(p, "w"), indent=2, default=float)
    print("\nwrote", p)
