"""XP2-C / XP2-D — true few-shot section audit and double-held-out deployability.

closure='protein'  : hold out protein groups only, ligands reused (XP1-comparable)
closure='double'   : hold out protein groups AND ligand scaffold components

The query ligand's main effect and loading always come from CHEMISTRY, never from
a fitted per-ligand table, so no ligand-ID lookup can carry the result.  Support
ligands are always drawn from TRAINING scaffolds of the held-out protein.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from xp2_core import (interaction_basis, paired_contrast, r2_vs_base,  # noqa: E402
                      ridge_cv, widest)
from xp2_panel import build, ligand_features  # noqa: E402

OUT = r"D:\MetaSieve\report\crossed_panel_deployability"
os.makedirs(OUT, exist_ok=True)

LAM_LR, LAM_V = 2.0, 1.0
N_FOLDS = 5
SEEDS = (0, 1, 2, 3, 4)
ARMS = ["P0", "LIG", "ADD", "ZERO", "SEC", "SEC-UHATSUP", "PERM", "FOREIGN",
        "RANDCORR", "ORACLE-TRSC"]


def solve_bv(resid, U, idx, lam):
    D = np.column_stack([np.ones(len(idx)), U[idx]])
    pen = np.diag([0.0] + [lam] * U.shape[1])
    return np.linalg.solve(D.T @ D + pen, D.T @ resid[idx])


def run(closure="double", rank=3, k=5, ligand_arm="L-ECFP", n_boot=2000, verbose=True):
    d, man = build()
    Y, M = d["Y"], d["M"]
    grp, scomp = d["group"], d["scaffold_component"]
    X = ligand_features(d["smiles"], kinds=("desc", "ecfp", "chemberta", "random"))[ligand_arm]

    pcomps, scomps = sorted(set(grp)), sorted(set(scomp))
    po = sorted(pcomps, key=lambda c: -int((grp == c).sum()))
    pfold = {c: i % N_FOLDS for i, c in enumerate(po)}
    so = list(scomps)
    np.random.default_rng(0).shuffle(so)
    sfold = {c: i % N_FOLDS for i, c in enumerate(so)}

    err_p = {a: {} for a in ARMS}
    err_s = {a: {} for a in ARMS}
    diag = []

    for f in range(N_FOLDS):
        te_p = np.isin(grp, [c for c in pcomps if pfold[c] == f])
        tr_p = ~te_p
        te_s = (np.isin(scomp, [c for c in scomps if sfold[c] == f])
                if closure == "double" else np.zeros(len(scomp), bool))
        tr_s = ~te_s
        if te_p.sum() == 0 or tr_p.sum() < 10:
            continue

        m = interaction_basis(Y[np.ix_(tr_s, tr_p)], M[np.ix_(tr_s, tr_p)],
                              rank, lam=LAM_LR, seed=1000 + f)
        mu, alpha_tr, U_tr = m["mu"], m["alpha"], m["U"]

        inner = scomp[tr_s]
        Uhat = np.zeros((len(Y), rank))
        Ahat = np.zeros(len(Y))
        p_tr, _ = ridge_cv(X[tr_s], U_tr, X[tr_s], inner)
        a_tr, _ = ridge_cv(X[tr_s], alpha_tr[:, None], X[tr_s], inner)
        Uhat[tr_s], Ahat[tr_s] = p_tr, a_tr.ravel()
        if te_s.any():
            p_te, _ = ridge_cv(X[tr_s], U_tr, X[te_s], inner)
            a_te, _ = ridge_cv(X[tr_s], alpha_tr[:, None], X[te_s], inner)
            Uhat[te_s], Ahat[te_s] = p_te, a_te.ravel()

        U_true = np.zeros((len(Y), rank))
        A_true = np.zeros(len(Y))
        U_true[tr_s], A_true[tr_s] = U_tr, alpha_tr

        te_idx = np.where(te_p)[0]
        for seed in SEEDS:
            rs = np.random.default_rng(10_000 * seed + f)
            nT = len(te_idx)
            der = rs.permutation(nT)
            t_ = 0
            while nT > 1 and np.any(der == np.arange(nT)) and t_ < 500:
                der = rs.permutation(nT)
                t_ += 1

            task = {}
            for t, j in enumerate(te_idx):
                s_pool = np.where(M[:, j] & tr_s)[0]
                q_pool = (np.where(M[:, j] & te_s)[0] if closure == "double"
                          else s_pool)
                if len(s_pool) < k + 5 or len(q_pool) < 5:
                    continue
                pick = rs.permutation(s_pool)
                sup = pick[:k]
                qry = q_pool if closure == "double" else pick[k:]
                if len(qry) < 5:
                    continue
                resid = Y[:, j] - mu - A_true
                task[t] = dict(j=j, sup=sup, qry=qry, resid=resid,
                               bv=solve_bv(resid, U_true, sup, LAM_V),
                               bv_hat=solve_bv(resid, Uhat, sup, LAM_V),
                               b_int=float(resid[sup].mean()))

            for t, T in task.items():
                if der[t] not in task:
                    continue
                j, sup, qry, resid = T["j"], T["sup"], T["qry"], T["resid"]
                yq = Y[qry, j]
                pc = str(grp[j])
                b4, v4 = T["bv"][0], T["bv"][1:]
                bh, vh = T["bv_hat"][0], T["bv_hat"][1:]

                # The support intercept is unpenalised, so it absorbs the mean of
                # the support loadings.  Only the CENTRED support design carries
                # identifiable interaction directions: at most min(k-1, d).  A
                # component of the query loading outside that span is multiplied
                # by a v-coordinate that the ridge prior set to zero, not by
                # anything the support identified.
                U_s = U_true[sup]
                Uc = U_s - U_s.mean(0, keepdims=True)
                sv = np.linalg.svd(Uc, compute_uv=False)
                sv = np.pad(sv, (0, max(0, rank - len(sv))))
                pos = sv[sv > 1e-8 * max(sv.max(), 1e-12)]
                srank = int(len(pos))
                cond = float(pos.max() / pos.min()) if srank else float("inf")
                P = Uc.T @ np.linalg.pinv(Uc.T) if srank else np.zeros((rank, rank))
                nq = np.linalg.norm(Uhat[qry], axis=1)
                cov = np.where(nq > 1e-12,
                               np.einsum("ij,ij->i", Uhat[qry] @ P.T, Uhat[qry])
                               / np.maximum(nq ** 2, 1e-24), 0.0)
                cov = np.clip(cov, 0, 1)

                rp = rs.permutation(sup)
                rperm = resid.copy()
                rperm[sup] = resid[rp]
                vp = solve_bv(rperm, U_true, sup, LAM_V)[1:]
                vfor = task[der[t]]["bv"][1:]
                rr = rs.normal(size=(len(qry), rank))
                rr *= (nq[:, None] / (np.linalg.norm(rr, axis=1, keepdims=True) + 1e-12))
                allsup = np.where(M[:, j] & tr_s)[0]
                bo, vo = solve_bv(resid, U_true, allsup, LAM_V)[0], \
                    solve_bv(resid, U_true, allsup, LAM_V)[1:]

                base = mu + Ahat[qry]
                pred = {
                    "P0": np.full(len(qry), mu),
                    "LIG": base,
                    "ADD": base + T["b_int"],
                    "ZERO": base + b4,
                    "SEC": base + b4 + Uhat[qry] @ v4,
                    "SEC-UHATSUP": base + bh + Uhat[qry] @ vh,
                    "PERM": base + b4 + Uhat[qry] @ vp,
                    "FOREIGN": base + b4 + Uhat[qry] @ vfor,
                    "RANDCORR": base + b4 + rr @ v4,
                    "ORACLE-TRSC": base + bo + Uhat[qry] @ vo,
                }
                sc_q = scomp[qry]
                for arm, pr in pred.items():
                    e = yq - pr
                    s, n = err_p[arm].get(pc, (0.0, 0))
                    err_p[arm][pc] = (s + float((e ** 2).sum()), n + e.size)
                    for sc in set(sc_q):
                        sel = sc_q == sc
                        s, n = err_s[arm].get(sc, (0.0, 0))
                        err_s[arm][sc] = (s + float((e[sel] ** 2).sum()), n + int(sel.sum()))

                diag.append({"fold": f, "seed": seed, "protein": str(d["kinase"][j]),
                             "group": pc, "k": k, "rank": rank,
                             "support_design_rank": srank,
                             "singular_values": [float(x) for x in sv],
                             "condition_number": cond,
                             "identified_section_dim": srank,
                             "query_coverage_mean": float(cov.mean()),
                             "n_query": int(len(qry)),
                             "underidentified": bool(srank < rank)})

    res = {"closure": closure, "rank": rank, "k": k, "ligand_arm": ligand_arm,
           "panel": man["panel"], "seeds": list(SEEDS), "n_folds": N_FOLDS,
           "arms": {}, "contrasts": {}}
    for a in ARMS:
        if not err_p[a]:
            continue
        s = sum(v[0] for v in err_p[a].values())
        n = sum(v[1] for v in err_p[a].values())
        res["arms"][a] = {"rmse": float((s / n) ** 0.5), "cells": int(n)}
        if a != "ADD":
            res["arms"][a]["r2_gamma_vs_ADD"] = widest(
                r2_vs_base(err_p[a], err_p["ADD"], n_boot, 11),
                r2_vs_base(err_s[a], err_s["ADD"], n_boot, 12))
    for name, (arm, base) in {
        "Delta_deploy__LIG_minus_SEC": ("SEC", "LIG"),
        "Delta_interaction__ADD_minus_SEC": ("SEC", "ADD"),
        "Delta_specific_zero__ZERO_minus_SEC": ("SEC", "ZERO"),
        "Delta_specific_foreign__FOREIGN_minus_SEC": ("SEC", "FOREIGN"),
        "Delta_specific_perm__PERM_minus_SEC": ("SEC", "PERM"),
        "Delta_randcorr__RANDCORR_minus_SEC": ("SEC", "RANDCORR"),
        "Delta_oracle__ADD_minus_ORACLE": ("ORACLE-TRSC", "ADD"),
        "Delta_strict__ADD_minus_SECUHATSUP": ("SEC-UHATSUP", "ADD"),
    }.items():
        if err_p[arm] and err_p[base]:
            res["contrasts"][name] = widest(
                paired_contrast(err_p[arm], err_p[base], n_boot, 21),
                paired_contrast(err_s[arm], err_s[base], n_boot, 22))
    if diag:
        cv = [x["query_coverage_mean"] for x in diag]
        res["diagnostics"] = {
            "n_tasks": len(diag),
            "mean_support_design_rank": float(np.mean([x["support_design_rank"] for x in diag])),
            "frac_underidentified": float(np.mean([x["underidentified"] for x in diag])),
            "median_condition_number": float(np.median([x["condition_number"] for x in diag])),
            "mean_query_coverage": float(np.mean(cv)),
            "coverage_percentiles_10_25_50_75_90":
                [float(np.percentile(cv, q)) for q in (10, 25, 50, 75, 90)],
            "sample": diag[:2],
        }
    return res


def pretty(r):
    print(f"\n--- closure={r['closure']} rank={r['rank']} k={r['k']} feat={r['ligand_arm']} ---")
    for a, v in r["arms"].items():
        g = v.get("r2_gamma_vs_ADD")
        gs = (f"{g['point']:+.4f} [{g['ci95'][0]:+.4f},{g['ci95'][1]:+.4f}]" if g
              else "(baseline)")
        print(f"  {a:14s} RMSE={v['rmse']:.4f}  R2_gamma vs ADD {gs}")
    for kk, v in r["contrasts"].items():
        print(f"  {kk:46s} {v['point']:+.5f} [{v['ci95'][0]:+.5f},{v['ci95'][1]:+.5f}]")
    dg = r.get("diagnostics", {})
    if dg:
        print(f"  tasks {dg['n_tasks']}, design rank {dg['mean_support_design_rank']:.2f}/{r['rank']},"
              f" under-identified {dg['frac_underidentified']:.2f},"
              f" median cond {dg['median_condition_number']:.1f},"
              f" mean coverage {dg['mean_query_coverage']:.3f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--closure", default="double")
    ap.add_argument("--rank", type=int, default=3)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--feat", default="L-ECFP")
    ap.add_argument("--tag", default="main")
    a = ap.parse_args()
    r = run(closure=a.closure, rank=a.rank, k=a.k, ligand_arm=a.feat)
    pretty(r)
    p = os.path.join(OUT, f"xp2cd_{a.tag}_{a.closure}_d{a.rank}_k{a.k}_{a.feat}.json")
    json.dump(r, open(p, "w"), indent=2, default=float)
    print("wrote", p)
