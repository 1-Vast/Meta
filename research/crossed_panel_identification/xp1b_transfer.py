"""XP1-B: does an unseen protein's interaction section transfer?

Every arm predicts the SAME held-out cells of the SAME held-out proteins, uses
the SAME ligand main effect fitted on training proteins only, and the SAME
low-rank ligand loading basis U fitted on training proteins only.  Arms differ
only in how the held-out protein's interaction coordinate v_j is obtained.

Registered in PREREG_XP1.md sections 5, 6, 7.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panels import load_klifs, load_metz, map_kinases  # noqa: E402
from lowrank import fit_interaction_basis  # noqa: E402

REPORT = r"D:\MetaSieve\report\crossed_panel_identification"
os.makedirs(REPORT, exist_ok=True)

AA = "ACDEFGHIKLMNPQRSTVWY"
PROP = {  # Kyte-Doolittle hydropathy, formal charge, volume, HB donors, HB acceptors
    "A": (1.8, 0, 88.6, 0, 0), "C": (2.5, 0, 108.5, 1, 0), "D": (-3.5, -1, 111.1, 0, 2),
    "E": (-3.5, -1, 138.4, 0, 2), "F": (2.8, 0, 189.9, 0, 0), "G": (-0.4, 0, 60.1, 0, 0),
    "H": (-3.2, 0.1, 153.2, 1, 1), "I": (4.5, 0, 166.7, 0, 0), "K": (-3.9, 1, 168.6, 1, 0),
    "L": (3.8, 0, 166.7, 0, 0), "M": (1.9, 0, 162.9, 0, 0), "N": (-3.5, 0, 114.1, 1, 1),
    "P": (-1.6, 0, 112.7, 0, 0), "Q": (-3.5, 0, 143.8, 1, 1), "R": (-4.5, 1, 173.4, 2, 0),
    "S": (-0.8, 0, 89.0, 1, 1), "T": (-0.7, 0, 116.1, 1, 1), "V": (4.2, 0, 140.0, 0, 0),
    "W": (-0.9, 0, 227.8, 1, 0), "Y": (-1.3, 0, 193.6, 1, 1),
}
DEFAULT = (0.0, 0.0, 0.0, 0.0, 0.0)


# ----------------------------------------------------------------- features
def pocket_identity_kernel(pockets):
    P = np.array([list(p) for p in pockets])
    return np.stack([(P == P[i]).mean(axis=1) for i in range(len(pockets))])


def pocket_onehot(pockets):
    X = np.zeros((len(pockets), 85 * len(AA)))
    for i, p in enumerate(pockets):
        for k, ch in enumerate(p):
            j = AA.find(ch)
            if j >= 0:
                X[i, k * len(AA) + j] = 1.0
    return X


def pocket_physchem(pockets):
    X = np.zeros((len(pockets), 85 * 5))
    for i, p in enumerate(pockets):
        for k, ch in enumerate(p):
            X[i, k * 5:(k + 1) * 5] = PROP.get(ch, DEFAULT)
    return (X - X.mean(0)) / (X.std(0) + 1e-9)


def group_onehot(labels):
    lv = sorted(set(labels))
    X = np.zeros((len(labels), len(lv)))
    for i, v in enumerate(labels):
        X[i, lv.index(v)] = 1.0
    return X


# ----------------------------------------------------------------- estimators
def ridge_fit(X, T, lam):
    return np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T @ T)


def ridge_predict(Xtr, Ttr, Xte, lam):
    """Ridge with an unpenalised intercept; dual form when d > n (d up to 1700)."""
    xm, tm = Xtr.mean(0), Ttr.mean(0)
    A, B, C = Xtr - xm, Ttr - tm, Xte - xm
    n, d = A.shape
    if d > n:
        K = A @ A.T
        return tm + C @ (A.T @ np.linalg.solve(K + lam * np.eye(n), B))
    return tm + C @ np.linalg.solve(A.T @ A + lam * np.eye(d), A.T @ B)


def solve_bv(resid, U, idx, lam_v):
    """Joint ridge solve of the intercept b and the interaction coordinate v."""
    D = np.column_stack([np.ones(len(idx)), U[idx]])
    pen = np.diag([0.0] + [lam_v] * U.shape[1])
    sol = np.linalg.solve(D.T @ D + pen, D.T @ resid[idx])
    return float(sol[0]), sol[1:]


# ----------------------------------------------------------------- inference
def _boot_indices(n_comp, n_boot, seed):
    rng = np.random.default_rng(seed)
    return rng.integers(0, n_comp, size=(n_boot, n_comp))


def contrast(err_a, err_b, comps, bidx):
    """Paired MSE gain of arm a over arm b, cluster-bootstrapped over components."""
    sa = np.array([err_a[c][0] for c in comps])
    sb = np.array([err_b[c][0] for c in comps])
    nn = np.array([err_a[c][1] for c in comps])
    tot_n = nn.sum()
    point = (sb.sum() - sa.sum()) / tot_n
    bs = (sb[bidx].sum(1) - sa[bidx].sum(1)) / nn[bidx].sum(1)
    return {"point": float(point),
            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "rmse_b": float((sb.sum() / tot_n) ** 0.5),
            "rmse_a": float((sa.sum() / tot_n) ** 0.5)}


def r2_gamma(err_a, err_base, comps, bidx):
    sa = np.array([err_a[c][0] for c in comps])
    sb = np.array([err_base[c][0] for c in comps])
    bs = 1.0 - sa[bidx].sum(1) / sb[bidx].sum(1)
    return {"point": float(1.0 - sa.sum() / sb.sum()),
            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]}


# ----------------------------------------------------------------- main
def build_closure(kin, hit, kind, Kid):
    fam = [hit[k]["family"] for k in kin]
    grp = [hit[k]["group"] for k in kin]
    if kind == "group":
        return np.array(grp), fam, grp
    if kind == "family":
        return np.array(fam), fam, grp
    if kind.startswith("pocket"):
        thr = float(kind.replace("pocket", "")) / 100.0
        n = len(kin)
        lab = np.arange(n)
        for i in range(n):
            for j in range(i + 1, n):
                if Kid[i, j] >= thr:
                    lab[lab == lab[j]] = lab[i]
        return np.array([f"pc{v}" for v in lab]), fam, grp
    raise ValueError(kind)


def run(closure="group", rank=8, k_support=16, seeds=(0, 1, 2, 3, 4),
        lam_lr=2.0, lam_v=1.0, density=0.60, n_boot=2000, verbose=True):
    Y, M, cid, kin = load_metz(density)
    hit, miss = map_kinases(kin, load_klifs())
    if miss:
        raise RuntimeError(f"unmapped kinases: {miss}")
    pockets = [hit[k]["pocket"] for k in kin]
    Kid = pocket_identity_kernel(pockets)
    comp, fam, grp = build_closure(kin, hit, closure, Kid)

    FEAT = {
        "pocket_identity_kernel": Kid,
        "pocket_onehot": pocket_onehot(pockets),
        "pocket_physchem": pocket_physchem(pockets),
        "group_onehot": group_onehot(grp),
        "family_onehot": group_onehot(fam),
    }
    # the production MetaSieve protein encoder, if its cache has been built
    fp = os.path.join(r"D:\MetaSieve\dataset\processed\crossed_panels",
                      f"metz{int(density*100)}_protein_features.npz")
    if os.path.exists(fp):
        d = np.load(fp, allow_pickle=True)
        if list(d["kinases"]) == list(kin):
            for key, col in (("esm2_t30_fullseq", "esm_full"),
                             ("esm2_t30_pocket85", "esm_pocket")):
                X = d[col].astype(float)
                FEAT[key] = (X - X.mean(0)) / (X.std(0) + 1e-9)
    cf = os.path.join(r"D:\MetaSieve\dataset\processed\crossed_panels",
                      f"metz{int(density*100)}_conformation_features.npz")
    if os.path.exists(cf):
        d = np.load(cf, allow_pickle=True)
        if list(d["kinases"]) == list(kin):
            X = d["X"].astype(float)
            X = np.column_stack([X, (X[:, 0] > 0).astype(float)])  # has-structure flag
            FEAT["klifs_conformation"] = (X - X.mean(0)) / (X.std(0) + 1e-9)
    comps_all = sorted(set(comp))
    order = sorted(comps_all, key=lambda c: -int((comp == c).sum()))
    folds = {c: i % 5 for i, c in enumerate(order)}
    if verbose:
        print(f"closure={closure}: {len(comps_all)} components / {len(kin)} kinases; "
              f"fold sizes {[sum((comp == c).sum() for c in comps_all if folds[c]==f) for f in range(5)]}")

    arms = ["A0", "A1", "A2", "A4", "A6", "A7", "A34", "AO1",
            "A3B::knn_pocket", "A5B::knn_pocket"] + \
           [f"A3::{f}" for f in FEAT] + [f"A5::{f}" for f in FEAT]
    err = {a: {} for a in arms}

    def add(arm, c, e):
        s, n = err[arm].get(c, (0.0, 0))
        err[arm][c] = (s + float((e ** 2).sum()), n + len(e))

    rng_g = np.random.default_rng(12345)

    for f in range(5):
        te_comp = [c for c in comps_all if folds[c] == f]
        te_idx = np.where(np.isin(comp, te_comp))[0]
        tr_idx = np.where(~np.isin(comp, te_comp))[0]
        if len(te_idx) == 0 or len(tr_idx) < 10:
            continue

        m = fit_interaction_basis(Y[:, tr_idx], M[:, tr_idx], rank, lam=lam_lr, seed=1000 + f)
        mu, alpha, U, Vtr = m["mu"], m["alpha"], m["U"], m["V"]

        # --- basis-free zero-shot: homolog-kernel average of training interaction
        #     residual COLUMNS (no low-rank basis, no ridge in v-space).
        add_tr = mu + alpha[:, None] + m["beta"][None, :]
        Gtr = np.where(M[:, tr_idx], Y[:, tr_idx] - add_tr, np.nan)
        Ktr = Kid[np.ix_(tr_idx, tr_idx)]
        Kte = Kid[np.ix_(te_idx, tr_idx)]

        def knn_column(krow, K, exclude=None):
            w = krow.copy()
            if exclude is not None:
                w[exclude] = -np.inf
            top = np.argsort(-w)[:K]
            Wt = np.clip(w[top], 0, None)
            if Wt.sum() <= 0:
                Wt = np.ones_like(Wt)
            G = Gtr[:, top]
            ok = np.isfinite(G)
            num = np.nansum(np.where(ok, G, 0) * Wt[None, :], axis=1)
            den = (ok * Wt[None, :]).sum(axis=1)
            return np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)

        best_K, best_e = 1, np.inf
        for K in (1, 3, 5, 10):
            e = 0.0
            for tt in range(len(tr_idx)):
                same = np.where(comp[tr_idx] == comp[tr_idx][tt])[0]
                pred = knn_column(Ktr[tt], K, exclude=same)
                ok = np.isfinite(Gtr[:, tt])
                e += float(((Gtr[ok, tt] - pred[ok]) ** 2).sum())
            if e < best_e:
                best_e, best_K = e, K
        knn_pred = {t: knn_column(Kte[t], best_K) for t in range(len(te_idx))}

        # zero-shot maps: protein features -> v, ridge, lambda by leave-one-component-out
        vpred, lam_chosen = {}, {}
        inner = comp[tr_idx]
        for fname, X in FEAT.items():
            if fname == "pocket_identity_kernel":
                Xtr, Xte = X[np.ix_(tr_idx, tr_idx)], X[np.ix_(te_idx, tr_idx)]
            else:
                Xtr, Xte = X[tr_idx], X[te_idx]
            best, best_e = 1.0, np.inf
            for lam in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0):
                e = 0.0
                for c in sorted(set(inner)):
                    a = inner != c
                    if a.sum() < 5:
                        continue
                    P = ridge_predict(Xtr[a], Vtr[a], Xtr[~a], lam)
                    e += float(((P - Vtr[~a]) ** 2).sum())
                if e < best_e:
                    best_e, best = e, lam
            vpred[fname] = ridge_predict(Xtr, Vtr, Xte, best)
            lam_chosen[fname] = best

        rnd = rng_g.normal(size=(Y.shape[1], 64))
        v_rand = ridge_predict(rnd[tr_idx], Vtr, rnd[te_idx], 10.0)

        for seed in seeds:
            rs = np.random.default_rng(1000 * seed + f)
            nT = len(te_idx)
            der = rs.permutation(nT)
            tries = 0
            while nT > 1 and np.any(der == np.arange(nT)) and tries < 500:
                der = rs.permutation(nT)
                tries += 1

            sup, tst, bv, b2 = {}, {}, {}, {}
            for t, j in enumerate(te_idx):
                obs = np.where(M[:, j])[0]
                if len(obs) < k_support + 20:
                    continue
                pick = rs.permutation(obs)
                sup[t], tst[t] = pick[:k_support], pick[k_support:]
                resid = Y[:, j] - mu - alpha
                bv[t] = solve_bv(resid, U, sup[t], lam_v)
                b2[t] = float(resid[sup[t]].mean())

            valid = [t for t in tst if der[t] in tst]
            for t in valid:
                j = te_idx[t]
                ti = tst[t]
                yte = Y[ti, j]
                c = comp[j]
                base = mu + alpha[ti]
                add("A0", c, yte - mu)
                add("A1", c, yte - base)
                a2 = base + b2[t]                       # additive + support location
                add("A2", c, yte - a2)
                for fname in FEAT:
                    add(f"A3::{fname}", c, yte - (a2 + U[ti] @ vpred[fname][t]))
                    add(f"A5::{fname}", c, yte - (a2 + U[ti] @ vpred[fname][der[t]]))
                add("A3B::knn_pocket", c, yte - (a2 + knn_pred[t][ti]))
                add("A5B::knn_pocket", c, yte - (a2 + knn_pred[der[t]][ti]))
                b4, v4 = bv[t]
                add("A4", c, yte - (base + b4 + U[ti] @ v4))
                b6, v6 = bv[der[t]]
                add("A6", c, yte - (base + b4 + U[ti] @ v6))
                add("A7", c, yte - (a2 + U[ti] @ v_rand[t]))
                vk = vpred["pocket_identity_kernel"][t]
                add("A34", c, yte - (base + b4 + U[ti] @ (0.5 * (vk + v4))))
                resid = Y[:, j] - mu - alpha
                bo, vo = solve_bv(resid, U, np.where(M[:, j])[0], lam_v)
                add("AO1", c, yte - (base + bo + U[ti] @ vo))

    comps = sorted(err["A2"])
    bidx = _boot_indices(len(comps), n_boot, seed=99)
    out = {"closure": closure, "rank": rank, "k_support": k_support,
           "density": density, "n_components": len(comps_all),
           "n_eval_components": len(comps), "seeds": list(seeds), "arms": {}}
    for a in arms:
        if not err[a]:
            continue
        s = sum(v[0] for v in err[a].values())
        n = sum(v[1] for v in err[a].values())
        out["arms"][a] = {"rmse": float((s / n) ** 0.5), "cells": int(n)}
        if a != "A2":
            out["arms"][a]["r2_gamma_vs_A2"] = r2_gamma(err[a], err["A2"], comps, bidx)
    A3 = "A3::pocket_identity_kernel"
    A5 = "A5::pocket_identity_kernel"
    out["contrasts"] = {
        "Delta_protein__A1_minus_A3": contrast(err[A3], err["A1"], comps, bidx),
        "Delta_protein__A1_minus_A4": contrast(err["A4"], err["A1"], comps, bidx),
        "Delta_interaction__A2_minus_A3": contrast(err[A3], err["A2"], comps, bidx),
        "Delta_interaction__A2_minus_A4": contrast(err["A4"], err["A2"], comps, bidx),
        "Delta_specific__A5_minus_A3": contrast(err[A3], err[A5], comps, bidx),
        "Delta_specific__A6_minus_A4": contrast(err["A4"], err["A6"], comps, bidx),
        "Delta_featurenull__A7_minus_A3": contrast(err[A3], err["A7"], comps, bidx),
        "Delta_oracle__A2_minus_AO1": contrast(err["AO1"], err["A2"], comps, bidx),
        "Delta_interaction__A2_minus_A3Bknn": contrast(err["A3B::knn_pocket"], err["A2"],
                                                       comps, bidx),
        "Delta_specific__A5Bknn_minus_A3Bknn": contrast(err["A3B::knn_pocket"],
                                                        err["A5B::knn_pocket"], comps, bidx),
    }
    return out


def pretty(out):
    print(f"\n--- closure={out['closure']}  rank={out['rank']}  k={out['k_support']} "
          f"({out['n_eval_components']} evaluated components) ---")
    print(f"{'arm':32s} {'RMSE':>8s}  {'R2_gamma vs A2 [95% CI]':>30s}")
    for a, v in out["arms"].items():
        r = v.get("r2_gamma_vs_A2")
        rs = (f"{r['point']:+.4f} [{r['ci95'][0]:+.4f},{r['ci95'][1]:+.4f}]"
              if r else "(baseline)")
        print(f"{a:32s} {v['rmse']:8.4f}  {rs:>30s}")
    print("  registered contrasts (dMSE > 0 = first-named arm better):")
    for k, v in out["contrasts"].items():
        print(f"   {k:34s} {v['point']:+.5f} [{v['ci95'][0]:+.5f},{v['ci95'][1]:+.5f}]"
              f"   RMSE {v['rmse_b']:.4f} -> {v['rmse_a']:.4f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--closure", default="group")
    ap.add_argument("--rank", type=int, default=8)
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--density", type=float, default=0.60)
    ap.add_argument("--tag", default="main")
    a = ap.parse_args()
    res = run(closure=a.closure, rank=a.rank, k_support=a.k, density=a.density)
    pretty(res)
    p = os.path.join(REPORT, f"xp1b_{a.tag}_{a.closure}_r{a.rank}_k{a.k}.json")
    with open(p, "w") as f:
        json.dump(res, f, indent=2, default=float)
    print("wrote", p)
