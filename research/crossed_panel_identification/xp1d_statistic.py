"""XP1-D: what is the interaction coordinate, and can any protein statistic reach it?

Three diagnostics, all on BLK-METZ-60:

  D1  direct predictability of the interaction coordinate V from each protein
      representation, under leave-one-component-out at three closure levels.
      Separates "features carry nothing" from "features carry only homolog
      interpolation".
  D2  learning curve in the number of training proteins, to separate a
      representation failure from a sample-size failure.
  D3  structure of V: how much of it is explained by kinase group, and which
      aligned KLIFS pocket positions carry univariate association.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lowrank import fit_interaction_basis  # noqa: E402
from panels import load_klifs, load_metz, map_kinases  # noqa: E402
from xp1b_transfer import (REPORT, build_closure, group_onehot,  # noqa: E402
                           pocket_identity_kernel, pocket_onehot,
                           pocket_physchem, ridge_predict)

CACHE = r"D:\MetaSieve\dataset\processed\crossed_panels"
AA = "ACDEFGHIKLMNPQRSTVWY"


def features(kin, pockets, fam, grp, density):
    F = {
        "pocket_identity_kernel": ("kernel", pocket_identity_kernel(pockets)),
        "pocket_onehot": ("dense", pocket_onehot(pockets)),
        "pocket_physchem": ("dense", pocket_physchem(pockets)),
        "group_onehot": ("dense", group_onehot(grp)),
        "family_onehot": ("dense", group_onehot(fam)),
    }
    p = os.path.join(CACHE, f"metz{int(density*100)}_protein_features.npz")
    if os.path.exists(p):
        d = np.load(p, allow_pickle=True)
        for key, col in (("esm2_t30_fullseq", "esm_full"),
                         ("esm2_t30_pocket85", "esm_pocket")):
            X = d[col].astype(float)
            F[key] = ("dense", (X - X.mean(0)) / (X.std(0) + 1e-9))
    p = os.path.join(CACHE, f"metz{int(density*100)}_conformation_features.npz")
    if os.path.exists(p):
        d = np.load(p, allow_pickle=True)
        X = d["X"].astype(float)
        X = np.column_stack([X, (X[:, 0] > 0).astype(float)])
        F["klifs_conformation"] = ("dense", (X - X.mean(0)) / (X.std(0) + 1e-9))
    F["random_gaussian"] = ("dense", np.random.default_rng(0).normal(size=(len(kin), 64)))
    return F


def loco_r2(X, V, comp, kind, lams=(0.1, 1, 10, 100, 1000, 1e4)):
    """Leave-one-closure-component-out R^2 for predicting V from X."""
    comps = sorted(set(comp))
    pred = np.zeros_like(V)
    for c in comps:
        te = comp == c
        tr = ~te
        if tr.sum() < 5:
            pred[te] = V[tr].mean(0)
            continue
        if kind == "kernel":
            Xtr, Xte = X[np.ix_(tr, tr)], X[np.ix_(te, tr)]
        else:
            Xtr, Xte = X[tr], X[te]
        inner = comp[tr]
        best, best_e = lams[0], np.inf
        for lam in lams:
            e = 0.0
            for c2 in sorted(set(inner)):
                a = inner != c2
                if a.sum() < 4:
                    continue
                P = ridge_predict(Xtr[a], V[tr][a], Xtr[~a], lam)
                e += float(((P - V[tr][~a]) ** 2).sum())
            if e < best_e:
                best_e, best = e, lam
        pred[te] = ridge_predict(Xtr, V[tr], Xte, best)
    sse = float(((V - pred) ** 2).sum())
    sst = float(((V - V.mean(0)) ** 2).sum())
    return 1.0 - sse / sst


def main(density=0.60, rank=8, lam_lr=2.0):
    Y, M, cid, kin = load_metz(density)
    hit, miss = map_kinases(kin, load_klifs())
    assert not miss, miss
    pockets = [hit[k]["pocket"] for k in kin]
    fam = [hit[k]["family"] for k in kin]
    grp = [hit[k]["group"] for k in kin]
    Kid = pocket_identity_kernel(pockets)
    F = features(kin, pockets, fam, grp, density)

    m = fit_interaction_basis(Y, M, rank, lam=lam_lr, seed=7)
    V = m["V"]                       # 82 x rank interaction coordinates
    # whiten so every coordinate contributes comparably
    V = (V - V.mean(0)) / (V.std(0) + 1e-9)
    out = {"rank": rank, "density": density, "n_proteins": int(len(kin))}

    print("=" * 78)
    print("D1  leave-one-component-out R^2 for predicting the interaction "
          "coordinate V from protein features")
    print(f"{'representation':26s} {'group':>9s} {'family':>9s} {'pocket60':>9s} "
          f"{'leave-1-protein':>16s}")
    d1 = {}
    closures = {}
    for cl in ("group", "family", "pocket60"):
        closures[cl] = build_closure(kin, hit, cl, Kid)[0]
    closures["single"] = np.array([f"p{i}" for i in range(len(kin))])
    for name, (kind, X) in F.items():
        row = {cl: loco_r2(X, V, closures[cl], kind) for cl in closures}
        d1[name] = row
        print(f"{name:26s} {row['group']:+9.4f} {row['family']:+9.4f} "
              f"{row['pocket60']:+9.4f} {row['single']:+16.4f}")
    out["D1_loco_r2_of_V"] = d1

    print("\n" + "=" * 78)
    print("D2  learning curve: does adding training proteins help the zero-shot map?")
    rng = np.random.default_rng(3)
    comp = closures["group"]
    d2 = {}
    for name in ("pocket_identity_kernel", "esm2_t30_fullseq", "klifs_conformation"):
        if name not in F:
            continue
        kind, X = F[name]
        curve = []
        for n_tr in (20, 35, 50, 65):
            vals = []
            for rep in range(12):
                comps = sorted(set(comp))
                rng.shuffle(comps)
                te_c = comps[:1]
                te = np.isin(comp, te_c)
                tr_pool = np.where(~te)[0]
                if len(tr_pool) < n_tr:
                    continue
                tr = rng.choice(tr_pool, n_tr, replace=False)
                if kind == "kernel":
                    Xtr, Xte = X[np.ix_(tr, tr)], X[np.ix_(np.where(te)[0], tr)]
                else:
                    Xtr, Xte = X[tr], X[te]
                p = ridge_predict(Xtr, V[tr], Xte, 100.0)
                sse = float(((V[te] - p) ** 2).sum())
                sst = float(((V[te] - V[tr].mean(0)) ** 2).sum())
                vals.append(1 - sse / max(sst, 1e-9))
            curve.append((n_tr, float(np.mean(vals)) if vals else np.nan))
        d2[name] = curve
        print(f"  {name:26s} " + "  ".join(f"n={n}: R2={r:+.3f}" for n, r in curve))
    out["D2_learning_curve"] = d2

    print("\n" + "=" * 78)
    print("D3  what is V?")
    G = group_onehot(grp)
    r2_group = 1 - ((V - G @ np.linalg.lstsq(G, V, rcond=None)[0]) ** 2).sum() / \
        ((V - V.mean(0)) ** 2).sum()
    Fam = group_onehot(fam)
    r2_fam = 1 - ((V - Fam @ np.linalg.lstsq(Fam, V, rcond=None)[0]) ** 2).sum() / \
        ((V - V.mean(0)) ** 2).sum()
    print(f"  in-sample R^2 of V on kinase GROUP  (8 levels): {r2_group:.4f}")
    print(f"  in-sample R^2 of V on kinase FAMILY (36 levels): {r2_fam:.4f}")
    P = np.array([list(p) for p in pockets])
    pos_r2 = []
    for k in range(85):
        col = P[:, k]
        if len(set(col)) < 2:
            pos_r2.append(0.0)
            continue
        D = group_onehot(list(col))
        r2 = 1 - ((V - D @ np.linalg.lstsq(D, V, rcond=None)[0]) ** 2).sum() / \
            ((V - V.mean(0)) ** 2).sum()
        pos_r2.append(float(r2))
    order = np.argsort(-np.array(pos_r2))[:10]
    print("  top KLIFS pocket positions by in-sample R^2 of V (1-indexed):")
    print("   " + ", ".join(f"pos{int(i)+1}={pos_r2[i]:.3f}" for i in order))
    print(f"  mean over all 85 positions: {np.mean(pos_r2):.3f} "
          f"(a single 20-level factor on {len(kin)} proteins costs ~"
          f"{19/len(kin):.3f} R^2 by chance)")
    out["D3"] = {"r2_group_insample": float(r2_group),
                 "r2_family_insample": float(r2_fam),
                 "pocket_position_r2": pos_r2,
                 "singular_spectrum": [float(x) for x in
                                       np.linalg.svd(m["U"] @ m["V"].T,
                                                     compute_uv=False)[:20]]}
    p = os.path.join(REPORT, "xp1d_statistic.json")
    json.dump(out, open(p, "w"), indent=2, default=float)
    print("\nwrote", p)


if __name__ == "__main__":
    main()
