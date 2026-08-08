"""XP2-F — external replication on the Klaeger 2017 kinobeads panel.

Independent laboratory, independent measurement technology (chemical proteomics
vs radiometric), independent compound set (clinical drugs vs Abbott series),
never used for any XP1 or XP2 model-selection decision.

Two conclusions are computed and kept separate, per PREREG_XP2 section 12:
  (1) DIRECTION transfer  - every source parameter frozen, no external fitting;
  (2) BASIS transfer      - plus one preregistered global affine calibration
                            fitted on external TRAINING components only.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "crossed_panel_identification"))
from xp2_core import (interaction_basis, paired_contrast, r2_vs_base,  # noqa: E402
                      ridge_cv, widest)
from xp2_panel import CACHE, build, ligand_features  # noqa: E402

RAW = r"D:\MetaSieve\dataset\raw\crossed_panels\kinase_panels"
OUT = r"D:\MetaSieve\report\crossed_panel_deployability"
os.makedirs(OUT, exist_ok=True)
KLAEGER_FLOOR = 5.0
MIN_HITS_KIN, MIN_HITS_DRUG = 15, 5
RANK, K, LAM_V = 3, 5, 1.0
SEEDS = (0, 1, 2, 3, 4)


def build_klaeger():
    import pandas as pd
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    from panels import load_klifs, map_kinases

    df = pd.read_csv(os.path.join(RAW, "klaeger_matrix.csv"), low_memory=False)
    drugs = df.iloc[:, 0].astype(str).str.strip().to_numpy()
    kin = np.array([c.strip().upper() for c in df.columns[1:]])
    W = df.iloc[:, 1:].to_numpy(float)
    hit = W > KLAEGER_FLOOR + 1e-9

    smi_map = json.load(open(os.path.join(CACHE, "klaeger_smiles.json")))
    smiles = np.array([(smi_map.get(d) or {}).get("smiles") for d in drugs], dtype=object)
    ok_d = np.array([isinstance(s, str) and Chem.MolFromSmiles(s) is not None
                     for s in smiles])
    hitk, miss = map_kinases(kin, load_klifs())
    ok_k = np.array([k in hitk for k in kin])

    keep_d, keep_k = ok_d.copy(), ok_k.copy()
    for _ in range(50):
        nd = hit[np.ix_(keep_d, keep_k)].sum(1)
        i = np.where(keep_d)[0]
        kd = keep_d.copy()
        kd[i[nd < MIN_HITS_DRUG]] = False
        nk = hit[np.ix_(kd, keep_k)].sum(0)
        j = np.where(keep_k)[0]
        kk = keep_k.copy()
        kk[j[nk < MIN_HITS_KIN]] = False
        if (kd == keep_d).all() and (kk == keep_k).all():
            break
        keep_d, keep_k = kd, kk

    di, ki = np.where(keep_d)[0], np.where(keep_k)[0]
    Y = W[np.ix_(di, ki)]
    M = hit[np.ix_(di, ki)]
    sm = smiles[di]
    scaf = [MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s)) for s in sm]
    uq = sorted(set(scaf))
    sc = np.array([f"kc{uq.index(s)}" for s in scaf])
    grp = np.array([hitk[k]["group"] for k in kin[ki]])
    man = {"panel": "BLK-KLAEGER-XP2", "n_drugs": int(len(di)), "n_kinases": int(len(ki)),
           "measured_cells": int(M.sum()), "density": float(M.mean()),
           "n_scaffolds": int(len(uq)), "floor": KLAEGER_FLOOR,
           "endpoint": "pKd_app (kinobeads); cells at the 5.0 floor are excluded, "
                       "mirroring the measured-only convention used on Metz",
           "drugs_unresolved_structure": int((~ok_d).sum()),
           "kinases_without_klifs": int((~ok_k).sum())}
    return dict(Y=Y, M=M, smiles=sm, drug=drugs[di], kinase=kin[ki], group=grp,
                scaffold_component=sc), man


def source_model(rank=RANK):
    """Freeze the Metz-side model: mu, chi_alpha, chi_u, all on the full XP2 panel."""
    d, man = build()
    X = ligand_features(d["smiles"], kinds=("ecfp",))["L-ECFP"]
    m = interaction_basis(d["Y"], d["M"], rank, lam=2.0, seed=7)
    inner = d["scaffold_component"]
    # frozen maps: fitted on ALL Metz ligands (Klaeger is the external set)
    def chi_u(Xq):
        return ridge_cv(X, m["U"], Xq, inner)[0]

    def chi_a(Xq):
        return ridge_cv(X, m["alpha"][:, None], Xq, inner)[0].ravel()
    return m["mu"], chi_u, chi_a, man


def solve_bv(resid, U, idx, lam):
    D = np.column_stack([np.ones(len(idx)), U[idx]])
    pen = np.diag([0.0] + [lam] * U.shape[1])
    return np.linalg.solve(D.T @ D + pen, D.T @ resid[idx])


def run(n_boot=2000):
    ext, eman = build_klaeger()
    print("external panel:", json.dumps(eman, indent=2))
    mu_src, chi_u, chi_a, sman = source_model()
    Xe = ligand_features(ext["smiles"], kinds=("ecfp",))["L-ECFP"]
    Uh = chi_u(Xe)
    Ah = chi_a(Xe)
    Y, M, sc = ext["Y"], ext["M"], ext["scaffold_component"]
    grp = ext["group"]

    arms = ["P0", "LIG", "ADD", "SEC", "FOREIGN", "PERM"]
    err_p = {a: {} for a in arms}
    err_s = {a: {} for a in arms}
    calib_pairs = []          # (prediction, truth) on external TRAINING components
    scomps = sorted(set(sc))
    np.random.default_rng(0).shuffle(scomps)
    cal_c = set(scomps[: max(1, len(scomps) // 3)])      # calibration third
    test_mask_scaffold = np.array([s not in cal_c for s in sc])

    for seed in SEEDS:
        rs = np.random.default_rng(seed)
        nK = Y.shape[1]
        der = rs.permutation(nK)
        t_ = 0
        while nK > 1 and np.any(der == np.arange(nK)) and t_ < 500:
            der = rs.permutation(nK)
            t_ += 1
        sols = {}
        for j in range(nK):
            obs = np.where(M[:, j])[0]
            if len(obs) < K + 5:
                continue
            # support from scaffold components disjoint from the query scaffold
            pick = rs.permutation(obs)
            sup = pick[:K]
            sup_sc = set(sc[sup])
            rest = np.array([q for q in pick[K:] if sc[q] not in sup_sc])
            # calibration queries and test queries are scaffold-disjoint by
            # construction: the calibration third of scaffold components is
            # reserved and never scored.
            qry = np.array([q for q in rest if test_mask_scaffold[q]], dtype=int)
            cal = np.array([q for q in rest if not test_mask_scaffold[q]], dtype=int)
            if len(qry) < 5:
                continue
            resid = Y[:, j] - mu_src - Ah
            sols[j] = (sup, qry, resid, solve_bv(resid, Uh, sup, LAM_V),
                       float(resid[sup].mean()), cal)
        for j, (sup, qry, resid, sol, b_int, cal) in sols.items():
            if der[j] not in sols:
                continue
            b4, v4 = sol[0], sol[1:]
            base = mu_src + Ah[qry]
            yq = Y[qry, j]
            rp = rs.permutation(sup)
            rperm = resid.copy()
            rperm[sup] = resid[rp]
            vp = solve_bv(rperm, Uh, sup, LAM_V)[1:]
            vf = sols[der[j]][3][1:]
            pred = {"P0": np.full(len(qry), mu_src), "LIG": base,
                    "ADD": base + b_int, "SEC": base + b4 + Uh[qry] @ v4,
                    "FOREIGN": base + b4 + Uh[qry] @ vf,
                    "PERM": base + b4 + Uh[qry] @ vp}
            pc, sq = str(grp[j]), sc[qry]
            for a, pr in pred.items():
                e = yq - pr
                s, n = err_p[a].get(pc, (0.0, 0))
                err_p[a][pc] = (s + float((e ** 2).sum()), n + e.size)
                for s_ in set(sq):
                    m_ = sq == s_
                    s2, n2 = err_s[a].get(s_, (0.0, 0))
                    err_s[a][s_] = (s2 + float((e[m_] ** 2).sum()), n2 + int(m_.sum()))
            if seed == SEEDS[0] and len(cal) >= 3:
                pc_cal = mu_src + Ah[cal] + b4 + Uh[cal] @ v4
                calib_pairs.append((pc_cal, Y[cal, j]))

    res = {"stage": "XP2-F", "external_panel": eman, "source_panel": sman["panel"],
           "conclusion_1_direction_transfer": {"arms": {}, "contrasts": {}}}
    D1 = res["conclusion_1_direction_transfer"]
    for a in arms:
        if not err_p[a]:
            continue
        s = sum(v[0] for v in err_p[a].values())
        n = sum(v[1] for v in err_p[a].values())
        D1["arms"][a] = {"rmse": float((s / n) ** 0.5), "cells": int(n)}
        if a != "ADD":
            D1["arms"][a]["r2_gamma_vs_ADD"] = widest(
                r2_vs_base(err_p[a], err_p["ADD"], n_boot, 11),
                r2_vs_base(err_s[a], err_s["ADD"], n_boot, 12))
    for nm, (a, b) in {"Delta_interaction__ADD_minus_SEC": ("SEC", "ADD"),
                       "Delta_specific_foreign__FOREIGN_minus_SEC": ("SEC", "FOREIGN"),
                       "Delta_specific_perm__PERM_minus_SEC": ("SEC", "PERM"),
                       "Delta_deploy__LIG_minus_SEC": ("SEC", "LIG")}.items():
        if err_p[a] and err_p[b]:
            D1["contrasts"][nm] = widest(
                paired_contrast(err_p[a], err_p[b], n_boot, 21),
                paired_contrast(err_s[a], err_s[b], n_boot, 22))

    # ---- conclusion 2: one global affine calibration on external TRAINING scaffolds
    P = np.concatenate([p for p, _ in calib_pairs]) if calib_pairs else np.array([])
    T = np.concatenate([t for _, t in calib_pairs]) if calib_pairs else np.array([])
    if P.size > 10:
        A = np.column_stack([np.ones(len(P)), P])
        coef = np.linalg.lstsq(A, T, rcond=None)[0]
        res["conclusion_2_basis_transfer"] = {
            "affine_a": float(coef[0]), "affine_b": float(coef[1]),
            "fitted_on": "external calibration scaffold components only",
            "note": "reported separately; never merged with conclusion 1"}
    return res


if __name__ == "__main__":
    r = run()
    D1 = r["conclusion_1_direction_transfer"]
    print("\n--- XP2-F conclusion 1: DIRECTION transfer (all source params frozen) ---")
    for a, v in D1["arms"].items():
        g = v.get("r2_gamma_vs_ADD")
        gs = (f"{g['point']:+.4f} [{g['ci95'][0]:+.4f},{g['ci95'][1]:+.4f}]" if g
              else "(baseline)")
        print(f"  {a:10s} RMSE={v['rmse']:.4f}  R2_gamma vs ADD {gs}")
    for k, v in D1["contrasts"].items():
        print(f"  {k:44s} {v['point']:+.5f} [{v['ci95'][0]:+.5f},{v['ci95'][1]:+.5f}]")
    if "conclusion_2_basis_transfer" in r:
        c = r["conclusion_2_basis_transfer"]
        print(f"\n--- conclusion 2: BASIS transfer affine a={c['affine_a']:.3f} "
              f"b={c['affine_b']:.3f} ---")
    p = os.path.join(OUT, "EXTERNAL_REPLICATION_RESULT.json")
    json.dump(r, open(p, "w"), indent=2, default=float)
    print("wrote", p)
