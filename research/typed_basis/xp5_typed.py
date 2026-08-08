"""XP5 — fixed named typed interaction basis, double held-out. Frozen by PREREG_XP5.md."""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "crossed_panel_deployability"))
from xp2_core import (additive_fit, paired_contrast, r2_vs_base,  # noqa: E402
                      ridge_predict, widest)
from xp2_panel import build  # noqa: E402

OUT = r"D:\MetaSieve\report\typed_basis"
os.makedirs(OUT, exist_ok=True)
N_FOLDS, LAMS = 5, (0.1, 1.0, 10.0, 100.0, 1e3, 1e4)

ACCEPTOR = set("DENQSTYH")
DONOR = set("KRNQSTYWH")
POS, NEG = set("KR"), set("DE")
HYDRO = set("AVLIMFWC")
AROM = set("FWYH")
POLAR = set("DEHKNQRSTY")
VOL = {"A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9, "G": 60.1,
       "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7, "M": 162.9, "N": 114.1,
       "P": 112.7, "Q": 143.8, "R": 173.4, "S": 89.0, "T": 116.1, "V": 140.0,
       "W": 227.8, "Y": 193.6}
GATEKEEPER = 44          # 0-indexed KLIFS position 45
HINGE = (45, 46, 47)     # 0-indexed KLIFS positions 46-48

CHANNELS = ["hbond_donor_compl", "hbond_acceptor_compl", "electrostatic_anionic",
            "electrostatic_cationic", "hydrophobic_burial", "aromatic_pi",
            "steric_gatekeeper", "steric_hinge", "polar_surface_compl",
            "size_compl"]


def protein_side(pockets):
    out = []
    for p in pockets:
        aa = [c for c in p if c in VOL]
        n = max(len(aa), 1)
        out.append([
            sum(c in ACCEPTOR for c in aa) / n,
            sum(c in DONOR for c in aa) / n,
            (sum(c in POS for c in aa) - sum(c in NEG for c in aa)) / n,
            (sum(c in NEG for c in aa) - sum(c in POS for c in aa)) / n,
            sum(c in HYDRO for c in aa) / n,
            sum(c in AROM for c in aa) / n,
            VOL.get(p[GATEKEEPER], 0.0) if len(p) > GATEKEEPER else 0.0,
            float(np.mean([VOL.get(p[i], 0.0) for i in HINGE if len(p) > i] or [0.0])),
            sum(c in POLAR for c in aa) / n,
            float(np.sum([VOL.get(c, 0.0) for c in aa])),
        ])
    return np.array(out, float)


def ligand_side(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors
    from rdkit.Chem import rdMolDescriptors as rd
    RDLogger.DisableLog("rdApp.*")
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        out.append([
            rd.CalcNumHBD(m), rd.CalcNumHBA(m),
            float(sum(a.GetFormalCharge() for a in m.GetAtoms() if a.GetFormalCharge() < 0)),
            float(sum(a.GetFormalCharge() for a in m.GetAtoms() if a.GetFormalCharge() > 0)),
            Crippen.MolLogP(m), rd.CalcNumAromaticRings(m),
            m.GetNumHeavyAtoms(), rd.CalcNumRotatableBonds(m),
            rd.CalcTPSA(m), Descriptors.MolWt(m),
        ])
    return np.array(out, float)


def run(n_boot=2000, seed=0):
    d, man = build()
    Y, M = d["Y"], d["M"]
    grp, scomp = d["group"], d["scaffold_component"]
    pockets = [str(x) for x in d["pocket"]]
    PS = protein_side(pockets)
    LS = ligand_side(list(d["smiles"]))
    nL, nP = Y.shape
    print(f"panel {nL} ligands x {nP} kinases, {int(M.sum())} cells; "
          f"{len(CHANNELS)} named channels")

    pc, sc = sorted(set(grp)), sorted(set(scomp))
    po = sorted(pc, key=lambda c: -int((grp == c).sum()))
    pfold = {c: i % N_FOLDS for i, c in enumerate(po)}
    so = list(sc)
    np.random.default_rng(seed).shuffle(so)
    sfold = {c: i % N_FOLDS for i, c in enumerate(so)}
    rng = np.random.default_rng(seed)

    ARMS = ["Z0", "TYPED", "RAND-C", "PERM-P", "FOREIGN-P", "SHUF-PAIR"]
    err_p = {a: {} for a in ARMS}
    err_s = {a: {} for a in ARMS}
    coefs = []

    for f in range(N_FOLDS):
        te_p = np.isin(grp, [c for c in pc if pfold[c] == f])
        te_s = np.isin(scomp, [c for c in sc if sfold[c] == f])
        tr_p, tr_s = ~te_p, ~te_s
        Mtr = M & tr_s[:, None] & tr_p[None, :]
        Mte = M & te_s[:, None] & te_p[None, :]
        if Mte.sum() < 30 or Mtr.sum() < 500:
            continue

        # interaction residual from TRAINING cells only
        Yt = np.where(Mtr, Y, 0.0)
        mu, a, b = additive_fit(Y, Mtr)
        G = Y - mu - a[:, None] - b[None, :]

        def design(Pm, Lm, mask):
            i, j = np.where(mask)
            X = Lm[i] * Pm[j]
            return X, i, j

        Xtr, itr, jtr = design(PS, LS, Mtr)
        mX, sX = Xtr.mean(0), Xtr.std(0) + 1e-9
        gtr = G[itr, jtr]
        RC = rng.normal(size=PS.shape)
        PSperm = PS[rng.permutation(nP)]

        inner = grp[jtr]
        best, be = LAMS[0], np.inf
        for lam in LAMS:
            e = 0.0
            for c in sorted(set(inner)):
                m = inner != c
                if m.sum() < 100 or (~m).sum() < 20:
                    continue
                pr = ridge_predict((Xtr[m] - mX) / sX, gtr[m][:, None],
                                   (Xtr[~m] - mX) / sX, lam)[:, 0]
                e += float(((pr - gtr[~m]) ** 2).sum())
            if e < be:
                be, best = e, lam
        lam = best

        ite, jte = np.where(Mte)
        gte = G[ite, jte]
        jfor = jte.copy()
        tp = np.unique(jte)
        if len(tp) > 1:
            shift = {t: tp[(k + 1) % len(tp)] for k, t in enumerate(tp)}
            jfor = np.array([shift[t] for t in jte])
        ishuf = rng.permutation(ite)

        def fit_predict(Pm, Lm, jt, it_):
            Xa = Lm[itr] * Pm[jtr]
            Xb = Lm[it_] * Pm[jt]
            ma, sa = Xa.mean(0), Xa.std(0) + 1e-9
            return ridge_predict((Xa - ma) / sa, gtr[:, None], (Xb - ma) / sa, lam)[:, 0]

        pred = {
            "Z0": np.zeros(len(gte)),
            "TYPED": fit_predict(PS, LS, jte, ite),
            "RAND-C": fit_predict(RC, LS, jte, ite),
            "PERM-P": fit_predict(PSperm, LS, jte, ite),
            "FOREIGN-P": fit_predict(PS, LS, jfor, ite),
            "SHUF-PAIR": fit_predict(PS, LS, jte, ishuf),
        }
        W = ridge_predict((Xtr - mX) / sX, gtr[:, None], np.eye(len(CHANNELS)), lam)
        coefs.append({"fold": f, "lambda": lam,
                      "channel_effect": {c: float(v) for c, v in
                                         zip(CHANNELS, (W[:, 0] - W[:, 0].mean()))}})

        for arm, pr in pred.items():
            e = gte - pr
            for g_ in set(grp[jte]):
                m = grp[jte] == g_
                s, n = err_p[arm].get(g_, (0.0, 0))
                err_p[arm][g_] = (s + float((e[m] ** 2).sum()), n + int(m.sum()))
            for c_ in set(scomp[ite]):
                m = scomp[ite] == c_
                s, n = err_s[arm].get(c_, (0.0, 0))
                err_s[arm][c_] = (s + float((e[m] ** 2).sum()), n + int(m.sum()))

    res = {"stage": "XP5", "panel": man["panel"], "channels": CHANNELS,
           "folds": N_FOLDS, "hyperparameters": coefs, "arms": {}, "contrasts": {}}
    for a in ARMS:
        if not err_p[a]:
            continue
        s = sum(v[0] for v in err_p[a].values())
        n = sum(v[1] for v in err_p[a].values())
        res["arms"][a] = {"rmse": float((s / n) ** 0.5), "cells": int(n)}
        if a != "Z0":
            res["arms"][a]["r2_gamma_vs_Z0"] = widest(
                r2_vs_base(err_p[a], err_p["Z0"], n_boot, 11),
                r2_vs_base(err_s[a], err_s["Z0"], n_boot, 12))
    for nm, (arm, base) in {
        "Delta_null__Z0_minus_TYPED": ("TYPED", "Z0"),
        "Delta_random__RANDC_minus_TYPED": ("TYPED", "RAND-C"),
        "Delta_permP__PERMP_minus_TYPED": ("TYPED", "PERM-P"),
        "Delta_specific__FOREIGNP_minus_TYPED": ("TYPED", "FOREIGN-P"),
        "Delta_pair__SHUFPAIR_minus_TYPED": ("TYPED", "SHUF-PAIR"),
    }.items():
        if err_p[arm] and err_p[base]:
            res["contrasts"][nm] = widest(
                paired_contrast(err_p[arm], err_p[base], n_boot, 21),
                paired_contrast(err_s[arm], err_s[base], n_boot, 22))
    return res


if __name__ == "__main__":
    r = run()
    print(f"\n{'arm':12s} {'RMSE':>8s}  {'R2_gamma vs Z0 [95% CI]':>30s}")
    for a, v in r["arms"].items():
        g = v.get("r2_gamma_vs_Z0")
        gs = (f"{g['point']:+.4f} [{g['ci95'][0]:+.4f},{g['ci95'][1]:+.4f}]" if g
              else "(baseline)")
        print(f"{a:12s} {v['rmse']:8.4f}  {gs:>30s}")
    print("\nregistered contrasts:")
    for k, v in r["contrasts"].items():
        print(f"  {k:44s} {v['point']:+.5f} [{v['ci95'][0]:+.5f},{v['ci95'][1]:+.5f}]")
    p = os.path.join(OUT, "XP5_TYPED_BASIS_RESULT.json")
    json.dump(r, open(p, "w"), indent=2, default=float)
    print("wrote", p)
