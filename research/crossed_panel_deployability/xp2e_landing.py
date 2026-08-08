"""XP2-E — biological landing comparison, without architecture expansion.

Does biology contribute to the section, or do the dataset axes contribute?
Same double closure, same cells, same ligand chemistry landing; the arms differ
only in where the PROTEIN-side coordinate v_j comes from.
"""
from __future__ import annotations

import json
import os
import ssl
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "crossed_panel_identification"))
from xp2_core import (interaction_basis, paired_contrast, r2_vs_base,  # noqa: E402
                      ridge_cv, ridge_predict, widest)
from xp2_panel import CACHE, build, ligand_features  # noqa: E402

OUT = r"D:\MetaSieve\report\crossed_panel_deployability"
os.makedirs(OUT, exist_ok=True)
LAM_LR, LAM_V = 2.0, 1.0
N_FOLDS, SEEDS = 5, (0, 1, 2, 3, 4)
ARMS = ["ADD", "SUP", "PROT-ESM", "PROT-POCKET", "PROT-GROUP",
        "PROT-ESM+SUP", "PROT-POCKET+SUP", "PROT-RANDOM"]


def protein_features(d):
    pockets = [str(x) for x in d["pocket"]]
    P = np.array([list(p) for p in pockets])
    Kid = np.stack([(P == P[i]).mean(1) for i in range(len(P))])
    grp = [str(x) for x in d["group"]]
    lv = sorted(set(grp))
    G = np.zeros((len(grp), len(lv)))
    for i, v in enumerate(grp):
        G[i, lv.index(v)] = 1.0
    E = _esm([str(u) for u in d["uniprot"]])
    E = (E - E.mean(0)) / (E.std(0) + 1e-9)
    R = np.random.default_rng(0).normal(size=(len(grp), 64))
    return {"PROT-ESM": ("dense", E), "PROT-POCKET": ("kernel", Kid),
            "PROT-GROUP": ("dense", G), "PROT-RANDOM": ("dense", R)}


def _esm(uniprot):
    path = os.path.join(CACHE, "esm2_t30_xp2_kinases.npz")
    if os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        m = {k: v for k, v in zip(z["keys"], z["emb"])}
        if all(u in m for u in uniprot):
            return np.stack([m[u] for u in uniprot])
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    sp = os.path.join(CACHE, "uniprot_seq_xp2.json")
    seqs = json.load(open(sp)) if os.path.exists(sp) else {}
    for u in uniprot:
        if u in seqs:
            continue
        req = urllib.request.Request(f"https://rest.uniprot.org/uniprotkb/{u}.fasta",
                                     headers={"User-Agent": "Mozilla/5.0"})
        t = urllib.request.urlopen(req, timeout=60, context=ctx).read().decode()
        seqs[u] = "".join(l.strip() for l in t.splitlines() if not l.startswith(">"))
    json.dump(seqs, open(sp, "w"))
    import torch
    from transformers import AutoTokenizer, EsmModel
    name = "facebook/esm2_t30_150M_UR50D"
    tok = AutoTokenizer.from_pretrained(name)
    mdl = EsmModel.from_pretrained(name).eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    mdl = mdl.to(dev)
    keys, embs = [], []
    with torch.no_grad():
        for i in range(0, len(uniprot), 4):
            ch = list(dict.fromkeys(uniprot))[i:i + 4] if False else uniprot[i:i + 4]
            enc = tok([seqs[c][:1024] for c in ch], return_tensors="pt", padding=True,
                      truncation=True, max_length=1024).to(dev)
            o = mdl(**enc).last_hidden_state
            mk = enc["attention_mask"].unsqueeze(-1).float()
            for c, e in zip(ch, ((o * mk).sum(1) / mk.sum(1)).cpu().numpy()):
                keys.append(c)
                embs.append(e)
    np.savez(path, keys=np.array(keys), emb=np.stack(embs))
    m = {k: v for k, v in zip(keys, embs)}
    return np.stack([m[u] for u in uniprot])


def solve_bv(resid, U, idx, lam):
    D = np.column_stack([np.ones(len(idx)), U[idx]])
    pen = np.diag([0.0] + [lam] * U.shape[1])
    return np.linalg.solve(D.T @ D + pen, D.T @ resid[idx])


def run(rank=3, k=5, ligand_arm="L-ECFP", n_boot=2000):
    d, man = build()
    Y, M = d["Y"], d["M"]
    grp, scomp = d["group"], d["scaffold_component"]
    X = ligand_features(d["smiles"], kinds=("ecfp",))[ligand_arm]
    PF = protein_features(d)

    pcomps, scomps = sorted(set(grp)), sorted(set(scomp))
    po = sorted(pcomps, key=lambda c: -int((grp == c).sum()))
    pfold = {c: i % N_FOLDS for i, c in enumerate(po)}
    so = list(scomps)
    np.random.default_rng(0).shuffle(so)
    sfold = {c: i % N_FOLDS for i, c in enumerate(so)}

    err_p = {a: {} for a in ARMS}
    err_s = {a: {} for a in ARMS}

    for f in range(N_FOLDS):
        te_p = np.isin(grp, [c for c in pcomps if pfold[c] == f])
        tr_p = ~te_p
        te_s = np.isin(scomp, [c for c in scomps if sfold[c] == f])
        tr_s = ~te_s
        if te_p.sum() == 0 or tr_p.sum() < 10:
            continue
        m = interaction_basis(Y[np.ix_(tr_s, tr_p)], M[np.ix_(tr_s, tr_p)],
                              rank, lam=LAM_LR, seed=1000 + f)
        mu, alpha_tr, U_tr, V_tr = m["mu"], m["alpha"], m["U"], m["V"]
        inner = scomp[tr_s]
        Uhat = np.zeros((len(Y), rank))
        Ahat = np.zeros(len(Y))
        Uhat[tr_s], _ = ridge_cv(X[tr_s], U_tr, X[tr_s], inner)
        a_, _ = ridge_cv(X[tr_s], alpha_tr[:, None], X[tr_s], inner)
        Ahat[tr_s] = a_.ravel()
        Uhat[te_s], _ = ridge_cv(X[tr_s], U_tr, X[te_s], inner)
        a_, _ = ridge_cv(X[tr_s], alpha_tr[:, None], X[te_s], inner)
        Ahat[te_s] = a_.ravel()
        U_true = np.zeros((len(Y), rank))
        A_true = np.zeros(len(Y))
        U_true[tr_s], A_true[tr_s] = U_tr, alpha_tr

        # zero-shot protein maps, trained on training proteins only
        tr_pi, te_pi = np.where(tr_p)[0], np.where(te_p)[0]
        gtr = grp[tr_p]
        vpred = {}
        for nm, (kind, F) in PF.items():
            if kind == "kernel":
                Xtr, Xte = F[np.ix_(tr_pi, tr_pi)], F[np.ix_(te_pi, tr_pi)]
            else:
                Xtr, Xte = F[tr_pi], F[te_pi]
            best, be = 1.0, np.inf
            for lam in (0.1, 1.0, 10.0, 100.0, 1e3, 1e4):
                e = 0.0
                for c in sorted(set(gtr)):
                    msk = gtr != c
                    if msk.sum() < 4:
                        continue
                    e += float(((ridge_predict(Xtr[msk], V_tr[msk], Xtr[~msk], lam)
                                 - V_tr[~msk]) ** 2).sum())
                if e < be:
                    be, best = e, lam
            vpred[nm] = ridge_predict(Xtr, V_tr, Xte, best)

        for seed in SEEDS:
            rs = np.random.default_rng(10_000 * seed + f)
            for t, j in enumerate(te_pi):
                s_pool = np.where(M[:, j] & tr_s)[0]
                q_pool = np.where(M[:, j] & te_s)[0]
                if len(s_pool) < k + 5 or len(q_pool) < 5:
                    continue
                sup = rs.permutation(s_pool)[:k]
                resid = Y[:, j] - mu - A_true
                sol = solve_bv(resid, U_true, sup, LAM_V)
                b4, v4 = sol[0], sol[1:]
                b_int = float(resid[sup].mean())
                base = mu + Ahat[q_pool]
                yq = Y[q_pool, j]
                pred = {
                    "ADD": base + b_int,
                    "SUP": base + b4 + Uhat[q_pool] @ v4,
                    "PROT-ESM": base + b_int + Uhat[q_pool] @ vpred["PROT-ESM"][t],
                    "PROT-POCKET": base + b_int + Uhat[q_pool] @ vpred["PROT-POCKET"][t],
                    "PROT-GROUP": base + b_int + Uhat[q_pool] @ vpred["PROT-GROUP"][t],
                    "PROT-RANDOM": base + b_int + Uhat[q_pool] @ vpred["PROT-RANDOM"][t],
                    "PROT-ESM+SUP": base + b4 + Uhat[q_pool] @ (0.5 * (v4 + vpred["PROT-ESM"][t])),
                    "PROT-POCKET+SUP": base + b4 + Uhat[q_pool] @ (0.5 * (v4 + vpred["PROT-POCKET"][t])),
                }
                pc = str(grp[j])
                sc_q = scomp[q_pool]
                for arm, pr in pred.items():
                    e = yq - pr
                    s, n = err_p[arm].get(pc, (0.0, 0))
                    err_p[arm][pc] = (s + float((e ** 2).sum()), n + e.size)
                    for sc in set(sc_q):
                        sel = sc_q == sc
                        s, n = err_s[arm].get(sc, (0.0, 0))
                        err_s[arm][sc] = (s + float((e[sel] ** 2).sum()), n + int(sel.sum()))

    res = {"stage": "XP2-E", "closure": "double", "rank": rank, "k": k,
           "ligand_arm": ligand_arm, "arms": {}}
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
    res["contrasts"] = {
        f"Delta__ADD_minus_{a}": widest(
            paired_contrast(err_p[a], err_p["ADD"], n_boot, 21),
            paired_contrast(err_s[a], err_s["ADD"], n_boot, 22))
        for a in ARMS if a != "ADD" and err_p[a]}
    return res


if __name__ == "__main__":
    r = run()
    print(f"\n--- XP2-E biological landing (double closure, d={r['rank']}, k={r['k']}) ---")
    for a, v in r["arms"].items():
        g = v.get("r2_gamma_vs_ADD")
        gs = (f"{g['point']:+.4f} [{g['ci95'][0]:+.4f},{g['ci95'][1]:+.4f}]" if g
              else "(baseline)")
        print(f"  {a:18s} RMSE={v['rmse']:.4f}  R2_gamma vs ADD {gs}")
    p = os.path.join(OUT, "XP2E_BIOLOGICAL_LANDING.json")
    json.dump(r, open(p, "w"), indent=2, default=float)
    print("wrote", p)
