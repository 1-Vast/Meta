"""XP1-C: independent-family replication of XP1-B on the PDSP Ki panel.

Different protein class (GPCR / transporter / channel), different assay
technology (radioligand displacement), different laboratory population, and a
sparse rather than complete crossing.  Same arm structure as XP1-B.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lowrank import fit_interaction_basis  # noqa: E402
from xp1b_transfer import (REPORT, _boot_indices, contrast, r2_gamma,  # noqa: E402
                           ridge_predict, solve_bv)

CACHE = r"D:\MetaSieve\dataset\processed\crossed_panels"


def esm_for_pdsp(targets, seqs):
    import torch
    from transformers import AutoTokenizer, EsmModel
    path = os.path.join(CACHE, "esm2_t30_pdsp.npz")
    if os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        m = {k: v for k, v in zip(d["keys"], d["emb"])}
    else:
        name = "facebook/esm2_t30_150M_UR50D"
        tok = AutoTokenizer.from_pretrained(name)
        mdl = EsmModel.from_pretrained(name).eval()
        dev = "cuda" if torch.cuda.is_available() else "cpu"
        mdl = mdl.to(dev)
        keys, embs = [], []
        with torch.no_grad():
            for i in range(0, len(targets), 4):
                ch = [str(t) for t in targets[i:i + 4]]
                enc = tok([seqs[c][:1024] for c in ch], return_tensors="pt",
                          padding=True, truncation=True, max_length=1024).to(dev)
                out = mdl(**enc).last_hidden_state
                mk = enc["attention_mask"].unsqueeze(-1).float()
                pooled = ((out * mk).sum(1) / mk.sum(1)).cpu().numpy()
                keys += ch
                embs += list(pooled)
        m = {k: v for k, v in zip(keys, embs)}
        np.savez(path, keys=np.array(keys), emb=np.stack(embs))
    return np.stack([m[str(t)] for t in targets])


def kmer_kernel(seqs, targets, k=3):
    from collections import Counter
    vocab = {}
    vecs = []
    for t in targets:
        s = seqs[str(t)]
        c = Counter(s[i:i + k] for i in range(len(s) - k + 1))
        for key in c:
            vocab.setdefault(key, len(vocab))
        vecs.append(c)
    X = np.zeros((len(targets), len(vocab)))
    for i, c in enumerate(vecs):
        for key, v in c.items():
            X[i, vocab[key]] = v
    X = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-9)
    return X @ X.T


def family_prefix(targets):
    out = []
    for t in map(str, targets):
        pre = "".join(ch for ch in t if ch.isalpha())
        out.append(pre[:4] if len(pre) > 4 else pre)
    return out


def onehot(labels):
    lv = sorted(set(labels))
    X = np.zeros((len(labels), len(lv)))
    for i, v in enumerate(labels):
        X[i, lv.index(v)] = 1.0
    return X


def run(rank=6, k_support=16, seeds=(0, 1, 2, 3, 4), lam_lr=5.0, lam_v=2.0,
        min_train_obs=3, n_boot=2000):
    d = np.load(os.path.join(CACHE, "pdsp_core.npz"), allow_pickle=True)
    Y, M = d["Y"], d["M"]
    targets, comp = d["targets"], d["cluster"]
    seqs = json.load(open(os.path.join(CACHE, "pdsp_sequences.json")))
    fam = family_prefix(targets)
    Kk = kmer_kernel(seqs, targets)
    E = esm_for_pdsp(targets, seqs)
    E = (E - E.mean(0)) / (E.std(0) + 1e-9)
    FEAT = {"esm2_t30_fullseq": E, "kmer3_kernel": Kk, "family_prefix": onehot(fam)}

    comps_all = sorted(set(comp))
    order = sorted(comps_all, key=lambda c: -int((comp == c).sum()))
    folds = {c: i % 5 for i, c in enumerate(order)}
    print(f"PDSP: {Y.shape[0]} ligands x {Y.shape[1]} targets, {int(M.sum())} cells, "
          f"{len(comps_all)} homology clusters (mmseqs 40%)")

    arms = ["A0", "A1", "A2", "A4", "A6", "AO1"] + \
           [f"A3::{f}" for f in FEAT] + [f"A5::{f}" for f in FEAT]
    err = {a: {} for a in arms}

    def add(arm, c, e):
        s, n = err[arm].get(c, (0.0, 0))
        err[arm][c] = (s + float((e ** 2).sum()), n + len(e))

    for f in range(5):
        te_comp = [c for c in comps_all if folds[c] == f]
        te_idx = np.where(np.isin(comp, te_comp))[0]
        tr_idx = np.where(~np.isin(comp, te_comp))[0]
        if len(te_idx) == 0 or len(tr_idx) < 8:
            continue
        Mtr = M[:, tr_idx]
        m = fit_interaction_basis(Y[:, tr_idx], Mtr, rank, lam=lam_lr, seed=100 + f)
        mu, alpha, U, Vtr = m["mu"], m["alpha"], m["U"], m["V"]
        good_lig = Mtr.sum(axis=1) >= min_train_obs      # U_i actually estimated

        vpred = {}
        inner = comp[tr_idx]
        for fname, X in FEAT.items():
            if fname.endswith("kernel"):
                Xtr, Xte = X[np.ix_(tr_idx, tr_idx)], X[np.ix_(te_idx, tr_idx)]
            else:
                Xtr, Xte = X[tr_idx], X[te_idx]
            best, best_e = 1.0, np.inf
            for lam in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0):
                e = 0.0
                for c in sorted(set(inner)):
                    a = inner != c
                    if a.sum() < 4:
                        continue
                    P = ridge_predict(Xtr[a], Vtr[a], Xtr[~a], lam)
                    e += float(((P - Vtr[~a]) ** 2).sum())
                if e < best_e:
                    best_e, best = e, lam
            vpred[fname] = ridge_predict(Xtr, Vtr, Xte, best)

        for seed in seeds:
            rs = np.random.default_rng(1000 * seed + f)
            nT = len(te_idx)
            der = rs.permutation(nT)
            tr_ = 0
            while nT > 1 and np.any(der == np.arange(nT)) and tr_ < 500:
                der = rs.permutation(nT); tr_ += 1
            sup, tst, bv, b2 = {}, {}, {}, {}
            for t, j in enumerate(te_idx):
                obs = np.where(M[:, j] & good_lig)[0]
                if len(obs) < k_support + 20:
                    continue
                pick = rs.permutation(obs)
                sup[t], tst[t] = pick[:k_support], pick[k_support:]
                resid = Y[:, j] - mu - alpha
                bv[t] = solve_bv(resid, U, sup[t], lam_v)
                b2[t] = float(resid[sup[t]].mean())
            valid = [t for t in tst if der[t] in tst]
            for t in valid:
                j, ti = te_idx[t], tst[t]
                yte, c = Y[ti, j], comp[j]
                base = mu + alpha[ti]
                add("A0", c, yte - mu)
                add("A1", c, yte - base)
                a2 = base + b2[t]
                add("A2", c, yte - a2)
                for fn in FEAT:
                    add(f"A3::{fn}", c, yte - (a2 + U[ti] @ vpred[fn][t]))
                    add(f"A5::{fn}", c, yte - (a2 + U[ti] @ vpred[fn][der[t]]))
                b4, v4 = bv[t]
                add("A4", c, yte - (base + b4 + U[ti] @ v4))
                _, v6 = bv[der[t]]
                add("A6", c, yte - (base + b4 + U[ti] @ v6))
                bo, vo = solve_bv(Y[:, j] - mu - alpha, U,
                                  np.where(M[:, j] & good_lig)[0], lam_v)
                add("AO1", c, yte - (base + bo + U[ti] @ vo))

    comps = sorted(err["A2"])
    bidx = _boot_indices(len(comps), n_boot, seed=99)
    out = {"panel": "BLK-PDSP-H", "rank": rank, "k_support": k_support,
           "n_eval_components": len(comps), "arms": {}}
    for a in arms:
        if not err[a]:
            continue
        s = sum(v[0] for v in err[a].values())
        n = sum(v[1] for v in err[a].values())
        out["arms"][a] = {"rmse": float((s / n) ** 0.5), "cells": int(n)}
        if a != "A2":
            out["arms"][a]["r2_gamma_vs_A2"] = r2_gamma(err[a], err["A2"], comps, bidx)
    out["contrasts"] = {
        "Delta_protein__A1_minus_A4": contrast(err["A4"], err["A1"], comps, bidx),
        "Delta_interaction__A2_minus_A4": contrast(err["A4"], err["A2"], comps, bidx),
        "Delta_specific__A6_minus_A4": contrast(err["A4"], err["A6"], comps, bidx),
        "Delta_interaction__A2_minus_A3esm": contrast(err["A3::esm2_t30_fullseq"],
                                                      err["A2"], comps, bidx),
        "Delta_specific__A5esm_minus_A3esm": contrast(err["A3::esm2_t30_fullseq"],
                                                      err["A5::esm2_t30_fullseq"],
                                                      comps, bidx),
        "Delta_oracle__A2_minus_AO1": contrast(err["AO1"], err["A2"], comps, bidx),
    }
    return out


if __name__ == "__main__":
    res = run()
    print(f"\n{'arm':28s} {'RMSE':>8s}  {'R2_gamma vs A2 [95% CI]':>30s}")
    for a, v in res["arms"].items():
        r = v.get("r2_gamma_vs_A2")
        rs = (f"{r['point']:+.4f} [{r['ci95'][0]:+.4f},{r['ci95'][1]:+.4f}]"
              if r else "(baseline)")
        print(f"{a:28s} {v['rmse']:8.4f}  {rs:>30s}")
    print("  contrasts:")
    for k, v in res["contrasts"].items():
        print(f"   {k:36s} {v['point']:+.5f} [{v['ci95'][0]:+.5f},{v['ci95'][1]:+.5f}]"
              f"   RMSE {v['rmse_b']:.4f} -> {v['rmse_a']:.4f}")
    p = os.path.join(REPORT, "xp1c_pdsp.json")
    json.dump(res, open(p, "w"), indent=2, default=float)
    print("wrote", p)
