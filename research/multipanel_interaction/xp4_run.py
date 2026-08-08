"""XP4 — does a deployment-observable basis predict the within-panel interaction?

Estimand: gamma, the within-panel double-centred interaction residual.  Because
main effects are removed inside each panel, panel assay offset, panel potency
scale and target druggability cannot be exploited by any arm.

Every arm predicts the same held-out cells under the double (panel-cluster x
scaffold-component) split.  Predictions are double-centred within panel before
scoring, so no arm can win by reintroducing a main effect.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "crossed_panel_deployability"))
from xp2_core import paired_contrast, r2_vs_base, widest  # noqa: E402
from xp4_build import CACHE, build  # noqa: E402

OUT = r"D:\MetaSieve\report\multipanel_interaction"
os.makedirs(OUT, exist_ok=True)
N_FOLDS = 5
PCA_DIM = 32
RANKS = (1, 2, 4, 8)
LAMS = (1.0, 10.0, 100.0, 1e3, 1e4, 1e5)
SEED = 0


def esm_embed(uniprot, seqs):
    path = os.path.join(CACHE, "esm2_t30_bdb.npz")
    if os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        m = {k: v for k, v in zip(z["keys"], z["emb"])}
        if all(u in m for u in uniprot):
            return np.stack([m[u] for u in uniprot])
    import torch
    from transformers import AutoTokenizer, EsmModel
    name = "facebook/esm2_t30_150M_UR50D"
    tok = AutoTokenizer.from_pretrained(name)
    mdl = EsmModel.from_pretrained(name).eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    mdl = mdl.to(dev)
    uq = sorted(set(uniprot))
    keys, embs = [], []
    with torch.no_grad():
        for i in range(0, len(uq), 4):
            ch = uq[i:i + 4]
            enc = tok([seqs[c][:1022] for c in ch], return_tensors="pt", padding=True,
                      truncation=True, max_length=1022).to(dev)
            o = mdl(**enc).last_hidden_state
            mk = enc["attention_mask"].unsqueeze(-1).float()
            for c, e in zip(ch, ((o * mk).sum(1) / mk.sum(1)).cpu().numpy()):
                keys.append(c)
                embs.append(e)
            print(f"  esm {min(i+4, len(uq))}/{len(uq)}", end="\r")
    print()
    np.savez(path, keys=np.array(keys), emb=np.stack(embs))
    m = {k: v for k, v in zip(keys, embs)}
    return np.stack([m[u] for u in uniprot])


def ecfp(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    uq = sorted(set(smiles))
    m = {s: np.array(list(gen.GetFingerprint(Chem.MolFromSmiles(s))), float) for s in uq}
    return np.stack([m[s] for s in smiles])


def within_panel_center(vals, pid, lid, tid, iters=200):
    """Double-centre inside each panel: remove panel mean, ligand and target effects."""
    out = np.zeros_like(vals)
    for p in np.unique(pid):
        m = pid == p
        y = vals[m].copy()
        li, ti = lid[m], tid[m]
        mu = y.mean()
        a = np.zeros(li.max() + 1)
        b = np.zeros(ti.max() + 1)
        for _ in range(iters):
            r = y - mu - b[ti]
            for u in np.unique(li):
                a[u] = r[li == u].mean()
            a -= a.mean()
            r = y - mu - a[li]
            for u in np.unique(ti):
                b[u] = r[ti == u].mean()
            b -= b.mean()
            mu = (y - a[li] - b[ti]).mean()
        out[m] = y - mu - a[li] - b[ti]
    return out


def pca_fit(X, dim, seed=0):
    mu = X.mean(0)
    Xc = X - mu
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return mu, Vt[:dim].T


def bilinear_fit(P, L, g, lam, rank):
    """Ridge on the outer-product features, then SVD truncation to `rank`."""
    n, dp = P.shape
    dl = L.shape[1]
    Z = (P[:, :, None] * L[:, None, :]).reshape(n, dp * dl)
    A = Z.T @ Z + lam * np.eye(dp * dl)
    w = np.linalg.solve(A, Z.T @ g)
    W = w.reshape(dp, dl)
    U, S, Vt = np.linalg.svd(W, full_matrices=False)
    r = min(rank, len(S))
    return (U[:, :r] * S[:r]) @ Vt[:r]


def predict(P, L, W):
    return np.einsum("ij,jk,ik->i", P, W, L)


def run(n_boot=2000):
    d, man = build()
    pmid, uni, smi = d["pmid"], d["uni"], d["smiles"]
    y = d["pki"].astype(float)
    scomp, pclus = d["scaffold_component"], d["protein_cluster"]
    seqs = json.load(open(os.path.join(CACHE, "bdb_sequences.json")))
    print(f"panels={man['panels']} cells={man['cells']} targets={man['targets']} "
          f"ligands={man['ligands']} scaffold_components={man['scaffold_components']} "
          f"protein_clusters={man['protein_clusters_mmseqs40']}")

    lig_ids = {s: i for i, s in enumerate(sorted(set(smi)))}
    tgt_ids = {t: i for i, t in enumerate(sorted(set(uni)))}
    lid = np.array([lig_ids[s] for s in smi])
    tid = np.array([tgt_ids[t] for t in uni])
    gamma = within_panel_center(y, pmid, lid, tid)
    print(f"gamma sd = {gamma.std(ddof=1):.4f} log units, raw pKi sd = {y.std(ddof=1):.4f}")

    Pfull = esm_embed(list(uni), seqs)
    Lfull = ecfp(list(smi))
    rng = np.random.default_rng(SEED)

    pclusters = sorted(set(pclus))
    scomps = sorted(set(scomp))
    o1 = list(pclusters); rng.shuffle(o1)
    o2 = list(scomps); rng.shuffle(o2)
    pfold = {c: i % N_FOLDS for i, c in enumerate(o1)}
    sfold = {c: i % N_FOLDS for i, c in enumerate(o2)}
    # a panel is held out if ANY of its targets is in a held-out protein cluster
    panel_fold = {}
    for p in sorted(set(pmid)):
        m = pmid == p
        panel_fold[p] = min(pfold[c] for c in set(pclus[m]))

    ARMS = ["Z0", "BILIN", "RAND-P", "RAND-L", "RAND-BOTH", "PERM-PAIR",
            "FOREIGN-P", "ORACLE-R"]
    err_p = {a: {} for a in ARMS}
    err_s = {a: {} for a in ARMS}
    chosen = []

    for f in range(N_FOLDS):
        te_panel = np.array([panel_fold[p] == f for p in pmid])
        te_sc = np.array([sfold[c] == f for c in scomp])
        te = te_panel & te_sc
        tr = (~te_panel) & (~te_sc)
        if te.sum() < 20 or tr.sum() < 200:
            continue
        # frozen encoders reduced by PCA fitted on TRAINING cells only
        mup, Vp = pca_fit(Pfull[tr], PCA_DIM)
        mul, Vl = pca_fit(Lfull[tr], PCA_DIM)
        Pp = (Pfull - mup) @ Vp
        Ll = (Lfull - mul) @ Vl
        Pp /= (Pp[tr].std(0) + 1e-9)
        Ll /= (Ll[tr].std(0) + 1e-9)
        RP = rng.normal(size=(len(pmid), PCA_DIM))
        RL = rng.normal(size=(len(pmid), PCA_DIM))

        # nested panel-grouped CV inside the training cells for (lam, rank)
        inner_panels = sorted(set(pmid[tr]))
        rng.shuffle(inner_panels)
        ib = {p: i % 4 for i, p in enumerate(inner_panels)}
        ibin = np.array([ib.get(p, -1) for p in pmid])
        best, be = (LAMS[0], RANKS[0]), np.inf
        for lam in LAMS:
            for rk in RANKS:
                e = 0.0
                for b in range(4):
                    a = tr & (ibin != b)
                    h = tr & (ibin == b)
                    if a.sum() < 100 or h.sum() < 20:
                        continue
                    W = bilinear_fit(Pp[a], Ll[a], gamma[a], lam, rk)
                    e += float(((predict(Pp[h], Ll[h], W) - gamma[h]) ** 2).sum())
                if e < be:
                    be, best = e, (lam, rk)
        lam, rk = best
        chosen.append({"fold": f, "lambda": lam, "rank": rk,
                       "train_cells": int(tr.sum()), "test_cells": int(te.sum())})

        W = bilinear_fit(Pp[tr], Ll[tr], gamma[tr], lam, rk)
        W_rp = bilinear_fit(RP[tr], Ll[tr], gamma[tr], lam, rk)
        W_rl = bilinear_fit(Pp[tr], RL[tr], gamma[tr], lam, rk)
        W_rb = bilinear_fit(RP[tr], RL[tr], gamma[tr], lam, rk)
        perm = rng.permutation(np.where(tr)[0])
        Pperm = Pp.copy()
        Pperm[np.where(tr)[0]] = Pp[perm]
        W_pp = bilinear_fit(Pperm[tr], Ll[tr], gamma[tr], lam, rk)

        # foreign protein at test time: rotate protein features among test panels
        te_idx = np.where(te)[0]
        tp = sorted(set(pmid[te]))
        shift = {p: tp[(i + 1) % len(tp)] for i, p in enumerate(tp)}
        Pfor = Pp.copy()
        for p in tp:
            src = np.where((pmid == shift[p]) & te)[0]
            dst = np.where((pmid == p) & te)[0]
            if len(src):
                Pfor[dst] = Pp[src[rng.integers(0, len(src), len(dst))]]

        # oracle: the held-out panel's own training-scaffold cells, same basis
        orc = np.zeros(len(pmid))
        for p in tp:
            own_tr = (pmid == p) & (~te_sc)
            own_te = (pmid == p) & te
            if own_tr.sum() >= 10:
                Wo = bilinear_fit(Pp[own_tr], Ll[own_tr], gamma[own_tr], lam, rk)
                orc[own_te] = predict(Pp[own_te], Ll[own_te], Wo)

        pred = {
            "Z0": np.zeros(len(te_idx)),
            "BILIN": predict(Pp[te], Ll[te], W),
            "RAND-P": predict(RP[te], Ll[te], W_rp),
            "RAND-L": predict(Pp[te], RL[te], W_rl),
            "RAND-BOTH": predict(RP[te], RL[te], W_rb),
            "PERM-PAIR": predict(Pp[te], Ll[te], W_pp),
            "FOREIGN-P": predict(Pfor[te], Ll[te], W),
            "ORACLE-R": orc[te],
        }
        # double-centre every prediction within panel so no arm can smuggle a main effect
        for a in pred:
            pred[a] = within_panel_center(pred[a], pmid[te], lid[te], tid[te], iters=50)

        gt = gamma[te]
        for a, pr in pred.items():
            e = gt - pr
            for p in set(pmid[te]):
                m = pmid[te] == p
                s, n = err_p[a].get(p, (0.0, 0))
                err_p[a][p] = (s + float((e[m] ** 2).sum()), n + int(m.sum()))
            for c in set(scomp[te]):
                m = scomp[te] == c
                s, n = err_s[a].get(c, (0.0, 0))
                err_s[a][c] = (s + float((e[m] ** 2).sum()), n + int(m.sum()))

    res = {"stage": "XP4", "panel": man, "hyperparameters_selected": chosen,
           "gamma_sd": float(gamma.std(ddof=1)), "arms": {}, "contrasts": {}}
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
        "Delta_null__Z0_minus_BILIN": ("BILIN", "Z0"),
        "Delta_bio__RANDBOTH_minus_BILIN": ("BILIN", "RAND-BOTH"),
        "Delta_protein__RANDP_minus_BILIN": ("BILIN", "RAND-P"),
        "Delta_ligand__RANDL_minus_BILIN": ("BILIN", "RAND-L"),
        "Delta_pair__PERMPAIR_minus_BILIN": ("BILIN", "PERM-PAIR"),
        "Delta_specific__FOREIGNP_minus_BILIN": ("BILIN", "FOREIGN-P"),
        "Delta_oracle__Z0_minus_ORACLE": ("ORACLE-R", "Z0"),
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
    p = os.path.join(OUT, "XP4_MULTIPANEL_RESULT.json")
    json.dump(r, open(p, "w"), indent=2, default=float)
    print("wrote", p)
