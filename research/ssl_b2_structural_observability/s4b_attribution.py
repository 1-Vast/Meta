"""S4b — attribution control: is the observable teacher signal protein-side or ligand-side?

S4 showed two channels are predictable from ESM+ECFP and beat random features,
but no channel beat the deranged-protein control.  This script settles the
attribution directly by adding ligand-only and protein-only arms on the SAME
cells and split.  It introduces no new hypothesis and no new hyperparameter
search: the ridge lambda selection is identical to S4.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "crossed_panel_deployability"))
from s2_teacher import CHANNELS  # noqa: E402
from s3s4_observability import CACHE, LAMS, N_FOLDS, OUT, build_dataset, esm  # noqa: E402
from xp2_core import paired_contrast, r2_vs_base, ridge_predict, widest  # noqa: E402


def main(n_boot=2000):
    d = build_dataset()
    Y, L, pclus, scaf = d["Y"], d["L"], d["pclus"], d["scaffold"]
    P = esm(d["seqs"])
    Ys = (Y - Y.mean(0)) / (Y.std(0) + 1e-9)
    rng = np.random.default_rng(0)
    pcs, scs = sorted(set(pclus)), sorted(set(scaf))
    o1, o2 = list(pcs), list(scs)
    rng.shuffle(o1)
    rng.shuffle(o2)
    pf = {c: i % N_FOLDS for i, c in enumerate(o1)}
    sf = {c: i % N_FOLDS for i, c in enumerate(o2)}

    ARMS = ["B0", "LIG-ONLY", "PROT-ONLY", "BOTH", "RAND"]
    err = {a: {c: {} for c in CHANNELS} for a in ARMS}
    errs = {a: {c: {} for c in CHANNELS} for a in ARMS}

    for f in range(N_FOLDS):
        te = np.array([pf[a] == f and sf[b] == f for a, b in zip(pclus, scaf)])
        tr = np.array([pf[a] != f and sf[b] != f for a, b in zip(pclus, scaf)])
        if te.sum() < 5 or tr.sum() < 50:
            continue
        Pz = (P - P[tr].mean(0)) / (P[tr].std(0) + 1e-9)
        Lz = (L - L[tr].mean(0)) / (L[tr].std(0) + 1e-9)
        feats = {"LIG-ONLY": Lz, "PROT-ONLY": Pz, "BOTH": np.hstack([Pz, Lz])}
        feats["RAND"] = rng.normal(size=feats["BOTH"].shape)
        inner = pclus[tr]
        pred = {"B0": np.repeat(Ys[tr].mean(0)[None, :], te.sum(), axis=0)}
        for a, X in feats.items():
            best, be = LAMS[0], np.inf
            for lam in LAMS:
                e = 0.0
                for c in sorted(set(inner))[:40]:
                    m = inner != c
                    if m.sum() < 30 or (~m).sum() == 0:
                        continue
                    e += float(((ridge_predict(X[tr][m], Ys[tr][m], X[tr][~m], lam)
                                 - Ys[tr][~m]) ** 2).sum())
                if e < be:
                    be, best = e, lam
            pred[a] = ridge_predict(X[tr], Ys[tr], X[te], best)
        for a, pr in pred.items():
            for ci, c in enumerate(CHANNELS):
                e = Ys[te][:, ci] - pr[:, ci]
                for u in set(pclus[te]):
                    m = pclus[te] == u
                    s, k = err[a][c].get(u, (0.0, 0))
                    err[a][c][u] = (s + float((e[m] ** 2).sum()), k + int(m.sum()))
                for u in set(scaf[te]):
                    m = scaf[te] == u
                    s, k = errs[a][c].get(u, (0.0, 0))
                    errs[a][c][u] = (s + float((e[m] ** 2).sum()), k + int(m.sum()))

    res = {"stage": "S4b", "question": "is the observable teacher signal ligand-side "
                                       "or protein-side?", "channels": {}}
    print(f"{'channel':24s} {'LIG-ONLY':>22s} {'PROT-ONLY':>22s} {'BOTH':>22s}"
          f" {'BOTH-minus-LIG':>24s}")
    for c in CHANNELS:
        if not err["B0"][c]:
            continue
        r = {}
        for a in ("LIG-ONLY", "PROT-ONLY", "BOTH"):
            r[a] = widest(r2_vs_base(err[a][c], err["B0"][c], n_boot, 11),
                          r2_vs_base(errs[a][c], errs["B0"][c], n_boot, 12))
        gain = widest(paired_contrast(err["BOTH"][c], err["LIG-ONLY"][c], n_boot, 21),
                      paired_contrast(errs["BOTH"][c], errs["LIG-ONLY"][c], n_boot, 22))
        res["channels"][c] = {**{k: v for k, v in r.items()},
                              "protein_gain_over_ligand_only": gain}
        print(f"{c:24s} {r['LIG-ONLY']['point']:+.3f}"
              f"[{r['LIG-ONLY']['ci95'][0]:+.3f},{r['LIG-ONLY']['ci95'][1]:+.3f}]"
              f" {r['PROT-ONLY']['point']:+.3f}"
              f"[{r['PROT-ONLY']['ci95'][0]:+.3f},{r['PROT-ONLY']['ci95'][1]:+.3f}]"
              f" {r['BOTH']['point']:+.3f}"
              f"[{r['BOTH']['ci95'][0]:+.3f},{r['BOTH']['ci95'][1]:+.3f}]"
              f" {gain['point']:+.4f}[{gain['ci95'][0]:+.4f},{gain['ci95'][1]:+.4f}]")
    json.dump(res, open(os.path.join(OUT, "S4B_ATTRIBUTION.json"), "w"), indent=2,
              default=float)
    print("\nwrote S4B_ATTRIBUTION.json")
    return res


if __name__ == "__main__":
    main()
