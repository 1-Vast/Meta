"""Phase 1 — marginal decomposition of B4 and B5 from SEALED per-pair predictions.

Phase 1 requires residue-marginal AP, atom-marginal AP and exact-pair AP. This
reads the sealed float16 score memmaps, so nothing is retrained and every number
is recomputable from the sealed artifacts.

residue-marginal AP : rank residues by max_j s_ij, label d_i > 0
atom-marginal AP    : rank atoms by max_i s_ij,    label e_j > 0
exact-pair AP       : the complete-matrix AP already sealed
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import component_macro, paired_bootstrap  # noqa: E402
from p0_seal_predictions import ap_from_order  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
S4 = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "sealed_preds"
S5 = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "sealed_preds_b5"


def marg_ap(scores_2d, y_2d, axis):
    """axis=1 -> residue marginal; axis=0 -> atom marginal."""
    s = scores_2d.max(axis=1 if axis == 1 else 0)
    lab = (y_2d.sum(axis=1 if axis == 1 else 0) > 0).astype(np.int8)
    if lab.sum() == 0 or lab.sum() == lab.size:
        return None
    tie = np.arange(s.size)
    return ap_from_order(np.lexsort((tie, -s)), lab)


def main():
    kept, _q, _c, _f = build()
    comp_of = protein_components(kept)
    _tr, _ha, held_A, _hb = make_split(kept, comp_of)
    idx4 = json.loads((S4 / "heldoutA_index.json").read_text())
    held_A = [r for r in held_A if r["source_key"] in idx4]
    total4 = max(v[0] + v[1] * v[2] for v in idx4.values())

    arms = {}
    for name, d in (("B4", S4), ("BL", S4), ("B5", S5), ("BX5", S5), ("BP5", S5)):
        p = d / f"heldoutA_{name}.f16.dat"
        if p.exists():
            arms[name] = np.memmap(p, dtype=np.float16, mode="r", shape=(total4,))

    res_ap = {n: {} for n in arms}
    atom_ap = {n: {} for n in arms}
    for rec in held_A:
        k = rec["source_key"]
        off, L, A = idx4[k]
        y = np.zeros((L, A), dtype=np.int8)
        for i, j in rec["edges"]:
            y[i, j] = 1
        for n, mm in arms.items():
            s = np.asarray(mm[off:off + L * A], dtype=np.float32).reshape(L, A)
            res_ap[n][k] = marg_ap(s, y, 1)
            atom_ap[n][k] = marg_ap(s, y, 0)

    out = {"schema": "MetaSieve.S7L2B.P1.MarginalDecomposition.v1",
           "created_utc": "2026-08-10",
           "source": "sealed per-pair float16 predictions; nothing retrained",
           "complexes": len(held_A),
           "residue_marginal_ap": {}, "atom_marginal_ap": {}}
    comps = {}
    for n in arms:
        cr, mr = component_macro(res_ap[n], comp_of)
        ca, ma = component_macro(atom_ap[n], comp_of)
        out["residue_marginal_ap"][n] = mr
        out["atom_marginal_ap"][n] = ma
        comps[n] = (cr, ca)
        print(f"  {n:4s} residue-marginal AP={mr:.6f}   atom-marginal AP={ma:.6f}",
              flush=True)

    if "B5" in comps and "B4" in comps:
        out["residue_marginal_B5_minus_B4"] = paired_bootstrap(comps["B5"][0],
                                                               comps["B4"][0])
        out["atom_marginal_B5_minus_B4"] = paired_bootstrap(comps["B5"][1],
                                                            comps["B4"][1])
    out["oracle_reference"] = {"Oracle_R_residue_marginal_pair_AP": 0.216329,
                               "Oracle_A_atom_marginal_pair_AP": 0.008496,
                               "note": "the oracle figures are PAIR AP obtained by "
                                       "broadcasting a true marginal; the marginal APs "
                                       "above are ranked over residues/atoms and are a "
                                       "different quantity, not directly comparable"}
    (OUT / "P1_MARGINAL_DECOMPOSITION.json").write_text(json.dumps(out, indent=2),
                                                        encoding="utf-8")
    print(json.dumps({k: out[k] for k in
                      ("residue_marginal_ap", "atom_marginal_ap",
                       "residue_marginal_B5_minus_B4", "atom_marginal_B5_minus_B4")},
                     indent=2))


if __name__ == "__main__":
    main()
