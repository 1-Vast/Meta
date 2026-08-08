"""Candidate mechanism: kinase conformational-state availability (DFG / aC-helix)
and subpocket accessibility, from KLIFS experimental structures.

Mechanistic motivation: type-II and allosteric binding requires a DFG-out /
aC-out accessible state.  This is a property of the protein's conformational
ensemble, not of its sequence at the aligned pocket positions, so it is a
candidate for exactly the information a sequence encoder cannot recover.

CAVEAT recorded with the result: KLIFS conformation labels are derived from
ligand-bound crystal structures, so the *availability* of a DFG-out structure is
partly a consequence of medicinal-chemistry effort on that kinase.  This arm is
therefore reported as a mechanism probe with a stated circularity risk, not as a
clean prospective feature.
"""
from __future__ import annotations

import json
import os
import ssl
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panels import SYMBOL_ALIASES, load_metz  # noqa: E402

ANN = r"D:\MetaSieve\dataset\raw\crossed_panels\protein_annotation"
CACHE = r"D:\MetaSieve\dataset\processed\crossed_panels"
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

SUBPOCKETS = ["front", "gate", "back", "fp_I", "fp_II", "bp_I_A", "bp_I_B",
              "bp_II_in", "bp_II_A_in", "bp_II_B_in", "bp_II_out", "bp_II_B",
              "bp_III", "bp_IV", "bp_V"]


def fetch(kid):
    url = f"https://klifs.net/api/structures_list?kinase_ID={kid}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    return json.loads(urllib.request.urlopen(req, timeout=90, context=CTX).read())


def main(density=0.60):
    rec = json.load(open(os.path.join(ANN, "klifs_kinase_information_human.json"),
                         encoding="utf-8"))
    by_hgnc = {(r.get("HGNC") or "").upper(): r for r in rec}
    by_name = {(r.get("name") or "").upper(): r for r in rec}
    _, _, _, kin = load_metz(density)

    cache_path = os.path.join(CACHE, "klifs_structures.json")
    store = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    rows, names = [], []
    for k in kin:
        s = SYMBOL_ALIASES.get(k, k)
        r = by_hgnc.get(s) or by_name.get(s)
        kid = str(r["kinase_ID"])
        if kid not in store:
            try:
                store[kid] = fetch(kid)
            except Exception as e:
                print(f"  !! {k}/{kid}: {type(e).__name__}")
                store[kid] = []
        st = [x for x in store[kid] if isinstance(x, dict)]
        n = len(st)
        def frac(f):
            return float(np.mean([bool(f(x)) for x in st])) if n else 0.0
        feat = [
            np.log1p(n),
            frac(lambda x: x.get("DFG") == "out"),
            frac(lambda x: x.get("DFG") == "in"),
            frac(lambda x: x.get("aC_helix") == "out"),
            frac(lambda x: x.get("aC_helix") == "in"),
            frac(lambda x: x.get("DFG") == "out" and x.get("aC_helix") == "out"),
            frac(lambda x: x.get("allosteric_ligand") not in (0, "0", None)),
            frac(lambda x: x.get("ligand") in (0, "0", None)),
        ] + [frac(lambda x, s=s_: bool(x.get(s))) for s_ in SUBPOCKETS]
        rows.append(feat)
        names.append(k)
    json.dump(store, open(cache_path, "w"))
    X = np.array(rows, float)
    cols = (["log_n_structures", "frac_DFG_out", "frac_DFG_in", "frac_aC_out",
             "frac_aC_in", "frac_DFGout_aCout", "frac_allosteric", "frac_apo"]
            + [f"frac_{s}" for s in SUBPOCKETS])
    out = os.path.join(CACHE, f"metz{int(density*100)}_conformation_features.npz")
    np.savez(out, kinases=np.array(names), X=X, columns=np.array(cols))
    print("wrote", out, X.shape)
    print("kinases with zero structures:",
          int((X[:, 0] == 0).sum()), "/", len(names))
    for i, c in enumerate(cols[:8]):
        print(f"   {c:22s} mean {X[:, i].mean():.3f}  sd {X[:, i].std():.3f}")


if __name__ == "__main__":
    main(0.60)
    main(0.70)
