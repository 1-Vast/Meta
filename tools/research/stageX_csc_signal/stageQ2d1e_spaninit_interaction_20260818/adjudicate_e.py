"""Stage Q2d-1e gate adjudication (frozen gate, never moved).

Rules (prereg SHA baf4bb72...):
- DOUBLE-COLD primary surface; median over 3 truth seeds.
- PASS requires correct: dz >= 0.70 AND sp >= 0.30 AND gap(dz) vs
  ligand_only >= 0.05.
- Every negative arm must FAIL the gate.
- correct dz > best negative-arm dz + 0.03.
- Reported per surface: correct beats all controls (informative).
Writes Q2D1E_GATE.json. PASS requires the gate to hold at EVERY level
A-E of M1 (NC1/NC2 must fail; M2/M3 reported).
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NEGATIVE_ARMS = ["ligand_only", "additive_only", "shuffled_protein",
                 "family_preserving_shuffle", "random_protein",
                 "no_interaction_head"]


def main():
    ladder = json.loads((HERE / "Q2D1E_LADDER.json").read_text(encoding="utf-8"))
    results = ladder["results"]
    out = {"schema": "MetaSieve.StageQ2d1d.GATE.v1",
           "preregistration_sha256": ladder["preregistration_sha256"],
           "rules": {"dz_min": 0.70, "sp_min": 0.30, "gap_vs_ligand_only": 0.05,
                     "margin_over_best_negative": 0.03,
                     "surface": "dc", "median_over_seeds": [0, 1, 2]},
           "levels": {}}
    for mech in results:
        for level in results[mech]:
            per_seed = {}
            for seed in results[mech][level]:
                arms = results[mech][level][seed]
                per_seed[int(seed)] = {a: arms[a]["dc"] for a in arms}
            seeds = sorted(per_seed)
            def med(arm, key):
                return float(np.median([per_seed[s][arm][key] for s in seeds]))
            correct_dz = med("correct", "dz")
            correct_sp = med("correct", "sp")
            lig_dz = med("ligand_only", "dz")
            neg_dz = [med(a, "dz") for a in NEGATIVE_ARMS]
            best_neg = max(neg_dz)
            every_neg_fails = all(
                not (med(a, "dz") >= 0.70 and med(a, "sp") >= 0.30 and
                     med(a, "dz") - lig_dz >= 0.05) for a in NEGATIVE_ARMS)
            core_pass = (correct_dz >= 0.70 and correct_sp >= 0.30 and
                         correct_dz - lig_dz >= 0.05)
            margin_pass = correct_dz > best_neg + 0.03
            pass_ = core_pass and every_neg_fails and margin_pass
            out["levels"][f"{mech}:{level}"] = {
                "correct": {"dz": round(correct_dz, 4), "sp": round(correct_sp, 4)},
                "ligand_only_dz": round(lig_dz, 4),
                "negative_dz": {a: round(med(a, "dz"), 4) for a in NEGATIVE_ARMS},
                "best_negative_dz": round(best_neg, 4),
                "every_negative_fails": every_neg_fails,
                "core_pass": core_pass, "margin_pass": margin_pass,
                "pass": pass_,
            }
            print(mech, level, "PASS" if pass_ else "FAIL",
                  "correct dz/sp", round(correct_dz, 3), round(correct_sp, 3),
                  "best_neg", round(best_neg, 3), flush=True)
    gate_pass = all(out["levels"][f"M1:{lv}"]["pass"] for lv in ("A", "B", "C", "D", "E"))
    nc_fail = all(not out["levels"][f"{m}:A"]["pass"] for m in ("NC1", "NC2"))
    out["GATE_PASS"] = bool(gate_pass)
    out["NC_FAILS"] = bool(nc_fail)
    json.dump(out, open(HERE / "Q2D1E_GATE.json", "w"), indent=1)
    print("GATE:", "PASS" if gate_pass else "FAIL", "NC fails:", nc_fail)
    return 0 if gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
