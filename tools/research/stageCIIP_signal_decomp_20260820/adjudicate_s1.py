"""CIIP-S1 frozen adjudicator (S1_ADDENDUM_THRESHOLDS_20260820.md rules).

Reads SEED{seed}_RESULT.json / RESULT.json, applies the frozen rules, writes
ADJUDICATION.json. No thresholds are read from anywhere but this file's
frozen constants (mirroring the addendum, sha beefb620...).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

TAU_R2 = 0.05          # tau_min for delta centered R2 contrasts
TAU_SPEAR = 0.10       # tau_min for the T0m delta-Spearman contrast
MDE_R2 = 0.566         # MDE(80%) from S0 power table
MDE_SPEAR = 0.208
NULL_TOL = 0.03        # NULL-ALL tolerance above floor cell (R2 estimands)
NULL_TOL_SPEAR = 0.10  # NULL-ALL tolerance for T0m Spearman cells


def judge_contrast(c):
    if c.get("n", 0) < 3:
        return {"verdict": "UNDEFINED", "reason": "insufficient pairs"}
    lo, pt = c["lo2.5"], c["point"]
    stable = c["lopo_sign_stable"]
    tau = TAU_SPEAR if c.get("estimand") == "T0m" else TAU_R2
    mde = MDE_SPEAR if c.get("estimand") == "T0m" else MDE_R2
    if lo > 0 and stable and pt >= tau:
        label = "PASS (power-labeled: below MDE at n=9)" if pt < mde else "PASS"
        return {"verdict": label, "point": pt, "lo2.5": lo}
    if lo > 0 and not stable:
        return {"verdict": "UNRESOLVED (LOPO sign unstable)", "point": pt, "lo2.5": lo}
    if lo <= 0:
        return {"verdict": "UNRESOLVED (CI crosses 0)", "point": pt, "lo2.5": lo}
    return {"verdict": "UNRESOLVED (below tau_min)", "point": pt, "lo2.5": lo}


def adjudicate(res):
    ev = res["evals"]
    seed = res["seed"]
    out = {"seed": seed, "contrasts": {}}
    for c in res["contrasts"]:
        out["contrasts"][c["name"]] = judge_contrast(c)
    out["severity_contrast"] = judge_contrast(res["severity_contrast"])
    # NULL-ALL screen (compliance fix 2026-08-20): EVERY estimand cell vs its
    # F7f floor cell (+tol); T0m cells use cross-pair Spearman vs the F7f-T0m
    # floor (or 0 if that floor is undefined/constant).
    floors_ok, cells = True, []
    for est in ["T0", "T1", "T2", "T3"]:
        f7 = ev.get(f"F7f-{est}-s{seed}", {}).get("agg", {}).get("mean_cr2")
        for tag, e in ev.items():
            if tag.startswith("F9") or f"-{est}-" not in tag:
                continue
            a = e.get("agg", {}).get("mean_cr2")
            cells.append({"arm": tag, "estimand": est, "cr2": a, "floor": f7})
            if a is not None and f7 is not None and a > f7 + NULL_TOL:
                floors_ok = False
    f7m = ev.get(f"F7f-T0m-s{seed}", {}).get("agg", {}).get("spearman")
    f7m_floor = f7m if f7m is not None and f7m == f7m else 0.0
    for tag, e in ev.items():
        if tag.startswith("F9") or "-T0m-" not in tag:
            continue
        a = e.get("agg", {}).get("spearman")
        cells.append({"arm": tag, "estimand": "T0m", "spearman": a,
                      "floor": f7m_floor})
        if a is not None and a == a and a > f7m_floor + NULL_TOL_SPEAR:
            floors_ok = False
    # permutation behavior: destroying labels must not HELP any arm
    perm_ok = True
    for tag, pe in res["perm_evals"].items():
        base = ev.get(tag, {}).get("agg", {}).get("mean_cr2")
        pm = pe.get("agg", {}).get("mean_cr2")
        if base is not None and pm is not None and pm > base + NULL_TOL:
            perm_ok = False
    any_pass = any(v["verdict"].startswith("PASS") for v in
                   list(out["contrasts"].values()) + [out["severity_contrast"]])
    out["null_all_screen"] = {"all_cells_at_floor_plus_tol": floors_ok,
                              "permutation_no_gain": perm_ok,
                              "cells": cells}
    if any_pass:
        out["overall"] = "PASS (see per-contrast verdicts)"
    elif floors_ok and perm_ok:
        out["overall"] = "S1-NULL-ALL"
    else:
        out["overall"] = "S1-UNRESOLVED (CI crosses 0; NULL-ALL conditions not met)"
    return out


def aggregate_seeds(adj_by_seed):
    """Seed stability: sign consistency of each contrast point across seeds."""
    names = set()
    for a in adj_by_seed.values():
        names.update(a["contrasts"])
        names.add("severity_contrast")
    agg = {}
    for n in sorted(names):
        pts, los, verdicts = [], [], []
        for s, a in sorted(adj_by_seed.items()):
            c = a["severity_contrast"] if n == "severity_contrast" else a["contrasts"].get(n)
            if c and "point" in c:
                pts.append(c["point"]); los.append(c["lo2.5"])
                verdicts.append(c["verdict"])
        if not pts:
            continue
        same_sign = all(np.sign(p) == np.sign(pts[0]) for p in pts)
        agg[n] = {"mean_point": float(np.mean(pts)), "points": pts,
                  "lo2.5_by_seed": los, "seed_sign_consistent": bool(same_sign),
                  "verdicts": verdicts}
    return agg


def main():
    import sys
    src = sys.argv[1] if len(sys.argv) > 1 else "RESULT.json"
    data = json.loads((HERE / src).read_text(encoding="utf-8"))
    if "seed" in data:
        out = adjudicate(data)
    else:
        out = {k: adjudicate(v) for k, v in data.items()}
        out["seed_aggregate"] = aggregate_seeds(out)
    (HERE / "ADJUDICATION.json").write_text(json.dumps(out, indent=1, default=str),
                                            encoding="utf-8")
    print(json.dumps(out if "seed" in data else
                     {"overall": {k: v["overall"] for k, v in out.items()
                                  if isinstance(v, dict) and "overall" in v},
                      "aggregate": out.get("seed_aggregate", {})}, indent=1, default=str))


if __name__ == "__main__":
    main()
