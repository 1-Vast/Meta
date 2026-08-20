"""Reachability table: each arm's gap to the F9 parent-profile ceiling on T1
(deliverable: parent-profile ceiling reachability). Mean + median views."""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
seeds = sys.argv[1:] if len(sys.argv) > 1 else ["1", "2", "3"]
rows = {}
for s in seeds:
    r = json.loads((HERE / f"SEED{s}_RESULT.json").read_text(encoding="utf-8"))
    ev = r["evals"]
    seed_rows = {}
    for tag, e in ev.items():
        if e.get("estimand") == "T0m" or "-T1-" not in tag:
            continue
        vals = [v["cr2"] for v in e["per_pair"].values()
                if not v.get("undefined") and v.get("cr2") is not None
                and np.isfinite(v["cr2"])]
        seed_rows[tag] = {
            "mean": float(np.mean(vals)) if vals else None,
            "median": float(np.median(vals)) if vals else None,
        }
    rows[s] = seed_rows
arms = sorted({t for sr in rows.values() for t in sr})
f9_mean = {s: rows[s].get(f"F9-T1-s{s}", {}).get("mean") for s in rows}
f9_med = {s: rows[s].get(f"F9-T1-s{s}", {}).get("median") for s in rows}
table = []
for a in arms:
    ent = {"arm": a.replace("-s1", "").replace("-s2", "").replace("-s3", ""),
           "mean_by_seed": [rows[s].get(a, {}).get("mean") for s in sorted(rows)],
           "median_by_seed": [rows[s].get(a, {}).get("median") for s in sorted(rows)]}
    ent["gap_to_F9_mean"] = [None if rows[s].get(a, {}).get("mean") is None or f9_mean[s] is None
                             else rows[s][a]["mean"] - f9_mean[s] for s in sorted(rows)]
    ent["gap_to_F9_median"] = [None if rows[s].get(a, {}).get("median") is None or f9_med[s] is None
                               else rows[s][a]["median"] - f9_med[s] for s in sorted(rows)]
    table.append(ent)
(HERE / "REACHABILITY.json").write_text(json.dumps(
    {"f9_mean_by_seed": f9_mean, "f9_median_by_seed": f9_med, "table": table},
    indent=1, default=str), encoding="utf-8")
for ent in table:
    gm = [f"{x:+.3f}" for x in ent["gap_to_F9_mean"] if x is not None]
    gd = [f"{x:+.3f}" for x in ent["gap_to_F9_median"] if x is not None]
    print(f"{ent['arm']:14s} gapF9_mean={gm} gapF9_median={gd}")
