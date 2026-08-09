"""S3 — power design for the structural test, computed from the realised split.

Simulation-based: inject a known signal of size R2_true into a synthetic channel
with the SAME cluster structure as the realised split, run the SAME cluster
bootstrap, and record the detection rate at the registered decision rule
(95% CI lower bound > 0).  This yields the minimum detectable effect, which is
then frozen as the S6 effect floor BEFORE the structural test set is opened.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CACHE = r"D:\MetaSieve\dataset\processed\ssl_b2"
OUT = r"D:\MetaSieve\report\ssl_b2_structural_observability"


def detect_rate(n_units, per_unit_n, r2_true, n_rep=400, n_boot=600, seed=0):
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_rep):
        sig = np.sqrt(max(r2_true, 0.0))
        sse_arm, sse_base, nn = [], [], []
        for u in range(n_units):
            k = per_unit_n[u % len(per_unit_n)]
            y = rng.normal(size=k)
            # arm explains a fraction r2_true of the variance, with a cluster effect
            pred = sig * (y + rng.normal(scale=0.35))
            sse_arm.append(float(((y - pred) ** 2).sum()))
            sse_base.append(float((y ** 2).sum()))
            nn.append(k)
        sa, sb = np.array(sse_arm), np.array(sse_base)
        idx = rng.integers(0, n_units, (n_boot, n_units))
        bs = 1.0 - sa[idx].sum(1) / np.maximum(sb[idx].sum(1), 1e-12)
        if np.percentile(bs, 2.5) > 0:
            hits += 1
    return hits / n_rep


def main():
    p = os.path.join(CACHE, "teacher_dataset.npz")
    if not os.path.exists(p):
        print("dataset not built yet")
        return
    d = np.load(p, allow_pickle=True)
    pclus, scaf = d["pclus"], d["scaffold"]
    n = len(d["Y"])
    units_p = len(set(pclus))
    units_s = len(set(scaf))
    sizes = np.bincount(np.unique(pclus, return_inverse=True)[1])
    # only the double-held-out diagonal is scored, ~1/5 of complexes per fold
    eff_units = max(int(round(units_p / 5)), 2)
    eff_sizes = [max(int(x), 1) for x in sizes[:50]] or [1]

    curve = {}
    for r2 in (0.02, 0.05, 0.10, 0.15, 0.20, 0.30):
        curve[str(r2)] = detect_rate(eff_units, eff_sizes, r2)
    mde = None
    for r2 in sorted(float(k) for k in curve):
        if curve[str(r2)] >= 0.80:
            mde = r2
            break

    out = {"schema": "MetaSieve.StructuralPowerAnalysis.v1",
           "complexes": int(n),
           "protein_clusters_total": int(units_p),
           "scaffolds_total": int(units_s),
           "effective_units_per_fold": eff_units,
           "decision_rule": "95% cluster-bootstrap CI lower bound > 0",
           "detection_rate_by_true_r2": curve,
           "minimum_detectable_r2_at_80pc_power": mde,
           "frozen_S6_effect_floor": mde if mde is not None else 0.30,
           "note": "computed from the realised split BEFORE the structural test "
                   "set was scored; the S6 floor may not be lowered afterwards"}
    json.dump(out, open(os.path.join(OUT, "STRUCTURAL_POWER_ANALYSIS.json"), "w"),
              indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
