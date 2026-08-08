"""XP1-E: destructive control for left-censoring.

BLK-METZ-60 keeps only cells above the pKi = 4.0 floor.  Conditioning on
`y > floor` is conditioning on the outcome, which can by itself induce apparent
non-additivity in a truly additive panel.  This control regenerates a panel that
is additive BY CONSTRUCTION with the real marginal effects and the real noise
scale, applies the identical floor and the identical XP1-B pipeline, and asks
whether the interaction arms still appear to work.

Falsification rule: if the synthetic additive panel reproduces a materially
positive A4 or AO1 `R2_gamma`, the XP1-B result is an artefact of truncation and
must be withdrawn.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xp1b_transfer as XB  # noqa: E402
from panels import additive_fit, load_metz  # noqa: E402

FLOOR = 4.0


def synthetic_additive(seed=0, density=0.60, interaction_sd=0.0, rank=8):
    """Additive (or additive + planted rank-r interaction) panel with the real
    marginals, the real noise scale and the same censoring floor."""
    Y, M, cid, kin = load_metz(density)
    mu, a, b, fit = additive_fit(Y, M)
    sd = float((Y - fit)[M].std(ddof=1))
    rng = np.random.default_rng(seed)
    n, p = Y.shape
    base = mu + a[:, None] + b[None, :]
    G = np.zeros_like(base)
    if interaction_sd > 0:
        U = rng.normal(size=(n, rank))
        V = rng.normal(size=(p, rank))
        G = U @ V.T
        G *= interaction_sd / G.std()
        sd = max(float(np.sqrt(max(sd ** 2 - interaction_sd ** 2, 1e-6))), 1e-3)
    Ys = base + G + rng.normal(0, sd, size=base.shape)
    # The additive fit is estimated on the *uncensored* cells of an already
    # peeled block, so its fitted values sit above the raw floor.  What must be
    # reproduced is the truncation FRACTION, not the nominal floor value, so the
    # synthetic floor is set to the quantile that matches the real block.
    keep_frac = float(M.mean())
    floor = float(np.quantile(Ys, 1.0 - keep_frac))
    Ms = Ys > floor
    return Ys, Ms, kin, dict(noise_sd=sd, planted_sd=interaction_sd,
                             real_resid_sd=float((Y - fit)[M].std(ddof=1)),
                             synth_floor=floor, nominal_floor=FLOOR,
                             synth_density=float(Ms.mean()),
                             real_density=keep_frac)


def run_on(Ys, Ms, kin, closure="group", rank=8, k_support=16, seeds=(0, 1, 2)):
    """Reuse the XP1-B evaluator by monkey-patching its panel loader."""
    orig = XB.load_metz
    XB.load_metz = lambda density=0.60: (Ys, Ms, np.arange(Ys.shape[0]), kin)
    try:
        return XB.run(closure=closure, rank=rank, k_support=k_support,
                      seeds=seeds, verbose=False, n_boot=1000)
    finally:
        XB.load_metz = orig


if __name__ == "__main__":
    # sd 0.456 = the reproducible interaction sd implied by XP1-A's rank curve
    # (additive CV RMSE 0.7174, rank-8 CV RMSE 0.5541).
    out = {}
    scenarios = [("additive_only", 0.0, 8),
                 ("planted_rank1_sd0.456", 0.456, 1),
                 ("planted_rank3_sd0.456", 0.456, 3),
                 ("planted_rank8_sd0.456", 0.456, 8)]
    for tag, isd, prank in scenarios:
        Ys, Ms, kin, meta = synthetic_additive(seed=11, interaction_sd=isd, rank=prank)
        meta["planted_rank"] = prank
        print(f"\n### {tag}: {meta}")
        res = run_on(Ys, Ms, kin)
        res["synthetic_meta"] = meta
        out[tag] = res
        for a in ("A1", "A4", "A6", "AO1", "A3::pocket_identity_kernel"):
            v = res["arms"].get(a, {})
            r = v.get("r2_gamma_vs_A2")
            print(f"   {a:30s} RMSE={v.get('rmse', float('nan')):.4f}  "
                  f"R2_gamma={r['point']:+.4f} [{r['ci95'][0]:+.4f},{r['ci95'][1]:+.4f}]"
                  if r else f"   {a:30s} RMSE={v.get('rmse', float('nan')):.4f}")
        for k in ("Delta_interaction__A2_minus_A4", "Delta_specific__A6_minus_A4",
                  "Delta_oracle__A2_minus_AO1"):
            c = res["contrasts"][k]
            print(f"   {k:34s} {c['point']:+.5f} [{c['ci95'][0]:+.5f},{c['ci95'][1]:+.5f}]")
    p = os.path.join(XB.REPORT, "xp1e_truncation_control.json")
    json.dump(out, open(p, "w"), indent=2, default=float)
    print("\nwrote", p)
