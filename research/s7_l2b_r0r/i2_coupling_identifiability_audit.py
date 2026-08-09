"""S7_L2B I-2 — coupling identifiability audit (label-side, no training, no GPU).

The question this answers, before any model is trained:

    Given the true residue marginal d_i and the true atom marginal e_j, is the
    SPECIFIC pairing pattern of the contact matrix predictable at all, or is it
    indistinguishable from a uniformly random matrix with the same margins?

If it is indistinguishable, then no pair-coupling model can identify coupling
from these labels, and the failure is LABEL/DATA insufficiency rather than
biological representation failure. That distinction is required by the task.

Objects, following the registered coupling definition:

    G_ij = mu + alpha_i + beta_j + C_ij
    C    = G - Proj_W(additive marginal space)   on the ACTUAL mask, with the
                                                 preregistered pair weights W

The projection is solved as a genuine weighted least-squares ANOVA on the mask.
Naive double centering is NOT used; it is only used as a SELF-TEST oracle in the
special case (complete mask, uniform W) where the two provably coincide.

Null model: degree-preserving bipartite rewiring by checkerboard swaps, which
holds every d_i and every e_j exactly fixed while randomising the pairing.
Rewiring is an EVALUATION control here and is never used as a training negative.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import average_precision, component_macro, paired_bootstrap  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
SEED = 20260815
N_REWIRE = 20            # independent rewirings per complex
SWAP_MULTIPLIER = 30     # swap attempts = SWAP_MULTIPLIER * n_positives


# ------------------------------------------------------------------ ANOVA
def weighted_additive_fit(Y, M, W, iters=200, tol=1e-12):
    """Proj_W(additive space) on mask M with weights W, by weighted ALS.

    Solves  min_{mu,alpha,beta}  sum_{M} W (Y - mu - alpha_i - beta_j)^2
    with the identifiability constraint that the W-weighted means of alpha and
    beta vanish. This is the general solver; it does NOT assume a complete mask
    or uniform weights.
    """
    Wm = W * M
    tot = Wm.sum()
    if tot <= 0:
        return np.zeros_like(Y), 0.0, np.zeros(Y.shape[0]), np.zeros(Y.shape[1])
    mu = (Wm * Y).sum() / tot
    a = np.zeros(Y.shape[0])
    b = np.zeros(Y.shape[1])
    rw = Wm.sum(1)
    cw = Wm.sum(0)
    for _ in range(iters):
        prev = (mu, a.copy(), b.copy())
        r = Y - mu - b[None, :]
        a = np.where(rw > 0, (Wm * r).sum(1) / np.maximum(rw, 1e-12), 0.0)
        r = Y - mu - a[:, None]
        b = np.where(cw > 0, (Wm * r).sum(0) / np.maximum(cw, 1e-12), 0.0)
        # re-centre so the decomposition is identified
        am = (rw * a).sum() / tot
        bm = (cw * b).sum() / tot
        a -= am
        b -= bm
        mu += am + bm
        if (abs(mu - prev[0]) < tol and np.abs(a - prev[1]).max() < tol
                and np.abs(b - prev[2]).max() < tol):
            break
    return mu + a[:, None] + b[None, :], mu, a, b


def selftest_anova():
    """In the special case of a complete mask and uniform weights, the weighted
    ALS solution must equal classical double centering. This validates the
    general solver without ever using double centering on real data."""
    rng = np.random.default_rng(0)
    Y = rng.normal(size=(23, 17))
    M = np.ones_like(Y)
    W = np.ones_like(Y)
    fit, _mu, _a, _b = weighted_additive_fit(Y, M, W)
    dc = (Y.mean(1)[:, None] + Y.mean(0)[None, :] - Y.mean())
    err = float(np.abs(fit - dc).max())
    # and a genuinely non-uniform case must NOT equal double centering
    W2 = rng.random(Y.shape) + 0.1
    fit2, *_ = weighted_additive_fit(Y, M, W2)
    diff2 = float(np.abs(fit2 - dc).max())
    return {"uniform_matches_double_centering_max_abs_err": err,
            "uniform_selftest_pass": bool(err < 1e-8),
            "weighted_differs_from_double_centering_max_abs": diff2,
            "weighted_is_genuinely_different": bool(diff2 > 1e-6)}


# ------------------------------------------------------- degree-preserving null
def rewire(Y, rng, swaps):
    """Checkerboard swaps: preserve every row and column degree exactly."""
    A = Y.copy()
    rows, cols = A.shape
    if rows < 2 or cols < 2:
        return A
    ones = np.argwhere(A == 1)
    if len(ones) < 2:
        return A
    done = 0
    for _ in range(swaps * 4):
        if done >= swaps:
            break
        i1, j1 = ones[rng.integers(len(ones))]
        i2, j2 = ones[rng.integers(len(ones))]
        if i1 == i2 or j1 == j2:
            continue
        if A[i1, j2] == 0 and A[i2, j1] == 0:
            A[i1, j1] = A[i2, j2] = 0
            A[i1, j2] = A[i2, j1] = 1
            ones = np.argwhere(A == 1)
            done += 1
    return A


def active_submatrix(Y):
    r = np.where(Y.sum(1) > 0)[0]
    c = np.where(Y.sum(0) > 0)[0]
    return Y[np.ix_(r, c)], len(r), len(c)


def coupling_stat(Y):
    """Leading-singular-value share of the marginal-orthogonal residual."""
    M = np.ones_like(Y, dtype=float)
    W = np.ones_like(Y, dtype=float)
    fit, *_ = weighted_additive_fit(Y.astype(float), M, W)
    C = Y - fit
    fro = float(np.linalg.norm(C))
    if fro < 1e-12 or min(C.shape) < 2:
        return None, fro
    s = np.linalg.svd(C, compute_uv=False)
    return float(s[0] / fro), fro


def main():
    st = selftest_anova()
    print("ANOVA self-test:", json.dumps(st), flush=True)
    if not st["uniform_selftest_pass"] or not st["weighted_is_genuinely_different"]:
        (OUT / "I2_COUPLING_IDENTIFIABILITY_AUDIT.json").write_text(
            json.dumps({"verdict": "ANOVA_SOLVER_CONTRACT_FAIL", "selftest": st},
                       indent=2), encoding="utf-8")
        return 1

    kept, _q, contract, _f = build()
    comp_of = protein_components(kept)
    _tr, _ha, held_A, _hb = make_split(kept, comp_of)
    print(f"held-out A complexes: {len(held_A)}", flush=True)

    rng = np.random.default_rng(SEED)
    ap_R, ap_A, ap_RA, ap_exact = {}, {}, {}, {}
    true_stat, null_stat, null_sd = [], [], []
    excess = []
    skipped = Counter()

    for n, rec in enumerate(held_A):
        L, Anum = rec["n_res"], rec["n_atoms"]
        Y = np.zeros((L, Anum), dtype=np.int8)
        for i, j in rec["edges"]:
            Y[i, j] = 1
        if Y.sum() == 0:
            skipped["no_positives"] += 1
            continue
        d = Y.sum(1).astype(float)
        e = Y.sum(0).astype(float)
        ri = np.repeat(np.arange(L), Anum)
        ai = np.tile(np.arange(Anum), L)
        yflat = Y.ravel()
        key = rec["source_key"]
        ap_R[key] = average_precision(np.repeat(d, Anum), yflat, ri, ai)
        ap_A[key] = average_precision(np.tile(e, L), yflat, ri, ai)
        fit, *_ = weighted_additive_fit(Y.astype(float), np.ones_like(Y, float),
                                        np.ones_like(Y, float))
        ap_RA[key] = average_precision(fit.ravel(), yflat, ri, ai)
        ap_exact[key] = average_precision(yflat.astype(float), yflat, ri, ai)

        # coupling beyond margins, on the ACTIVE submatrix
        Ya, nr, nc = active_submatrix(Y)
        if nr < 3 or nc < 3 or Ya.sum() < 4:
            skipped["active_submatrix_too_small"] += 1
            continue
        t_stat, _fro = coupling_stat(Ya.astype(float))
        if t_stat is None:
            skipped["degenerate_residual"] += 1
            continue
        nulls = []
        for _ in range(N_REWIRE):
            Yr = rewire(Ya, rng, SWAP_MULTIPLIER * int(Ya.sum()))
            s, _ = coupling_stat(Yr.astype(float))
            if s is not None:
                nulls.append(s)
        if len(nulls) < 5:
            skipped["insufficient_nulls"] += 1
            continue
        true_stat.append(t_stat)
        null_stat.append(float(np.mean(nulls)))
        null_sd.append(float(np.std(nulls)))
        excess.append((t_stat - np.mean(nulls)) / max(np.std(nulls), 1e-9))
        if (n + 1) % 400 == 0:
            print(f"  {n+1}/{len(held_A)}", flush=True)

    cm_R, macro_R = component_macro(ap_R, comp_of)
    cm_A, macro_A = component_macro(ap_A, comp_of)
    cm_RA, macro_RA = component_macro(ap_RA, comp_of)
    _cm_E, macro_E = component_macro(ap_exact, comp_of)

    excess = np.array(excess)
    true_stat = np.array(true_stat)
    null_stat = np.array(null_stat)
    frac_above = float((excess > 0).mean())
    med_z = float(np.median(excess))

    out = {
        "schema": "MetaSieve.S7L2B.I2.CouplingIdentifiabilityAudit.v1",
        "created_utc": "2026-08-09",
        "anova_selftest": st,
        "held_out_A_complexes": len(held_A),
        "complexes_in_coupling_test": int(len(true_stat)),
        "skipped": dict(skipped),
        "oracle_marginal_ap_macro_over_components": {
            "Oracle_R_true_residue_marginal": macro_R,
            "Oracle_A_true_atom_marginal": macro_A,
            "Oracle_RA_additive_projection": macro_RA,
            "exact_pair_upper_bound": macro_E,
            "measured_B4_for_reference": 0.022949,
            "measured_ligand_only_for_reference": 0.004499,
        },
        "coupling_beyond_margins": {
            "statistic": "leading singular value share of the marginal-orthogonal "
                         "residual on the active submatrix",
            "null": f"degree-preserving checkerboard rewiring, {N_REWIRE} per complex, "
                    f"{SWAP_MULTIPLIER}x positives swap attempts",
            "true_mean": float(true_stat.mean()),
            "null_mean": float(null_stat.mean()),
            "difference": float(true_stat.mean() - null_stat.mean()),
            "median_z_vs_null": med_z,
            "fraction_of_complexes_above_own_null": frac_above,
        },
    }
    # headroom reading
    out["headroom"] = {
        "oracle_RA_minus_B4": macro_RA - 0.022949,
        "oracle_R_minus_B4": macro_R - 0.022949,
        "reading": "how much AP a model could still gain if it recovered the TRUE "
                   "marginals; this bounds the marginal-recovery problem separately "
                   "from the coupling problem",
    }
    (OUT / "I2_COUPLING_IDENTIFIABILITY_AUDIT.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in
                      ("oracle_marginal_ap_macro_over_components",
                       "coupling_beyond_margins", "headroom")}, indent=2))
    print(f"\nwrote {OUT / 'I2_COUPLING_IDENTIFIABILITY_AUDIT.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
