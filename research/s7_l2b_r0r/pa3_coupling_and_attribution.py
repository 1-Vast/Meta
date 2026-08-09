"""Phase 2A / Phases 3+4 — weighted marginal/coupling decomposition and the
matched attribution battery, including the degree-preserving rewiring null.

Registered by PREREG_S7_L2B_PHASE2A.md (sha 4e01401d...), sections 7 and 8,
plus computational amendment 01 (exact tied-AP expectation) and amendment 02
(closed-form harmonic evaluation of that expectation, and ridge-IRLS budget).

Three blocks are kept strictly separate and separately labelled:
  DEPLOYABLE   — decomposition of the sealed prediction logits (label-free)
  ORACLE       — label-fitted ceilings (Rasch additive null); never a model arm
  NULL         — evaluation-only nulls (degree-preserving rewiring, shuffles)

Nothing is trained. No affinity source is opened.
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import component_macro, paired_bootstrap  # noqa: E402
from i2_coupling_identifiability_audit import weighted_additive_fit, selftest_anova  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
S4 = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "sealed_preds"
S5 = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "sealed_preds_b5"
PREREG_SHA = "4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e"

SEED_REWIRE = 20260817
SEED_BOOT = 20260818
SEED_SHUFFLE = 20260821
N_REWIRE = 20
BURN_IN_MULT = 100
BETWEEN_MULT = 30
ORTHO_TOL = 1e-8
COUPLING_MARGIN = 0.01
TC_MEDIAN_Z = 2.0
TC_FRAC_ABOVE = 0.60
MIX_CHECKPOINTS = (0, 1, 5, 10, 30)
RASCH_NEWTON = 6
RASCH_ALS = 12
RASCH_RIDGE = 1e-6


# ---------------------------------------------------------------- tied AP
# Exact E[AP] under a uniformly random ordering inside every tied block.
# Per block of size n with k positives, preceded by a items containing b
# positives, the expected sum of precisions over that block's positives is
#     (k/n) * [ (b+1) * S1 + (k-1)/(n-1) * (n - (a+1) * S1) ],
#     S1 = sum_{j=1..n} 1/(a+j) = H(a+n) - H(a).
# Using a cumulative harmonic table this is fully vectorised over blocks, with
# no Python loop, which is what makes the rewiring null tractable.
def ap_prep(scores: np.ndarray):
    order = np.argsort(-scores, kind="stable")
    s = scores[order]
    bounds = np.flatnonzero(np.r_[True, s[1:] != s[:-1], True])
    N = scores.size
    harm = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, N + 1, dtype=np.float64))))
    return order, bounds, harm


def ap_apply(prep, labels: np.ndarray):
    order, bounds, harm = prep
    y = labels[order].astype(np.float64)
    P = y.sum()
    if P == 0 or P == y.size:
        return None
    csum = np.concatenate(([0.0], np.cumsum(y)))
    lo = bounds[:-1].astype(np.int64)
    hi = bounds[1:].astype(np.int64)
    n = (hi - lo).astype(np.float64)
    k = csum[hi] - csum[lo]
    a = lo.astype(np.float64)
    b = csum[lo]
    m = k > 0
    if not m.any():
        return 0.0
    n, k, a, b = n[m], k[m], a[m], b[m]
    s1 = harm[hi[m]] - harm[lo[m]]
    grow = np.where(n > 1, (k - 1.0) / np.maximum(n - 1.0, 1.0), 0.0)
    contrib = (k / n) * ((b + 1.0) * s1 + grow * (n - (a + 1.0) * s1))
    return float(contrib.sum() / P)


def selftest_ap():
    """The closed form must reproduce plain AP when there are no ties, and the
    Monte-Carlo expectation when there are."""
    rng = np.random.default_rng(3)
    worst_plain = worst_mc = 0.0
    for _ in range(8):
        N = int(rng.integers(200, 800))
        s = rng.normal(size=N)
        y = (rng.random(N) < 0.05).astype(np.int8)
        if y.sum() == 0:
            continue
        o = np.argsort(-s, kind="stable")
        yy = y[o].astype(float)
        plain = float((np.cumsum(yy) / np.arange(1, N + 1) * yy).sum() / yy.sum())
        worst_plain = max(worst_plain, abs(ap_apply(ap_prep(s), y) - plain))
        sb = (rng.random(N) < 0.2).astype(np.float64)
        vals = []
        for _ in range(400):
            oo = np.lexsort((rng.random(N), -sb))
            z = y[oo].astype(float)
            vals.append(float((np.cumsum(z) / np.arange(1, N + 1) * z).sum() / z.sum()))
        worst_mc = max(worst_mc, abs(ap_apply(ap_prep(sb), y) - float(np.mean(vals))))
    return {"max_abs_diff_untied_vs_plain_ap": worst_plain,
            "max_abs_diff_tied_vs_monte_carlo": worst_mc,
            "pass": bool(worst_plain < 1e-9 and worst_mc < 0.03)}


# ------------------------------------------- efficient degree-preserving swap
def rewire_stream(edges, occupied, rng, n_swaps, max_attempts_mult=8):
    """Checkerboard swaps preserving every row and column degree exactly."""
    m = len(edges)
    if m < 2 or n_swaps <= 0:
        return 0, 0
    done = attempts = 0
    cap = max_attempts_mult * n_swaps
    while done < n_swaps and attempts < cap:
        attempts += 1
        i1 = int(rng.integers(m))
        i2 = int(rng.integers(m))
        if i1 == i2:
            continue
        r1, a1 = edges[i1]
        r2, a2 = edges[i2]
        if r1 == r2 or a1 == a2:
            continue
        if (r1, a2) in occupied or (r2, a1) in occupied:
            continue
        occupied.discard((r1, a1))
        occupied.discard((r2, a2))
        occupied.add((r1, a2))
        occupied.add((r2, a1))
        edges[i1] = (r1, a2)
        edges[i2] = (r2, a1)
        done += 1
    return done, attempts


def coupling_stat(Y):
    """Leading singular-value share of the marginal-orthogonal residual on the
    active submatrix. Identical statistic to I-2, kept unchanged on purpose."""
    M = np.ones_like(Y, dtype=float)
    fit, *_ = weighted_additive_fit(Y.astype(float), M, M)
    C = Y - fit
    fro = float(np.linalg.norm(C))
    if fro < 1e-12 or min(C.shape) < 2:
        return None
    return float(np.linalg.svd(C, compute_uv=False)[0] / fro)


def ridge_additive_fit(Y, W, lam, iters=30, tol=1e-10):
    """Weighted additive fit with a RIDGE ON THE COEFFICIENTS (not on the
    weights). Needed because the unpenalised Rasch MLE does not exist here:
    the matrix is ~0.07% positive, so almost every residue row and atom column
    contains no positive at all and its coefficient is completely separated."""
    tot = W.sum()
    mu = float((W * Y).sum() / tot)
    a = np.zeros(Y.shape[0])
    b = np.zeros(Y.shape[1])
    rw, cw = W.sum(1), W.sum(0)
    for _ in range(iters):
        prev_a, prev_b = a.copy(), b.copy()
        a = (W * (Y - mu - b[None, :])).sum(1) / (rw + lam)
        b = (W * (Y - mu - a[:, None])).sum(0) / (cw + lam)
        if max(np.abs(a - prev_a).max(), np.abs(b - prev_b).max()) < tol:
            break
    return mu + a[:, None] + b[None, :]


def rasch_fit(Yf, newton=RASCH_NEWTON, als=RASCH_ALS, ridge=RASCH_RIDGE):
    """Label-fitted additive logistic null  logit P(Y=1) = mu + alpha_r + beta_a,
    by damped IRLS with a coefficient ridge. ORACLE ONLY; never a model arm.

    The unpenalised MLE is non-existent (complete separation), so the reported
    object is explicitly the ridge-penalised fit and its achieved step size is
    reported rather than assumed converged."""
    L, A = Yf.shape
    p = float(np.clip(Yf.mean(), 1e-6, 1 - 1e-6))
    eta = np.full((L, A), np.log(p / (1 - p)))
    delta = float("nan")
    for _ in range(newton):
        mu_ = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        w = np.clip(mu_ * (1 - mu_), 1e-6, None)
        zz = eta + (Yf - mu_) / w
        # relative ridge: scaled to the mean per-coefficient weight mass
        lam = ridge_scale(w, L, A)
        fit = np.clip(ridge_additive_fit(zz, w, lam, iters=als), -30.0, 30.0)
        step = 0.5 * (fit - eta)                       # damping
        delta = float(np.abs(step).max())
        eta = eta + step
        if delta < 1e-6:
            break
    return eta, delta


def ridge_scale(w, L, A):
    """Ridge large enough to keep separated coefficients finite: 1e-2 of the
    typical per-row/per-column weight mass."""
    return 1e-2 * float(w.sum()) / max(L + A, 1)


def main():
    t0 = time.time()
    st_anova = selftest_anova()
    st_ap = selftest_ap()
    print("ANOVA self-test:", json.dumps(st_anova), flush=True)
    print("AP self-test:", json.dumps(st_ap), flush=True)
    if not (st_anova["uniform_selftest_pass"] and st_anova["weighted_is_genuinely_different"]
            and st_ap["pass"]):
        (OUT / "PHASE2A_MARGINAL_COUPLING_AUDIT.json").write_text(
            json.dumps({"verdict": "ESTIMATOR_CONTRACT_FAIL",
                        "anova": st_anova, "ap": st_ap}, indent=2), encoding="utf-8")
        return 1

    kept, _q, _c, _f = build()
    comp_of = protein_components(kept)
    _tr, _ha, held_A, _hb = make_split(kept, comp_of)
    idx = json.loads((S4 / "heldoutA_index.json").read_text(encoding="utf-8"))
    held_A = [r for r in held_A if r["source_key"] in idx]
    total = max(v[0] + v[1] * v[2] for v in idx.values())

    arms = {}
    for name, d in (("BL", S4), ("B4", S4), ("B5", S5), ("BP5", S5), ("BX5", S5)):
        arms[name] = np.memmap(d / f"heldoutA_{name}.f16.dat", dtype=np.float16,
                               mode="r", shape=(total,))
    DECOMP = ("B5", "B4", "BX5", "BP5", "BL")
    NULLARMS = ("B5", "B4", "BX5")

    per = {a: {w: {} for w in ("full", "res", "atom", "add", "coup")} for a in DECOMP}
    coup_null = {a: {} for a in NULLARMS}
    coup_shuf = {a: {} for a in NULLARMS}
    teacher_true, teacher_null_mean, teacher_z = {}, {}, {}
    oracle_rasch, oracle_exact = {}, {}
    oracle_ls_add, oracle_R, oracle_A = {}, {}, {}
    ortho_max = dc_max = 0.0
    rasch_delta = []
    skipped = Counter()
    mixing = {c: [] for c in MIX_CHECKPOINTS}
    mix_between = []
    degree_violations = 0
    nonswitchable = set()

    rng_rw = np.random.default_rng(SEED_REWIRE)
    rng_mix = np.random.default_rng(SEED_REWIRE + 1)
    rng_sh = np.random.default_rng(SEED_SHUFFLE)

    for n_i, rec in enumerate(held_A):
        key = rec["source_key"]
        off, L, A = idx[key]
        N = L * A
        yflat = np.zeros(N, dtype=np.int8)
        for i, j in rec["edges"]:
            yflat[i * A + j] = 1
        if yflat.sum() == 0:
            skipped["no_positives"] += 1
            continue

        cprep = {}
        for a in DECOMP:
            G = np.asarray(arms[a][off:off + N], dtype=np.float64).reshape(L, A)
            ones = np.ones_like(G)
            fit, mu, alpha, beta = weighted_additive_fit(G, ones, ones)
            C = G - fit
            xtc = np.concatenate(([C.sum()], C.sum(1), C.sum(0)))
            ortho_max = max(ortho_max, float(np.linalg.norm(xtc) / (1.0 + np.linalg.norm(C))))
            if a == "B5":
                dc = G.mean(1)[:, None] + G.mean(0)[None, :] - G.mean()
                dc_max = max(dc_max, float(np.abs(fit - dc).max()))
            per[a]["full"][key] = ap_apply(ap_prep(G.ravel()), yflat)
            per[a]["res"][key] = ap_apply(ap_prep(np.repeat(alpha, A)), yflat)
            per[a]["atom"][key] = ap_apply(ap_prep(np.tile(beta, L)), yflat)
            per[a]["add"][key] = ap_apply(ap_prep(fit.ravel()), yflat)
            pc = ap_prep(C.ravel())
            per[a]["coup"][key] = ap_apply(pc, yflat)
            if a in NULLARMS:
                cprep[a] = pc
                coup_shuf[a][key] = ap_apply(ap_prep(rng_sh.permutation(C.ravel())), yflat)

        # ------------------------------- NULL: degree-preserving label rewiring
        act_r = np.flatnonzero(yflat.reshape(L, A).sum(1) > 0)
        act_a = np.flatnonzero(yflat.reshape(L, A).sum(0) > 0)
        E = int(yflat.sum())
        Ysub = yflat.reshape(L, A)[np.ix_(act_r, act_a)].astype(float)
        switchable = (E >= 4 and act_r.size >= 3 and act_a.size >= 3)
        if switchable:
            edges = [(int(i), int(j)) for i, j in rec["edges"]]
            occupied = set(edges)
            d_row0 = np.bincount([e[0] for e in edges], minlength=L)
            d_col0 = np.bincount([e[1] for e in edges], minlength=A)
            orig = set(edges)
            me, mo, prev = list(edges), set(edges), 0
            for cp in MIX_CHECKPOINTS:
                if cp * E > prev:
                    rewire_stream(me, mo, rng_mix, cp * E - prev)
                    prev = cp * E
                mixing[cp].append(len(orig & set(me)) / E)
            done, _att = rewire_stream(edges, occupied, rng_rw, BURN_IN_MULT * E)
            if done < 0.5 * BURN_IN_MULT * E:
                switchable = False
        if not switchable:
            nonswitchable.add(key)
        else:
            tstat = coupling_stat(Ysub)
            null_ap = {a: [] for a in NULLARMS}
            tnull, last = [], None
            for _s in range(N_REWIRE):
                rewire_stream(edges, occupied, rng_rw, BETWEEN_MULT * E)
                if not (np.array_equal(np.bincount([e[0] for e in edges], minlength=L), d_row0)
                        and np.array_equal(np.bincount([e[1] for e in edges], minlength=A), d_col0)):
                    degree_violations += 1
                yr = np.zeros(N, dtype=np.int8)
                for i, j in edges:
                    yr[i * A + j] = 1
                for a in NULLARMS:
                    v = ap_apply(cprep[a], yr)
                    if v is not None:
                        null_ap[a].append(v)
                ts = coupling_stat(yr.reshape(L, A)[np.ix_(act_r, act_a)].astype(float))
                if ts is not None:
                    tnull.append(ts)
                cur = set(edges)
                if last is not None:
                    mix_between.append(len(last & cur) / E)
                last = cur
            for a in NULLARMS:
                if null_ap[a]:
                    coup_null[a][key] = float(np.mean(null_ap[a]))
            if tstat is not None and len(tnull) >= 5:
                teacher_true[key] = tstat
                teacher_null_mean[key] = float(np.mean(tnull))
                teacher_z[key] = float((tstat - np.mean(tnull)) / max(np.std(tnull), 1e-9))

        # ------------------------------------- ORACLE (label-fitted ceilings)
        Yf = yflat.reshape(L, A).astype(np.float64)
        eta, dlt = rasch_fit(Yf)
        rasch_delta.append(dlt)
        oracle_rasch[key] = ap_apply(ap_prep(eta.ravel()), yflat)
        oracle_exact[key] = ap_apply(ap_prep(yflat.astype(np.float64)), yflat)
        ls_fit, _m, la, lb = weighted_additive_fit(Yf, np.ones_like(Yf), np.ones_like(Yf))
        oracle_ls_add[key] = ap_apply(ap_prep(ls_fit.ravel()), yflat)
        oracle_R[key] = ap_apply(ap_prep(np.repeat(Yf.sum(1), A)), yflat)
        oracle_A[key] = ap_apply(ap_prep(np.tile(Yf.sum(0), L)), yflat)

        if (n_i + 1) % 200 == 0:
            print(f"  {n_i+1}/{len(held_A)}  {time.time()-t0:.0f}s", flush=True)

    # ------------------------------------------------------------- summarise
    def mac(d):
        return component_macro({k: v for k, v in d.items() if v is not None}, comp_of)

    block_deploy, comp_tables = {}, {}
    for a in DECOMP:
        block_deploy[a], comp_tables[a] = {}, {}
        for w in ("full", "res", "atom", "add", "coup"):
            cm, m = mac(per[a][w])
            block_deploy[a][w] = m
            comp_tables[a][w] = cm

    _c1, m_rasch = mac(oracle_rasch)
    _c2, m_ex = mac(oracle_exact)
    _c3, m_ls = mac(oracle_ls_add)
    _c4, m_or = mac(oracle_R)
    _c5, m_oa = mac(oracle_A)
    null_tables = {a: mac(coup_null[a]) for a in NULLARMS}
    shuf_tables = {a: mac(coup_shuf[a]) for a in NULLARMS}

    tz = np.array(list(teacher_z.values())) if teacher_z else np.array([])
    tc_median_z = float(np.median(tz)) if tz.size else float("nan")
    tc_frac = float((tz > 0).mean()) if tz.size else float("nan")
    TC = bool(tz.size and tc_median_z >= TC_MEDIAN_Z and tc_frac >= TC_FRAC_ABOVE)

    bc_null = paired_bootstrap(comp_tables["B5"]["coup"], null_tables["B5"][0],
                               n_boot=10000, seed=SEED_BOOT)
    bc_null["pass"] = bool(bc_null["delta"] >= COUPLING_MARGIN and bc_null["lcb95_one_sided"] > 0)
    bc_bx = paired_bootstrap(comp_tables["B5"]["coup"], comp_tables["BX5"]["coup"],
                             n_boot=10000, seed=SEED_BOOT)
    bc_bx["pass"] = bool(bc_bx["delta"] >= COUPLING_MARGIN and bc_bx["lcb95_one_sided"] > 0)
    BC = bool(bc_null["pass"] and bc_bx["pass"])

    # registered requirement: recompute the point estimate on the switchable
    # subset so the non-switchable exclusion cannot move a verdict silently
    sw_comp = set(null_tables["B5"][0])
    b5_coup_switchable = float(np.mean([v for k, v in comp_tables["B5"]["coup"].items()
                                        if k in sw_comp])) if sw_comp else float("nan")
    bx_coup_switchable = float(np.mean([v for k, v in comp_tables["BX5"]["coup"].items()
                                        if k in sw_comp])) if sw_comp else float("nan")

    res = {
        "schema": "MetaSieve.S7L2B.P2A.MarginalCouplingAndAttribution.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA,
        "amendments": ["PREREG_S7_L2B_PHASE2A_AMENDMENT_01.md",
                       "PREREG_S7_L2B_PHASE2A_AMENDMENT_02.md"],
        "repo_commit": "623602e76b7d4f445af069014782278163183d59",
        "estimator_selftests": {"anova": st_anova, "tie_aware_ap": st_ap},
        "numerics": {
            "max_orthogonality_ratio": ortho_max,
            "orthogonality_tolerance": ORTHO_TOL,
            "orthogonality_pass": bool(ortho_max <= ORTHO_TOL),
            "max_abs_diff_weighted_ALS_vs_double_centering_B5": dc_max,
            "double_centering_admissible": "Phase 0 C5 proved the mask complete and "
                                           "uniformly weighted; weighted ALS is used "
                                           "regardless and double centering only as a "
                                           "numerical self-test",
            "rasch_irls_max_final_delta": float(np.max(rasch_delta)) if rasch_delta else None,
            "rasch_irls_median_final_delta": float(np.median(rasch_delta)) if rasch_delta else None,
        },
        "DEPLOYABLE_prediction_decomposition": {
            "metric": "exact tie-aware E[AP], component-macro over protein closure "
                      "components, complete residue x atom matrix",
            "components": len(comp_tables["B5"]["full"]),
            "arms": block_deploy,
        },
        "ORACLE_label_fitted_ceilings": {
            "warning": "label-fitted; NEVER a deployable model arm and never a Gate",
            "separation_note": "the unpenalised Rasch MLE does not exist on this "
                               "corpus: the matrix is ~0.07% positive, so almost every "
                               "residue row and atom column contains no positive and "
                               "its coefficient is completely separated. The reported "
                               "logistic object is therefore the ridge-penalised fit, "
                               "and the least-squares additive projection of Y is "
                               "reported alongside it as the well-posed ceiling.",
            "rasch_additive_null_ridge_penalised": m_rasch,
            "rasch_converged": bool(rasch_delta and np.max(rasch_delta) < 1e-3),
            "rasch_usable_as_a_ceiling": bool(rasch_delta and np.max(rasch_delta) < 1e-3
                                              and m_rasch >= m_ls),
            "rasch_interpretation": "REPORTED FOR TRANSPARENCY ONLY. The damped "
                                    "ridge-IRLS did not converge (see "
                                    "rasch_irls_max_final_delta) because the design is "
                                    "completely separated. Its AP is therefore NOT a "
                                    "valid ceiling and must not be quoted as one. The "
                                    "well-posed label-fitted additive ceiling on this "
                                    "corpus is least_squares_additive_projection_of_Y.",
            "least_squares_additive_projection_of_Y": m_ls,
            "true_residue_margin_only": m_or,
            "true_atom_margin_only": m_oa,
            "exact_pair_upper_bound": m_ex,
        },
        "NULL_evaluation_only": {
            "degree_preserving_rewiring_coupling_ap": {a: null_tables[a][1] for a in NULLARMS},
            "within_complex_shuffle_of_coupling_ap": {a: shuf_tables[a][1] for a in NULLARMS},
            "rewiring_spec": {
                "swap": "checkerboard; row and column degrees exactly preserved",
                "burn_in_swaps": f"{BURN_IN_MULT} x positives",
                "swaps_between_samples": f"{BETWEEN_MULT} x positives",
                "independent_rewires": N_REWIRE, "seed": SEED_REWIRE,
                "degree_preservation_violations": degree_violations,
                "complexes_total": len(held_A),
                "complexes_non_switchable": len(nonswitchable),
                "complexes_switchable": len(held_A) - len(nonswitchable),
                "components_with_a_rewiring_null": len(null_tables["B5"][0]),
                "components_total": len(comp_tables["B5"]["coup"]),
                "edge_overlap_vs_original_by_swaps_per_edge": {
                    str(c): float(np.mean(v)) for c, v in mixing.items() if v},
                "mean_overlap_between_successive_samples": (
                    float(np.mean(mix_between)) if mix_between else None)},
            "note": "rewiring is an attribution null only; never a non-binder label "
                    "and never a training negative",
        },
        "TC_teacher_edge_coupling": {
            "statistic": "leading singular-value share of the marginal-orthogonal "
                         "residual of Y on the active submatrix (identical to I-2)",
            "complexes": int(tz.size),
            "true_mean": float(np.mean(list(teacher_true.values()))) if teacher_true else None,
            "null_mean": float(np.mean(list(teacher_null_mean.values()))) if teacher_null_mean else None,
            "median_z": tc_median_z, "fraction_above_own_null": tc_frac,
            "thresholds": {"median_z": TC_MEDIAN_Z, "fraction": TC_FRAC_ABOVE},
            "TC": TC,
        },
        "BC_b5_edge_coupling": {
            "b5_coupling_ap": block_deploy["B5"]["coup"],
            "b5_coupling_ap_switchable_subset": b5_coup_switchable,
            "bx5_coupling_ap_switchable_subset": bx_coup_switchable,
            "vs_degree_preserving_null": bc_null,
            "vs_wrong_ligand_BX5_coupling": bc_bx,
            "margin": COUPLING_MARGIN, "BC": BC,
        },
        "matched_attribution_battery": {
            "1_full_B5": block_deploy["B5"]["full"],
            "2_residue_marginal": block_deploy["B5"]["res"],
            "3_atom_marginal": block_deploy["B5"]["atom"],
            "4_weighted_additive": block_deploy["B5"]["add"],
            "5_marginal_orthogonal_coupling": block_deploy["B5"]["coup"],
            "6_wrong_ligand_full": block_deploy["BX5"]["full"],
            "6b_wrong_ligand_coupling": block_deploy["BX5"]["coup"],
            "7_wrong_protein_full": block_deploy["BP5"]["full"],
            "7b_wrong_protein_coupling": block_deploy["BP5"]["coup"],
            "8_within_complex_shuffle_of_coupling": shuf_tables["B5"][1],
            "9_degree_preserving_rewiring_null": null_tables["B5"][1],
        },
        "skipped": dict(skipped),
        "elapsed_sec": round(time.time() - t0, 1),
    }
    (OUT / "PHASE2A_MARGINAL_COUPLING_AUDIT.json").write_text(
        json.dumps(res, indent=2), encoding="utf-8")
    (OUT / "PHASE2A_COMPONENT_TABLES.json").write_text(
        json.dumps({"deployable": comp_tables,
                    "rewiring_null": {a: null_tables[a][0] for a in NULLARMS},
                    "shuffle": {a: shuf_tables[a][0] for a in NULLARMS},
                    "teacher_z": teacher_z}, indent=1), encoding="utf-8")
    print(json.dumps({k: res[k] for k in
                      ("numerics", "DEPLOYABLE_prediction_decomposition",
                       "ORACLE_label_fitted_ceilings", "NULL_evaluation_only",
                       "TC_teacher_edge_coupling", "BC_b5_edge_coupling")},
                     indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
