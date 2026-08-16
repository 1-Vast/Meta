# K-LBP v2 R3 — synthetic estimator sensitivity under the empirical V_t (preregistration)

**Frozen:** 2026-07-27, **before any R3 statistic was computed.**
**Amendment A1 (2026-07-27, pre-run):** the frozen inference text said "component bootstrap of `γ̂`
(10,000 draws over components)" per replicate. At the frozen grid (200 replicates × 5 regimes × 5
fold refits per draw) that is computationally infeasible. **Before any run**, the per-replicate
interval evaluated for coverage in G-R3-1/G-R3-2 is replaced by the **fold-pooled t-CI** (5 fold
estimates → `mean ± t(0.975,4)·sd/√5`), whose coverage is exactly what the 200 replicates measure;
the 10,000-draw component bootstrap is reserved for R4's single real-data interval, where it is
computed once. All gates, thresholds, regimes, replicate counts, and verdict codes are **unchanged**.
This amendment is recorded here, before any result, and may not be revisited after a result.
**Amendment A2 (2026-07-28, pre-run estimator and regime repair):** no R3 result artifact exists and no
registered R3 stage statistic has been computed (unit-test smoke calculations are non-stage checks).
Read-only implementation review found that the original phrase
"subtracts the `V_t`-induced bias in the γ numerator" did not define an implementable scalar subtraction
for a direction adaptively selected from the same noisy coefficients. The following operational
definitions therefore control wherever the frozen text below is ambiguous or conflicts:

1. **Variant E is cross-fitted signed held-fold GLS.** For each frozen component fold, only its complement
   learns the adaptive unit directions `(a_f,c_f)` with the registered non-negative alternating-GLS fit.
   Holding those directions fixed, the held fold jointly profiles its own GLS intercept `w̄_f` and an
   **unconstrained signed** coefficient `γ_f`. Its closed form is the Schur-complement GLS solution using
   the held targets' empirical `V_t^-1`. The R3 point estimate is `mean_f γ_f`; A1's interval is the
   five-fold t interval around that mean. Cross-fitting, rather than a standalone trace subtraction,
   removes the adaptive direction-search bias: held noise never chooses `(a_f,c_f)`. Signed `γ_f` is an
   estimating coefficient and may be negative under the null; the structural model parameter remains
   `γ >= 0`, and positive recovery remains the alternative being certified. Variant N remains the
   non-negative in-sample adaptive fit.
2. **S2 component scales are resampled.** Form the finite empirical pool of component-average
   `tr(V_t)`, divided by its median. In every synthetic replicate, draw with replacement one factor per
   frozen component and multiply every target covariance in that component by the shared drawn factor.
   This preserves the empirical heavy tail without deterministically multiplying a component by its own
   scale a second time.
3. **S3 has an exact, verified degenerate set.** Select exactly
   `floor(0.225*n_targets + 0.5)` targets without replacement using seed `1729+55`, once before the
   replicates. S3 first applies the S2 resampling. For each selected target and replicate, its covariance
   scale is increased, if needed, until `tr(V_t) >= 4||β_t^truth||²`; Gaussian coefficient noise from
   that effective covariance is then redrawn conditional on the upstream audit statistic
   `||β̂_t||²-tr(V_t) <= 0`. Every selected target is asserted non-positive in every replicate and the
   maximum, count, and redraw total are recorded. This is a registered stress condition, not an
   unconditional Gaussian arm; S1/S2 retain the unconditional Gaussian construction.
4. **G-R3-3 uses the registered paired median.** Bootstrap replicates resample the 200 paired
   `γ̂_N-γ̂_E` values and recompute their median. Its percentile LCB95, not a bootstrap mean, must exceed
   zero.
5. **Execution validity.** Alternating GLS retains the last accepted monotone iterate if a proposed
   update is non-monotone. Convergence, iterations, and non-monotone stops are counted for E, N, and N2.
   Every replicate must produce exactly five E fold fits and one N/N2 fit, every returned fit must
   converge, the seed must equal 1729, and the replicate count must equal 200 for a certifying result;
   otherwise the run fails closed to `R3_ESTIMATOR_INSENSITIVE_NO_DECISION`. A non-200 `--replicates`
   value is a smoke run only and must truthfully control the number of generated replicates. Because one
   non-converged fit is sufficient to make G-R3-0 false, the registered run stops immediately at the
   first such fit and records all remaining replicates/regimes as not run; invalid partial estimates are
   not completed or interpreted.
6. **Scale-relative numerical stabilization.** A pre-result end-to-end smoke run exposed that empirical
   precision matrices make an absolute `1e-8 I` jitter numerically zero at the scale of the GLS normal
   equations. Every positive-semidefinite normal solve therefore uses the deterministic sequence
   `(1e-10, 1e-8, 1e-6) * max(1, max|diag(G)|)` and fails closed if all three solves fail. This changes
   only floating-point stabilization; the GLS objective, constraints, directions, gates, and data are
   unchanged.
7. **Positive-definite covariance inversion.** The same smoke audit found small numerical negative
   eigenvalues in the empirical sandwich matrices. `V_t` is first symmetrized and a candidate jitter is
   accepted only after Cholesky succeeds; its inverse is formed through the Cholesky factors and
   symmetrized. Failure through `1e-4 I` stops the run. The previous `solve`-only check could accept an
   indefinite matrix and the previous identity fallback is prohibited.

All scientific thresholds, seeds, regime names, and the A1 interval are unchanged. This amendment is
recorded before any R3 result and may not be revisited after a result.
**Program:** `task.md` Part 9 (K-LBP v2). **Stage:** R3. **Gating:** for estimator certification only.
R3 certifies that the rank-1 scalar-gated estimator can recover a known γ under realistic noise. It
reads **no real affinity label beyond what PARC M0 already read** (Metz TRAIN cells only, for the
empirical V_t); it fits no real coordinate to real labels; it authorizes no scientific claim in either
direction. A null `γ̂` on real data is only interpretable after R3 passes — this is the repair of the
disclosed PARC M0 G0 defect (`V_t = σ²I` there; the **empirical sandwich `V_t`** here).
**Design rationale:** `reports/active/model_blueprint_reconstruction_2026-07-27.md` §5, §10.

---

## 1. The single question

> **Can the Stage-2 estimator — `(w̄, γ, a, c)` fitted to per-target empirical coefficients
> `(β̂_t, V_t)` — recover a known γ* within its component-bootstrap interval, across the noise regimes
> the real substrate actually exhibits?**

The model under test (task.md §9.3, on the coefficient scale):

```text
beta_t = w_bar + gamma * (a^T k_t) * c + eps_t,   eps_t ~ N(0, V_t)
||a|| = ||c|| = 1, gamma >= 0
```

`beta_t` is the per-target ridge interaction direction (64-d ligand PCA ambient) with sandwich
covariance `V_t`, estimated by `research/hqgbma_stage_d.target_coefficients` on the exact projected
residual `M_X^W y` (the PARC M0 / Stage D machinery, unchanged).

## 2. What is read

* ChEMBL-37 Metz dense pKi panel, **TRAIN cells only**, via `research.panel_gate_pa.load_panel_train()`,
  restricted to the PARC M0 eligible set (111 targets). Only used to obtain the **empirical**
  `(V_t, n_t, signal_t)` noise distribution and the frozen component folds. No development,
  confirmation, Davis, or sealed label is read. `sealed_test_consumed=false`.
* The frozen component fold map (`research.panel_power.component_folds(DualCold.panel())`).
* A candidate coordinate `k_t` for generating structured truth: the R1 `det_proxy_card` if built, else
  the R1 `klifs_pocket_composition`; plus a `random_coordinate` arm (Gaussian, matched dimension) as
  the no-structure control. **No R2/LLM card is used.**

## 3. The estimator under test (frozen)

Alternating GLS, 200 outer iterations max, convergence at relative objective change < 1e-8, seed 1729:

1. Given `(a, c)`: `w̄` and `γ` by GLS — closed form per target block, `w̄` from
   `Σ_t V_t⁻¹ (β̂_t − γ s_t c) ` with `s_t = aᵀk_t`, and `γ = max(0, num/den)` from the same GLS.
2. Given `(w̄, γ, c)`: `a` from GLS over coordinates (normalized, `||a||=1`).
3. Given `(w̄, γ, a)`: `c` from GLS over the 64-d ambient (normalized, `||c||=1`).
4. Initialization: `a = 0` (γ = 0, the nested null) and spectral init from the leading V-weighted
   cross-covariance direction; **both** are run and the lower-objective fit retained (reported).
5. **Variant E (errors-in-variables, the registered form):** objective subtracts the `V_t`-induced
   bias in the γ numerator (the EC-Helmert-EB correction). **Variant N (naive):** no correction.
   Both are run on every synthetic replicate; the bias of N relative to E is itself a registered
   readout (it must show N biased **upward**, documenting why E is mandatory).
6. Inference: component bootstrap of `γ̂` (10,000 draws over components, seed 1729) and a per-component
   sign/signed-rank decomposition of the correction term's weighted fit improvement.

## 4. Synthetic regimes (frozen grid; truth known by construction)

For each regime: 200 replicates, seed 1729 + regime index × 1000 + replicate.

| regime | truth | noise |
| --- | --- | --- |
| S1 well-conditioned | `γ* ∈ {0.0, 0.5, 1.0}` × coordinate `det_proxy_card` | `V_t` empirical, unmodified |
| S2 heteroscedastic | `γ* = 0.5` | `V_t` empirical scaled by a component-level factor drawn from the empirical cross-component scale distribution (heavy tail preserved) |
| S3 degenerate-signal | `γ* = 0.5` | as S2, with 22.5% of targets' signal driven non-positive (matching the real fraction disclosed in the §2.8 audit) |
| S4 rank-misspecification | truth rank-2 (`γ₁*, γ₂* = 0.5, 0.3`, orthogonal directions) | as S1 — measures rank-1 attenuation (the Oracle-registered sensitivity) |
| S5 null coordinate | `γ* = 0.5` but `k_t := random_coordinate` | as S1 — false-positive rate of `γ̂` under no structure |

`a*, c*` drawn uniformly on their unit spheres per replicate (seeded); `w̄*` drawn Gaussian matched to
the empirical per-target coefficient scale; the true `β̂_t` is then generated as
`w̄* + γ*(a*ᵀ k_t)c* + chol(V_t) z_t`, `z_t ~ N(0, I)`.

## 5. Frozen gates

| gate | requirement |
| --- | --- |
| **G-R3-1 (recovery, variant E)** | In S1/S2/S3 at `γ* ≥ 0.5`: the component-bootstrap 95% CI of `γ̂` covers `γ*` in ≥ **90%** of replicates (nominal 95%, tolerance for Monte Carlo error at n=200 replicates), and median `γ̂/γ* ∈ [0.8, 1.25]` |
| **G-R3-2 (null control, variant E)** | In S1 at `γ* = 0` and in S5: the CI covers 0 in ≥ 90% of replicates; median `γ̂ ≤ 0.10` |
| **G-R3-3 (E vs N bias)** | In S2/S3: variant N's median `γ̂` exceeds variant E's by a positive amount whose bootstrap LCB95 > 0 (the upward-bias signature that motivates E) |
| **G-R3-4 (degenerate regime)** | S3 recovery degrades gracefully: coverage within 15 points of S1, not a cliff |
| **G-R3-5 (rank sensitivity, non-gating)** | S4: rank-1 `γ̂` attenuation reported (`γ̂/||γ*||_F`), expected ≈ 0.6–0.9; **non-gating**, it informs the R4 rank-2 diagnostic |

## 6. Frozen verdict rule

```text
G-R3-1 and G-R3-2 and G-R3-3 pass      -> R3_ESTIMATOR_CERTIFIED        (R4 may be attempted for surviving coordinates)
any of G-R3-1..G-R3-3 fails            -> R3_ESTIMATOR_INSENSITIVE_NO_DECISION
                                          (no R4; the estimator is redesigned and re-registered,
                                           NOT re-run on a relaxed threshold)
```

## 7. Declared expected outcome (stated before running)

Variant E recovers `γ*` within CI in the well-conditioned and heteroscedastic regimes; variant N shows
the registered upward bias, largest in S3 where 22.5% of targets carry non-positive signal. The
genuinely uncertain gate is G-R3-4: if the real heavy-tailed `V_t` distribution breaks GLS recovery
even with the correction, the estimator needs a robust variant — that would be a design finding, not a
threshold failure.

## 8. Prohibited rescues

No threshold change after a result. No regime may be dropped after seeing its outcome. No replicate
count increase after seeing results. The empirical `V_t` may not be replaced by `σ²I` anywhere in R3
(that replacement is the defect this stage repairs). No real coordinate is fitted to real labels in
R3. R3 authorizes no affinity claim.

## 9. Artifacts

```text
research/klbp_r3_synthetic.py                 runner, deterministic seed 1729
reports/active/klbp_r3.json                   machine-readable result, parses with allow_nan=False
reports/active/klbp_r3_decision.md            verdict + what was NOT shown
tests/test_klbp_r3.py                         GLS closed forms vs torch reference, V_t correctness,
                                              regime construction, seed discipline
```
