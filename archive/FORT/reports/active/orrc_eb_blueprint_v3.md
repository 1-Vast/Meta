# ORRC-EB Blueprint v3 — spectral-adaptive ORRC and error-corrected Helmert empirical Bayes

Supersedes `orrc_eb_blueprint_v2.md`. Written 2026-07-25 after the corrected PD-M stage and under
the `OPEN_DATA_ONLY` amendment. No model was trained for this revision, no threshold was relaxed, no
experimental result was altered, and no sealed or external label was read.

---

## 1. Audit of the current state

### 1.1 `task.md`

Three statements in `task.md` are now obsolete and are corrected here rather than quietly dropped.

| statement | status | correction |
|---|---|---|
| "the mathematics is closed" (v2 verdict) | **obsolete** | PD-M exposed two open items: the penalty sits at the frozen grid boundary and the spectral effective rank is 22-26, not `<= 8` |
| "effective rank at most 8" (ORRC contract, section E.1 of v2) | **withdrawn** | the figure was inherited from the old adapter's dimension, which is provably non-identifiable. Rank is an estimated quantity requiring permutation and bootstrap support |
| "`panel_davis` is the independent confirmation source" | **withdrawn** | registered and sealed, but underpowered: `MDE80 = 0.1596` against the `0.0614` reference |

### 1.2 ORRC code

`research/orrc_v2.py` and `research/panel_pdm.py` implement v2 sections A.1-A.4 correctly: weighted
projector `M_X^W`, latent as an inert observed-space decomposition, frozen `eps`-ridge for
uniqueness, no rank constraint in the objective or the selection grid, nested per-fold selection,
Holm correction. Audits: label projection relative KKT `1.17e-13`, edge-design `2.97e-15`, solver
dual feasibility `1.28e-4`, complementarity `2.46e-6`, deterministic repeat disagreement exactly
`0.0`. `research/orrc.py`, `research/panel_gate_pd0.py` and `research/panel_gate_pd1.py` remain
`SUPERSEDED -- DO NOT RUN`.

Two implementation gaps that v3 must close:

* the estimator has exactly one scalar knob (`lam_B`), and PD-M drove it to the boundary of the
  frozen grid in all five folds. A single global shrinkage level cannot simultaneously suit a
  leading singular value of `0.199` and a twentieth of `~0.02`;
* the transferable-subspace claim rests on two permutation nulls. The mandate requires nine
  qualification criteria, of which four (sampling-mask-preserving permutation, document-block
  sensitivity, bootstrap principal-angle stability, held-component predictive transfer) are not yet
  implemented.

### 1.3 PD-M result

`reports/active/panel_pdm.json`. Directions 1, 2 and 4 are feature-explainable after Holm
correction over a 16-test family, with component-bootstrap LCB95 above zero, stable to removing the
top 1% residual-energy ligands. Energy splits `0.405` feature / `0.595` latent.

**What it does not establish.** PD-M is train-only. It shows that the projected interaction residual
contains directions the frozen features explain and that transfer across held-out *training*
components. It does not show predictive superiority over B0 on any held-out evaluation set, and
under the `OPEN_DATA_ONLY` verdict no such set currently exists. Any claim of an affinity model from
PD-M or PD-H is prohibited.

---

## 2. Mathematical literature table

Every imported idea, with the audit the mandate requires. Nothing is adopted by name.

### 2.1 Adaptive weighted nuclear norm

1. **Original estimand and assumptions.** Reweighted convex relaxation (Candes, Wakin & Boyd 2008
   for `l1`; Gu et al. 2014 and Chen & Chi for the matrix case): minimise
   `sum_j omega_j sigma_j(B)` with `omega_j` inversely related to a pilot estimate of `sigma_j`, to
   approximate the rank penalty while keeping convexity for fixed weights. Assumes a pilot estimate
   independent of the data used for the final fit, and weights bounded away from zero.
2. **DTA object.** `B` is the bilinear target-by-ligand interaction coefficient on the exactly
   projected design `Z_perp`; `sigma_j` are its interaction singular values.
3. **Assumptions that hold.** Convexity for fixed weights; the projected loss is a proper quadratic;
   a pilot is available from cross-fitted training components.
4. **Assumptions that fail.** Independence of the pilot is *not* automatic — a pilot computed on the
   same rows makes the procedure non-convex in disguise and invalidates the KKT certificate. The
   usual spiked-model justification also assumes homoscedastic noise, which an incomplete panel with
   heterogeneous per-cell replication violates.
5. **Adaptation.** Weights are generated **only** from cross-fitted training information: for fold
   `f`, `omega_j^{(f)} = 1 / (sigma_j^{pilot,(-f)} + gamma)` with `gamma` frozen and the pilot fitted
   on the complement of fold `f`. The final fit then solves a convex problem with fixed weights, and
   the nuclear-norm KKT certificate remains valid with `lam * omega_j` per direction.
6. **Expected benefit.** Removes the single-knob boundary pathology: leading directions are shrunk
   less than trailing ones, so the estimator no longer has to choose one global level for a spectrum
   spanning an order of magnitude.
7. **Falsification.** If adaptive weights are only a reparametrisation of the grid, the selected
   solution and its evidence-qualified subspace will coincide with fixed-`lam` ORRC within numerical
   tolerance; and the cross-fitted pilot must beat a *shuffled* pilot (weights permuted across
   directions) on held-component transfer.
8. **Replacement baseline.** Fixed nuclear-norm ORRC at the same effective degrees of freedom.

### 2.2 Monotone singular-value shrinkage

1. **Original estimand and assumptions.** Optimal shrinkage of singular values under a spiked
   covariance model (Gavish & Donoho 2014, 2017; Donoho, Gavish & Johnstone 2018; Nadakuditi 2014
   OptShrink): for `Y = X + noise` with i.i.d. Gaussian noise and known aspect ratio, the asymptotically
   optimal estimator applies a *scalar, monotone* function to each observed singular value.
2. **DTA object.** The map from observed `sigma_j(B_hat_raw)` to a shrunk `eta(sigma_j)` before the
   subspace enters the few-shot posterior.
3. **Assumptions that hold.** Monotonicity and the constraint `0 <= eta(s) <= s` are structural, not
   distributional; a low-parameter monotone family is estimable from training components alone.
4. **Assumptions that fail.** The i.i.d. Gaussian noise model and the fixed aspect ratio do **not**
   hold: the panel is incomplete, the projection `M_X^W` induces correlated residuals, and the
   effective aspect ratio differs per fold. The closed-form Gavish-Donoho threshold is therefore
   **not** transportable and must not be used as a formula.
5. **Adaptation.** Keep only the structure. `eta_theta` is a monotone I-spline with a small frozen
   number of knots, constrained to `0 <= eta_theta(s) <= s` and non-decreasing, with `theta` fitted by
   nested component-held-out training error or a train-only marginal objective. The closed-form
   optimal shrinker is used only as a *reference curve* in reporting.
6. **Expected benefit.** Fewer parameters than a per-direction weight vector, with a shape that
   cannot invert the spectral order, and no dependence on an inapplicable asymptotic regime.
7. **Falsification.** The fitted `eta_theta` must beat both the identity (`eta(s) = s`, i.e. no
   shrinkage) and the best single soft threshold on held-out training components; if it does not, the
   extra flexibility is unjustified and the fixed-`lam` estimator stands.
8. **Replacement baseline.** Fixed hard-rank reduced-rank regression at the bootstrap-stable rank.

### 2.3 Random-effects meta-analysis with known estimation error

1. **Original estimand and assumptions.** Multivariate random-effects meta-analysis
   (DerSimonian & Laird 1986; REML): observed study effects `w_hat_t ~ N(mu, Sigma_0 + V_t)` where
   `V_t` is the *known* within-study covariance. Assumes `V_t` known (or well estimated), studies
   independent, and the between-study covariance `Sigma_0` shared.
2. **DTA object.** `w_hat_t` is the cross-fitted per-target interaction coefficient, `V_t` its
   estimation covariance, and `Sigma_0` the prior covariance the few-shot posterior needs.
3. **Assumptions that hold.** `V_t` is genuinely computable, because each `w_hat_t` comes from a
   linear solve whose covariance we can propagate; homology components give the independence unit.
4. **Assumptions that fail.** `mu = 0` is imposed rather than estimated (a non-zero prior mean is
   locked); `V_t` is estimated, not known; and targets within a homology component are not
   independent, so the unit must be the component, not the target.
5. **Adaptation.** Fit `Sigma_0` from `w_hat_t ~ N(0, Sigma_0 + V_t)` by regularised marginal
   likelihood with the factor-analytic parameterisation `Sigma_0 = L L^T + diag(exp(delta))`, using
   one `w_hat` per component (or a component-weighted likelihood), and never treating `w_hat_t` as
   noise-free.
6. **Expected benefit.** This is the precise correction the previous adapter lacked. Treating
   `w_hat_t` as noise-free inflates `Sigma_0` by the average estimation error, which makes the
   posterior over-trust support data exactly where support is weakest.
7. **Falsification.** The corrected estimator must beat the uncorrected one (which ignores `V_t`) on
   train-only held-component marginal likelihood, and a **shuffled** `V_t` (error covariances
   permuted across targets) must destroy that advantage.
8. **Replacement baseline.** Diagonal `Sigma_0`; Ledoit-Wolf/OAS shrinkage; unregularised full
   covariance.

### 2.4 Within-transformation / orthogonal deviations

1. **Original estimand and assumptions.** Panel-data fixed-effects elimination (Arellano & Bover
   1995 forward orthogonal deviations; the Helmert transform): premultiply by an orthonormal
   contrast matrix `H` with `H1 = 0` to remove an additive unit effect while preserving a spherical
   error covariance, avoiding the serial correlation that first differencing induces.
2. **DTA object.** `z_s = H_k (y_s - b_s)` removes the unknown target-level offset, including B0
   calibration error, from the support residuals.
3. **Assumptions that hold.** The offset is exactly additive within a target; `H_k H_k^T = I` keeps
   the transformed noise spherical when `D_s = sigma^2 I`.
4. **Assumptions that fail.** Heteroscedastic support noise makes `C_s = H_k D_s H_k^T` non-diagonal,
   so the transformed likelihood must carry the full `C_s`, not a scalar.
5. **Adaptation.** Retain `C_s` in full and Cholesky-solve against it; collapse duplicate support
   ligands *before* forming `H_k`, so the contrast dimension is exactly `k-1` and never `k(k-1)/2`.
6. **Expected benefit.** Ranking evidence is separated from calibration evidence by construction,
   which is exactly the confound that made the BM0/PC adapter adapt the mean instead of the order.
7. **Falsification.** A synthetic target with a pure offset and no interaction must yield `log BF = 0`
   and zero ranking correction.
8. **Replacement baseline.** Intercept-only posterior; raw (untransformed) residual posterior.

### 2.5 Factor-analytic covariance and shrinkage

1. **Original estimand and assumptions.** `Sigma = L L^T + Psi` with low-rank `L` and diagonal `Psi`
   (factor analysis; Tipping & Bishop 1999 PPCA); Ledoit & Wolf 2004 and Chen et al. 2010 (OAS)
   provide well-conditioned shrinkage targets when `n` is comparable to `p`.
2. **DTA object.** The prior covariance `Sigma_0` over interaction coefficients, at the estimated
   dimension `d ~ 22-26` with roughly 100 component-level observations.
3. **Assumptions that hold.** Positive definiteness by construction; the regime `n ~ 100`, `d ~ 25`
   is exactly where a full unstructured covariance is unstable and a factor model is standard.
4. **Assumptions that fail.** The factor rank is unknown and must not be chosen on external data;
   PPCA assumes homoscedastic residual variance, which `diag(exp(delta))` relaxes.
5. **Adaptation.** Factor rank selected by nested training components only; `L`, `delta` fitted by
   regularised marginal likelihood with an eigenvalue floor and a condition-number cap.
6. **Expected benefit.** A `25 x 25` covariance from ~100 observations is otherwise ill-conditioned;
   this is the concrete blocker v2 Amendment 2 flagged.
7. **Falsification.** The factor model must beat diagonal and Ledoit-Wolf/OAS baselines on train-only
   held-component marginal likelihood.
8. **Replacement baseline.** Diagonal; Ledoit-Wolf; OAS; unregularised full.

### 2.6 Neyman-orthogonal / partialled-out estimation

1. **Original estimand and assumptions.** Partialling out nuisance parameters so the target estimator
   is first-order insensitive to nuisance error (Frisch-Waugh-Lovell; Neyman orthogonality;
   Chernozhukov et al. 2018 double machine learning). Assumes the nuisance is correctly specified as
   a linear projection, and cross-fitting to avoid own-observation bias.
2. **DTA object.** `M_X^W` removes target and ligand main effects exactly, so the interaction
   coefficient is orthogonal to them in the observed edge space.
3. **Assumptions that hold.** The nuisance is a genuinely linear incidence design, so the projection
   is exact rather than estimated by a learner.
4. **Assumptions that fail.** Exactness holds only *in-sample*; a held-out target's main effect is
   unavailable at deployment, which is why the readout permits a label-blind per-target constant.
5. **Adaptation.** Already implemented: projection for estimation, product readout with training-only
   reference measures for prediction, with the algebraic identity `M_X^W Z_centered = M_X^W Z`
   guaranteeing the two describe one estimand.
6. **Expected benefit.** No soft centering or correlation penalty is needed, so no real interaction
   is suppressed by a penalty term (the CFRI defect).
7. **Falsification.** `max |X^T W M_X^W v|` must vanish to `1e-8` relative, and a pure main-effect
   vector must project to zero.
8. **Replacement baseline.** Soft-centred CFRI (retained only as a control arm).

---

## 3. Innovation 1 — OSA-ORRC (Orthogonal Spectral-Adaptive ORRC): preregistration

**Preserved unchanged.** `r = M_X^W y`, `Z_perp = M_X^W Z`, exact observed-edge nuisance
orthogonality, the inert latent decomposition, the frozen product readout with training-only
reference measures, and the `eps`-ridge uniqueness convention. Soft centering and correlation
penalties stay forbidden.

**Replaced.** The single fixed nuclear-norm level and the inherited rank convention.

### 3.1 Two candidates, both selected on training components only

```text
(A)  min_B  0.5 ||W^(1/2)(r - Z_perp vec(B))||_2^2 + lambda * sum_j omega_j sigma_j(B) + (eps/2)||B||_F^2
     omega_j^(f) = 1 / (sigma_j^{pilot,(-f)} + gamma),  pilot fitted on the complement of fold f

(B)  B_hat = P diag(eta_theta(sigma_j)) Q^T   from the unshrunk projected solution P diag(sigma) Q^T
     subject to  0 <= eta_theta(s) <= s  and  eta_theta non-decreasing in s
     eta_theta a monotone I-spline with a frozen knot set; theta by nested component-held-out error
```

Candidate choice is made on training components only and **never** on any external panel. The frozen
grids are `lambda in {0.50, 0.25, 0.10, 0.05, 0.02, 0.01} x lambda_max` (widened *before* any run of
this specification, and registered here explicitly because PD-M hit the old boundary) and
`gamma in {0.1, 0.25, 0.5} x sigma_1^{pilot}`. If the selected `lambda` again lands on a boundary,
that is reported as `OSA_ORRC_GRID_BOUNDARY_STOP` and the grid is **not** widened afterwards.

### 3.2 Four ranks, reported separately

| rank | definition |
|---|---|
| algebraic | `#{ sigma_j > 1e-10 }` |
| spectral effective | `#{ sigma_j >= 0.05 sigma_1 }` |
| bootstrap-stable | largest `d` whose leading-`d` subspace has component-bootstrap principal-angle stability above the frozen bound |
| **evidence-qualified transferable** | directions passing all nine criteria in 3.3 |

Only the evidence-qualified transferable subspace may enter the few-shot posterior.

### 3.3 Evidence qualification — all nine required

1. exposure-matched target-feature permutation;
2. exposure-matched ligand-feature permutation;
3. **sampling-mask-preserving permutation** (permute values within the observed mask so the
   missingness pattern is held fixed);
4. **publication/document-block sensitivity** (leave-document-out refit; the direction must survive);
5. multiplicity correction over the full frozen family (Holm across all directions x all nulls);
6. homology-component bootstrap LCB95 above zero;
7. top-energy-ligand removal;
8. **bootstrap principal-angle stability** of the direction across component bootstraps;
9. **held-component predictive transfer** of that direction alone.

### 3.4 Required baselines, parameter-matched

fixed nuclear-norm ORRC; fixed hard-rank reduced-rank regression; ordinary ridge bilinear
regression; random target features; random ligand features; unconstrained interaction control (I0);
ligand-only B0.

### 3.5 Stop verdicts

`OSA_ORRC_GRID_BOUNDARY_STOP` if the penalty again selects a boundary point;
`OSA_ORRC_NO_STABLE_SUBSPACE_STOP` if no direction passes all nine criteria. Neither may be followed
by widening a grid or dropping a criterion.

---

## 4. Innovation 2 — EC-Helmert-EB: preregistration

### 4.1 Model

```text
z_s = H_k (y_s - b_s),   X_s = H_k Psi_s,   H_k 1 = 0,  H_k H_k^T = I_{k-1}
z_s | w_t ~ N(X_s w_t, C_s),      C_s = H_k D_s H_k^T
w_t       ~ N(0, Sigma_0)                       (zero mean; a target-conditioned mean stays locked)
Sigma_post = (Sigma_0^{-1} + X_s^T C_s^{-1} X_s)^{-1}
m_post     = Sigma_post X_s^T C_s^{-1} z_s
```

`D_s = diag(sigma_meas^2(i) + v_base(i))` separates train-only measurement variance from the base's
own predictive variance. Duplicate support ligands (identical parent connectivity or Tanimoto
`>= 0.95`) are collapsed **before** `H_k` is formed.

### 4.2 The error correction

For every cross-fitted training target retain `w_hat_t`, its estimation covariance `V_t`, and its
component and fold provenance. Fit

```text
w_hat_t ~ N(0, Sigma_0 + V_t),      Sigma_0 = L L^T + diag(exp(delta))
```

by regularised marginal likelihood, one observation per homology component, factor rank chosen by
nested training components. Positive definiteness is structural; an eigenvalue floor and a
condition-number cap of `1e6` are applied and the raw condition number is reported. Treating
`w_hat_t` as noise-free is exactly the error this module exists to remove.

### 4.3 Required covariance controls

diagonal; Ledoit-Wolf; OAS; uncorrected (ignores `V_t`); unregularised full; low-rank-plus-diagonal;
shuffled coefficient-error covariance.

### 4.4 Evidence and invariants

```text
M0: z_s ~ N(0, C_s)          M1: z_s ~ N(0, C_s + X_s Sigma_0 X_s^T)
log BF = log p(z|M1) - log p(z|M0),   pi_t = rho BF / (1 + rho BF),  rho = 1 frozen
```

Mandatory, each an executable test: exact B0 mean at `k = 0`; no ranking correction at `k <= 1`;
support-permutation invariance; duplicate collapse before contrasts; no support label entering any
encoder; separate intercept and ranking posteriors, with the intercept excluded from `log BF`;
evidence-based adaptation with no learned free-form gate; finite positive posterior covariance;
epistemic contraction under repeated consistent evidence; wrong-target and label-permuted support
falsifications; FP32 Cholesky solves throughout.

### 4.5 Trainable global parameters

The system must not collapse to a fixed handcrafted score. Trainable globally: the ligand base `b(d)`
and its variance head, the shared interaction feature maps, the spectral shrinker parameters
`theta`, and the covariance parameters `(L, delta)`. Target-specific inference remains exact and
closed-form. This is checked by a test asserting a non-empty trainable-parameter set whose
perturbation changes predictions.

---

## 5. Gate sequence and its current state

| gate | content | state |
|---|---|---|
| OPEN-S0 | licence, availability, provenance, reproducibility | **done** |
| OPEN-S1 | endpoint, missingness, overlap, firewall, power shape | **done — failed for every candidate** |
| PD-M2 | train-only spectral regularisation and transferable-subspace closure | specified, admissible, not run |
| PD-H | train-only error-corrected covariance and Helmert closure | specified, admissible, not run |
| OPEN-P | arm-blind external panel power audit | **blocked: no admissible panel** |
| one-seed external gate | first predictive claim | **blocked** |
| three-seed | — | blocked |
| final single-use confirmation | — | blocked |

---

## 6. Verdict

```text
NO_OPEN_POWERED_INDEPENDENT_PANEL
```

The mathematics of v3 is specified and falsifiable, and the two innovations are within the budget.
But the mandate's own rule applies before any of it can be tested predictively: no open, adequately
powered, independent panel exists. ChEMBL is training-only and its Metz development rows are spent;
`panel_davis` is registered, sealed and underpowered by a factor of `2.6`; BindingDB-native curation
yields 38 targets with a `>= 40`-ligand assay-controlled block against a requirement of ~100
components, before any firewall; the Klaeger `Kdapp` matrix is not openly recoverable; PKIS2 is
percent-inhibition and barred from continuous-affinity confirmation.

PD-M2 and PD-H may proceed as train-only mathematical development. **No predictive claim about
dual-cold affinity ranking is available from this program until an open, powered, independent panel
exists**, and none may be manufactured from train-only results.
