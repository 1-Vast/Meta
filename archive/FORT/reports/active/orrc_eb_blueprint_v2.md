# ORRC-EB blueprint, revision 2 — corrected mathematical and governance specification

Written 2026-07-25 under the audit verdicts `BLUEPRINT_REQUIRES_MAJOR_REVISION` (mathematics and
architecture) and `INSUFFICIENT_FILES_FOR_EMPIRICAL_VERDICT` (empirics). No training was run for this
revision, no threshold was changed, no experimental result was modified, and no confirmation or
sealed label was read. Gate PB remains a preregistered failure and Gate PC remains exploratory; this
document does not reinterpret either as positive.

**Status of the halted implementation.** A pre-revision Gate PD0 run had already been launched when
the audit arrived and it finished at 19:07 before the process could be stopped. Its artifact
`reports/active/panel_gate_pd0.json` is preserved byte-identically and is marked **VOID** by
`reports/active/panel_gate_pd0_status.md`: it implements the *superseded* contract (unweighted
projector, full-grid two-sided latent constraint, rank-restricted candidate grid) and therefore
cannot be a PD0 result under this revision. Its statistics are not used, cited or interpreted
anywhere below. One thing from it is retained, as a bug diagnosis rather than a result: the fitted
auxiliary-feature energy fell monotonically as the latent penalty was relaxed, which is exactly the
leakage that section C proves the old constraint permits.

---

## A. Corrected mathematical specification

### A.1 Edge space, weights and the weighted nuisance projector

Let `Omega` be the observed target-ligand cell set of one endpoint-consistent panel, `n = |Omega|`,
and index edges `e = (t,l)`. Fix a diagonal weight matrix `W = diag(w_e)`, `w_e > 0`, chosen
label-free and frozen before any fit. The inner product and norm are

```text
<a,b>_W = a^T W b ,      ||a||_W^2 = a^T W a .
```

Default `W = I`. A non-uniform `W` is admissible only if it is a *validated* inverse-variance model;
the record count `n_records` is explicitly not such a model (PA and HOAB registered unit weights for
this reason), so any non-uniform choice requires its own train-only reliability audit and must be
registered before use.

Let `X in R^{n x (T+L)}` be the target/ligand incidence design (one indicator block per axis). The
weighted nuisance projector onto the `W`-orthogonal complement of `col(X)` is

```text
M_X^W = I - X (X^T W X)^+ X^T W .
```

`M_X^W` is idempotent and `W`-self-adjoint (`<M a, b>_W = <a, M b>_W`); with `W = I` it reduces to
the exact Hodge projector Gate PA audited. The interaction residual and design are

```text
r      = M_X^W y
Z      in R^{n x pq},   row_e(Z) = vec(u_t v_l^T)
Z_perp = M_X^W Z .
```

Centering the features first changes nothing: `M_X^W Z_centered = M_X^W Z`, because centering shifts
`Z` by columns that lie in `col(X)`. This is the algebraic reason the training estimand and the
out-of-sample readout can be unified in section D.

### A.2 The convex estimation problem

```text
minimise over B in R^{p x q}:
    F(B) = 0.5 || r - Z_perp vec(B) ||_W^2 + lam_B ||B||_*  + (eps/2) ||B||_F^2
```

* squared loss in the `W` inner product, nuclear norm for reduced rank, and a frozen ridge
  `eps = 1e-8 * lam_B` whose only purpose is uniqueness (section C.3);
* there is **no** hard rank constraint. `rank(B) <= 8` is not imposed, because it is not convex.
  Effective rank is a post-hoc frozen spectral rule (section A.4);
* `lam_B` is selected only by nested train-component cross-validation on a frozen scale-free grid
  `lam_B in {0.50, 0.25, 0.10, 0.05} x lam_B_max`, `lam_B_max = ||U^T diag(W r) V||_op` (the smallest
  penalty for which `B = 0` is optimal). The grid is **not** filtered by rank; filtering it by rank
  would smuggle the non-convex constraint back in through the selection rule.

### A.3 The latent block, corrected

The latent term is an **observed-edge variable** `l in R^n`, never a full matrix and never an
imputed cell. It is required to be orthogonal in the *observed* space to both the nuisance design and
the interaction design:

```text
(C1)  X^T W l      = 0
(C2)  Z_perp^T W l = 0
```

Both are affine (linear) constraints, so the problem remains convex. Because `(C1)` implies
`M_X^W l = l`, the pair `(C1) and (C2)` is equivalent to the single statement

```text
l in range(Pi),    Pi = I - A (A^T W A)^+ A^T W ,   A = [X , Z] .
```

Row- and column-space orthogonality of a full-grid latent matrix (`U^T L = 0`, `L V = 0`) is **not**
sufficient and is withdrawn; section C.1 gives an explicit counterexample on a 3x3 panel with one
missing cell.

**Consequence (inertness).** Under `(C1)-(C2)` the residual decomposes `W`-orthogonally as
`r = P_{Z_perp} r + Pi r`, the objective separates, and the estimate of `B` is *identical* whether or
not a latent block is fitted, for every latent penalty. The latent block therefore cannot bias, mask
or inflate the feature coefficient. It is retained only as a reporting object:

```text
feature energy  = ||P_{Z_perp} r||_W^2 / ||r||_W^2
latent energy   = ||Pi r||_W^2         / ||r||_W^2
```

and, if a low-rank summary of the latent part is wanted, it is computed *after* `B` by a nuclear-norm
fit on `Pi r` alone. It never enters any prediction, for any target or ligand, seen or unseen.

### A.4 Effective rank, frozen after solving

```text
effective_rank(B) = # { i : sigma_i(B) >= 0.05 * sigma_1(B) }
```

applied to the solved coefficient. The numerical rank (`sigma_i > 1e-10`) is reported beside it. At
most the leading `8` directions are *tested*, which is a multiplicity budget, not a model constraint;
per-direction p-values are Holm-corrected across the tested directions. The rank claim must be
supported by permutation and component-bootstrap evidence per direction and may not be inherited
from the dimension of the previous Bayesian adapter.

### A.5 Numerical tolerances (frozen)

| quantity | tolerance |
|---|---|
| projection relative KKT `max\|X^T W M_X^W v\| / max\|X^T W v\|` | `< 1e-8` |
| projection idempotence `\|\|M^2 v - M v\|\| / \|\|M v\|\|` | `< 1e-7` |
| independent solver disagreement (LSMR vs LSQR vs dense pseudo-inverse) | `< 1e-6` |
| primal feasibility of `(C1)`,`(C2)`: `max\|A^T W l\| / max\|l\|` | `<= 1e-8` |
| nuclear-norm dual feasibility `max(0, \|\|grad\|\|_op / lam - 1)` | `<= 1e-3` |
| complementarity `\|<grad,B> - lam\|\|B\|\|_*\| / (1 + lam\|\|B\|\|_*)` | `<= 1e-3` |
| deterministic repeat disagreement (bitwise-capable path) | `<= 1e-6` |
| empirical-Bayes covariance condition number | `<= 1e6` |

---

## B. Assumptions required for identifiability

1. `Omega` is fixed and label-blind; no cell is added, removed or imputed on the basis of a label.
2. `w_e > 0`, frozen, label-free. Non-uniform weights require a registered variance model.
3. The target and ligand feature maps are frozen before any fit, fitted on training entities only,
   and no affinity label selects a basis or a dimension.
4. `eps > 0` in A.2. Without it the nuclear-norm minimiser may be a set rather than a point.
5. `B` is reported in `N^perp`, the orthogonal complement of `null(Z_perp)` (section D.3); the ridge
   in A.2 selects exactly this representative.
6. Reference measures `u_bar_W`, `v_bar_W` are computed from training edges only (section D.2).
7. The latent block satisfies `(C1)-(C2)`; otherwise the split is determined by the penalty ratio
   rather than by the data (section C.1).
8. For the empirical-Bayes stage: `D_s` is SPD, `Sigma_{0,t}` is SPD after the eigenvalue floor, and
   `k >= 2` so that at least one Helmert contrast exists.
9. Support rows are de-duplicated (section E.6); the contrast space has dimension exactly `k-1` and
   is never treated as `k(k-1)/2` independent pair differences.
10. Every statistic is aggregated to the sequence-homology component before inference; edges,
    ligands and episodes are never resampling units.

---

## C. Uniqueness of the `B` / `L` decomposition: counterexample and proof

### C.1 Counterexample — full-grid feature orthogonality is not observed-space orthogonality

Take 3 targets with feature `u = (1, 0, -1)`, 3 ligands with `v = (1, 0, -1)`, `U = [1, u]`,
`V = [1, v]`, and the latent matrix `L = a b^T` with `a = b = (1, -2, 1)`, which satisfies the old
constraint exactly:

```text
max |U^T L| = 0 ,   max |L V| = 0 ,   <Z, L>_full-grid = 0 .
```

Now delete the single cell `(1,1)`, leaving `n = 8` observed edges and unit weights. Then, for
`l = P_Omega(L)`:

```text
X^T l          = (-1, 0, 0, -1, 0, 0)^T   (not 0)
<Z_perp, l>_W  = -2.25
cosine(Z_perp, M_X^W l) = -0.293
```

One missing cell out of nine is enough to align the "orthogonal" latent block 29% with the
interaction direction. On the real panel, 81% of cells are missing. Hence `(C2)` is necessary, and
any decomposition claiming orthogonality from full-grid row/column constraints is not identifiable.
(Verification script content is reproduced in `research/orrc.py` tests; the numbers above are exact.)

### C.2 Proposition — inertness of the latent block under `(C1)-(C2)`

Let `S = range(Pi)` as in A.3. For `l in S`, `<Z_perp vec(B), l>_W = 0` for every `B`. Therefore

```text
|| r - Z_perp vec(B) - l ||_W^2
    = || P_{Z_perp} r - Z_perp vec(B) ||_W^2 + || Pi r - l ||_W^2 + const ,
```

the objective is separable, and `argmin_B` does not depend on `l` or on the latent penalty. The
converse also holds: if `(C2)` is dropped, the cross term `-2 <Z_perp vec(B), l>_W` is generally
non-zero (C.1), the two blocks compete, and the split is set by `lam_B / lam_L` rather than by the
data. This is precisely the behaviour observed in the void pre-revision grid.

### C.3 Proposition — existence and uniqueness of `B`

`F` in A.2 is the sum of a convex quadratic, a convex norm and a strictly convex quadratic
(`eps > 0`), hence strictly convex and coercive on `R^{p x q}`; the minimiser exists and is unique.
Its fitted vector `Z_perp vec(B)` is unique even at `eps = 0`, because the quadratic is strictly
convex in the fitted vector. For `eps = 0` the argmin may be a convex set whenever `null(Z_perp)`
intersects the subdifferential geometry of the nuclear norm; the ridge removes this and selects the
minimum-Frobenius-norm element, which lies in `N^perp`. Optimality is certified by the two frozen
KKT residuals in A.5, which are the exact conditions
`||grad||_op <= lam_B` and `<grad, B> = lam_B ||B||_*` for the nuclear norm.

### C.4 Proposition — the current adapter is non-identifiable and cannot re-orient

With `phi_j(t,d) = s_j(t) z_j(d)` and `w_j ~ N(0, lambda_j(t)^{-1})`, write `S = diag(s(t))`,
`Lambda = diag(lambda(t))`, `Phi_s = Z_s S`. The posterior mean correction at a query is

```text
g_q = z_q^T S (Lambda + S Z_s^T Z_s S / sigma^2)^{-1} S Z_s^T r / sigma^2
    = z_q^T (diag(1/v) + Z_s^T Z_s / sigma^2)^{-1} Z_s^T r / sigma^2 ,
      with  v_j = s_j(t)^2 / lambda_j(t) ,
```

and the model evidence is `N(0, sigma^2 I + Z_s diag(v) Z_s^T)`. Both depend on `(s, lambda)` only
through `v`. The gate and the precision are therefore **not separately identifiable**: the network
has two positive pathways parametrising one positive effective-variance vector.

Furthermore, if the support design is `C`-orthogonal (`Z_s^T Z_s` diagonal, e.g. after whitening),
then `m_j = (1/v_j + a_j)^{-1} b_j` with `a_j > 0`, so `sign(m_j) = sign(b_j)` for **every** positive
`v`. A positive, coordinate-wise target conditioning can shrink coordinates but cannot flip a sign or
rotate the basis — it cannot change *which way* the correction points. This is a structural
explanation of the Gate PC outcome and the reason a **full covariance** `Sigma_{0,t}` (section E) is
the identifiable and expressive replacement. The existing `model/bayesian_meta.py` is reclassified as
a positive-rescaling ablation; the identifiable normalisation for any future use is `s(t) === 1`,
with all target conditioning carried by the covariance.

---

## D. The out-of-sample estimand

### D.1 One model, one estimand

```text
y(t,l) = a_t + b_l + g(t,l) + eps ,      g(t,l) = (u_t - u_bar_W)^T B (v_l - v_bar_W)
```

with `a_t`, `b_l` free nuisance parameters. Profiling `(a, b)` out by weighted least squares is
*exactly* applying `M_X^W`, and by A.1 the projected design is the same for raw and centered
features. Training and prediction therefore share one estimand: `B` is the coefficient of the
bilinear interaction term in this model, and the readout is that same term.

### D.2 Training-only reference measures

```text
u_bar_W = sum_{e in Omega} w_e u_{t(e)} / sum_e w_e
v_bar_W = sum_{e in Omega} w_e v_{l(e)} / sum_e w_e
```

These are edge-weighted marginals of the *training* design, not unweighted entity means, and they
are frozen with the fit. Under the induced reference product measure, `g` is exactly double-centered:
`E_t[g(t,l)] = 0` for every `l` and `E_l[g(t,l)] = 0` for every `t`. That is the formal sense in
which the readout is an interaction and contains no target-only or ligand-only component.

### D.3 Null space and normalisation

`N = null(Z_perp) = { vec(B) : Z B-image lies in col(X) on Omega }`. Elements of `N` are coefficient
matrices whose bilinear function is additively decomposable *on the observed edges* but not
necessarily off them — they are exactly the directions that training cannot see and that would
extrapolate arbitrarily. The reported coefficient is normalised to `N^perp`; the ridge of A.2
performs this automatically. Any component of a candidate `B` in `N` is set to zero and its
Frobenius weight is reported as `null_space_energy`.

### D.4 Behaviour outside the training feature span

Target coordinates are PCA scores in a basis fitted on training targets, so an unseen target has no
component outside the basis by construction; only its coordinate magnitude can be unusual. Ligand
coordinates use a map fitted on training ligands and *applied* to unseen ligands. Policy:

* no clipping, no re-fitting, no re-centering on the query set — any of these would be a tuned
  nonlinearity introduced after seeing the query distribution;
* the leverage `||v_l - v_bar_W||^2` and the fraction of query ligands beyond the training coordinate
  range are **reported** with every evaluation;
* a registered pre-specified extrapolation stratum (top decile of leverage) is reported separately,
  and is descriptive, never a gate.

### D.5 Query centering and candidate-set dependence

A per-target constant may be subtracted from `g` over that target's fixed query set. It is
label-blind, and within-target Spearman is exactly invariant to it, so it cannot manufacture ranking.
It does change RMSE, therefore: the primary ranking metric is reported without and with centering
(numerically identical by construction, and asserted in a test), and RMSE is reported **both** ways.
Predictions may depend on the *identity* of the query candidate set only through this centering
constant; no other transductive dependence is permitted, and no development label may enter it.

---

## E. Helmert-contrast empirical Bayes: complete specification

### E.1 Transferable interaction coordinates

From the frozen fit `B = sum_i sigma_i a_i c_i^T`, retain the `d` directions that survived the
train-only direction audit. For any ligand and target

```text
psi_l = ( sigma_i * c_i^T (v_l - v_bar_W) )_{i=1..d}      (computable for any ligand)
c_t   = ( a_i^T (u_t - u_bar_W) )_{i=1..d}                (unknown for a new target)
g(t,l) = c_t^T psi_l
```

The few-shot problem is inference of `c_t` from `k` support measurements, with `psi` frozen.

### E.2 Contrast construction

`H_k in R^{(k-1) x k}` is any fixed orthonormal Helmert matrix: `H_k 1 = 0`, `H_k H_k^T = I_{k-1}`.

```text
z_s = H_k (y_s - b_s) ,     X_s = H_k Psi_s in R^{(k-1) x d}
```

`H_k` annihilates the target-level offset, so `z_s` contains ranking information only and the base's
absolute calibration error cannot leak into it.

### E.3 Noise model

```text
D_s = diag( sigma_meas^2(i) + v_base(i) ) ,     C_s = H_k D_s H_k^T
```

`sigma_meas^2` comes from the train-only replicate-spread audit of the panel (never from query
outcomes); `v_base(i)` is the ligand-only base's own predictive variance at that support ligand. This
is the required separation of measurement noise from base uncertainty. Homoscedastic
`D_s = sigma^2 I` is the default until the reliability audit supports otherwise; both are frozen
before use. `C_s` is SPD whenever `D_s` is.

### E.4 Prior and posterior

```text
c_t ~ N(0, Sigma_{0,t})                                   (zero mean, locked)
Sigma_post = ( Sigma_{0,t}^{-1} + X_s^T C_s^{-1} X_s )^{-1}
m_post     = Sigma_post X_s^T C_s^{-1} z_s
g_hat(q)   = pi_t * psi_q^T m_post
```

### E.5 Covariance estimator (the identifiable replacement for the positive gate)

```text
Sigma_{0,t} = sum_j kappa(t,j) c_hat_j c_hat_j^T + eta I ,     kappa(t,j) >= 0, sum_j kappa = 1
```

* `c_hat_j` are **cross-fitted** coefficient vectors of *training* targets: `c_hat_j` comes from a
  fit whose training set excluded `j`'s entire homology component. No held-out target coefficient and
  no development label may enter.
* `kappa(t,j) = softmax_j( cos(ESM(t), ESM(j)) / tau )` uses label-blind ESM similarity only, with
  the bandwidth `tau` fixed by the median-similarity heuristic on training targets and never tuned.
* eigenvalue floor `eta = max(1e-6, 1e-3 * tr(Sigma_bar)/d)`; eigenvalues are then clipped so the
  condition number is at most `1e6`, and the raw condition number is reported.
* the estimator is a second moment, so the prior mean stays exactly zero while the *orientation* is
  target-conditioned. This is what a positive diagonal cannot do (C.4).

### E.6 Evidence, abstention and duplicates

```text
M0 : z_s ~ N( 0, C_s )
M1 : z_s ~ N( 0, C_s + X_s Sigma_{0,t} X_s^T )
log BF = log p(z_s | M1) - log p(z_s | M0)
pi_t   = rho * BF / (1 + rho * BF) ,      prior odds rho = 1 (frozen, not tuned)
```

* **k = 0 and k = 1**: no contrast exists (`k-1 = 0`), so `pi_t = 0` by construction and the
  prediction is *exactly* the B0 predictive mean. Hard abstention, asserted bitwise in tests.
* **duplicates**: support ligands with identical parent connectivity or Tanimoto `>= 0.95` are
  collapsed to one row before `H_k` is formed, with the averaged label and noise scaled by `1/m`.
* **the scalar intercept posterior** (a separate `N(0, tau^2)` model on the support residual mean
  with noise `sigma^2/k`) adjusts absolute calibration only. It is reported and ablated separately,
  and it never enters `log BF` or the ranking prediction.

---

## F. Revised gate sequence — development separated from confirmation

The previous sequence PD0 -> PD1 -> PD2 is withdrawn because PD1 would have been scored on the Gate
PB development rows, and ORRC-EB was designed *after* those rows were observed. Those rows can no
longer provide independent confirmation for this route.

| stage | status | data | may authorise |
|---|---|---|---|
| **PD-M** method development | train-only, non-gating | panel TRAIN cells only | nothing outside itself |
| **PD-X** exploratory check | exploratory, labelled | panel development rows | hypotheses only |
| **PD-C** confirmation | preregistered gate | a newly registered independent panel, or approved sealed access | multi-seed review |

* **PD-M** implements section A, runs the numerical audits of A.5, and reports the train-only
  direction audit (component-level statistics, graph-exposure-matched permutations on both axes,
  Holm correction, component bootstrap, top-1% ligand stability). It produces no development
  prediction and no verdict that authorises a predictive gate. Its purpose is to establish that the
  mathematics behaves as specified on real data.
* **PD-X** may be run only after PD-M, must be titled `EXPLORATORY` in its report and in
  `history.md`, and its numbers may never be quoted as confirmation, may never authorise a later
  gate, and may never be used to select any hyper-parameter, threshold or direction count.
* **PD-C** requires an evaluation source that this route has never seen. Two admissible sources:
  (i) a newly registered endpoint-consistent panel — the pKd candidates `CHEMBL1908390` (179x70 core)
  and `CHEMBL3991601` (101x156 core) are *candidates only*; each must pass its own registry audit,
  keep pKi and pKd separate with panel-specific intercepts, ligand bases, noise scales and endpoint
  scales, and demonstrate its own component count and power before any label is read; or
  (ii) separately approved access to a genuinely untouched confirmation set. Neither is registered
  today.
* **Power unit**: the sequence-homology component, always. The threshold for PD-C is
  `max(0.03, MDE80)` where MDE80 is computed from *arm heterogeneity* on the new source, not from
  same-arm retraining noise. `reports/active/panel_power_pd1.json` remains a valid arithmetic record
  of arm heterogeneity in the PB contrasts (`0.0614`), but the gate it was written for is withdrawn,
  so it is a reference value and not a live threshold.

---

## G. Minimal code and tests required once the mathematics closes

Only after PD-C's evaluation source is registered. Modules:

1. `research/orrc.py` — weighted projector `M_X^W`, edge design, observed-edge latent projector `Pi`,
   convex solver with the `eps` ridge, KKT/feasibility audits, frozen readout.
2. `research/orrc_contrast.py` — Helmert construction, `C_s`, cross-fitted `Sigma_{0,t}`, posterior,
   Bayes factor, hard abstention, duplicate collapsing, separate intercept posterior.
3. `research/panel_pdm.py` — the train-only development stage of section F.

Tests that must pass before any stage is run:

* projection: exact main-effect annihilation, weighted-KKT, idempotence, three-solver agreement;
* **the C.1 counterexample as a regression test**: full-grid two-sided orthogonality must be shown
  *insufficient*, and `(C1)-(C2)` sufficient, on the 3x3-minus-one-cell panel;
* **inertness**: `B` identical with and without the latent block, for at least three latent penalties;
* uniqueness: two solver runs agree bitwise; `B` has zero component in `null(Z_perp)`;
* readout: `g` is double-centered under the training reference measures; Spearman is invariant to
  query centering; the ligand map is unchanged by adding unseen molecules;
* adapter identifiability: two `(s, lambda)` pairs with equal `s^2/lambda` give bitwise-identical
  predictions and evidence (the C.4 statement as an executable test);
* contrast layer: `H_k 1 = 0`, `H_k H_k^T = I`; `k <= 1` returns the B0 mean bitwise; permuting
  support rows leaves the posterior unchanged; duplicated support does not increase evidence;
  `Sigma_{0,t}` is SPD with condition number within the frozen limit; the intercept posterior does
  not change `log BF`.

---

## H. Verdict

```text
ORRC_BLUEPRINT_REQUIRES_REVISION
```

Sections A–E close the mathematics: the observed-space decomposition is specified and proved
identifiable (A.3, C.1–C.3), the out-of-sample mapping is unified with the training estimand and
normalised (D), the covariance estimator is fully specified and SPD by construction (E.5), and the
power unit is fixed at the homology component (F). Item 4 of the audit is satisfied by removing the
rank constraint from both the objective and the selection grid.

The single remaining gap is the one the verdict rule makes decisive: **the independent evaluation
source is specified only as a requirement set, not as a registered, power-verified panel.** Until a
new endpoint-consistent panel is built, audited and power-checked — or untouched confirmation access
is separately approved — no predictive gate for this route is admissible, so `ORRC_PD0_READY` cannot
honestly be claimed. The next deliverable is therefore the registration of that source, not any
model run.

Locked for the duration: signed target-conditioned prior mean, Hierarchical MoT, multi-seed runs,
long training, confirmation and sealed evaluation.

---

## Amendment 1, 2026-07-25 — the independent evaluation source was searched for and is not available

Section F's open item was closed procedurally and answered negatively.
`reports/active/panel_davis_registration.md` registers `CHEMBL1908390` (Davis 2011 pKd) as a
confirmation candidate: 102 independent homology components against the development panel's 101, a
different publication, endpoint and assay platform, no censoring floor, zero overlap with the main
registry's development or confirmation cells, and no ORRC-route gate has ever read one of its
labels. It is built by `tools/panel_registry.py`, which `tests/test_panel_registry_equivalence.py`
proves reproduces the frozen Metz registry cell-for-cell.

The pre-registered arm-blind power audit then failed its own rule. Davis has 12 query ligands per
target against Metz's 43, so its per-component retraining spread is `0.2314` against `0.0948` and its
retraining MDE80 is `0.0688`. Under the rule frozen before the audit ran,
`MDE80_PD-C = max(0.03, 2.32 x 0.0688) = 0.1596`, far above the `0.0614` arm-heterogeneity reference.
**No confirmation gate is run on Davis**; the panel stays sealed with `consumed=false`.

Pooling every remaining candidate in the local extract by inverse variance reaches only
`MDE80 ~ 0.082`. An admissible source needs at least ~100 independent components **and** at least ~40
query ligands per target after the scaffold-disjoint split — a second Metz-scale panel, which the
frozen extract does not contain.

Consequences for this blueprint:

* section F stage **PD-C remains unavailable** for want of a source, not for want of a specification;
* stage **PD-M** (train-only method development on the Metz panel) is unaffected and remains the only
  authorized next execution step, because it produces no development prediction;
* the verdict below is unchanged at `ORRC_BLUEPRINT_REQUIRES_REVISION`, and the blocking item is now
  quantified rather than open-ended: acquiring or constructing a Metz-scale independent panel.

---

## Amendment 2, 2026-07-25 — stage PD-M executed; the corrected mathematics recovers transferable directions

`research/orrc_v2.py` + `research/panel_pdm.py` implement sections A.1-A.4 and run the train-only,
non-gating development stage. Report: `reports/active/panel_pdm.json`.

Audits all pass: label projection relative KKT `1.17e-13`, edge-design projection `2.97e-15`, solver
dual feasibility `1.28e-4` and complementarity `2.46e-6` against the frozen `1e-3`, deterministic
repeat disagreement exactly `0.0`. Exact residual energy decomposition: feature `0.4052`, latent
`0.5948`.

After 4,096 exposure-matched permutations on both axes and Holm correction across 16 tests, three of
eight tested directions are feature-explainable with component-bootstrap `LCB95 > 0` (directions 1,
2 and 4: mean component correlations `0.0984`, `0.0758`, `0.0521`), and the finding survives removal
of the top 1% residual-energy ligands. Finding: `transferable_feature_explainable_direction`.

Two facts require amendments to this document rather than celebration:

1. **The rank ceiling in section E.1 is wrong for this substrate.** Nested selection chose effective
   rank `22-26`, not `<= 8`. The `d <= 8` figure was inherited from the old Bayesian adapter's
   dimension, which section C.4 shows was never an identifiable quantity. `d` must be treated as an
   estimated property with permutation and bootstrap support, not a design constant; the Helmert
   layer of section E must be specified for the `d` that survives the direction audit, and the
   computational cost of `Sigma_{0,t}` at `d ~ 25` re-checked (a `25 x 25` SPD covariance from ~100
   cross-fitted training coefficients is at the edge of what that estimator supports, so a
   registered shrinkage or a registered `d` truncation rule is required before PD-C).
2. **The selected penalty sits at the frozen grid boundary** (`0.05 x lam_max` in all five folds).
   The grid was not extended after seeing this. Any future stage must register a wider grid *before*
   running, and must report whether the boundary is still active.

Neither fact changes the verdict. PD-M authorises nothing, and the blocking item remains a
Metz-scale independent evaluation source.
