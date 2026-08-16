# R-MAON G0 preregistration

**Frozen:** 2026-07-28, before an R-MAON runner, smoke result, or formal result existed.
**Seed:** 1729. **Environment:** `D:\anaconda\envs\drug\python.exe`.
**Parent design:** `reports/active/post_r3_multiagent_deep_route_selection_2026-07-28.md`.

**Implementation-audit amendment A1 (2026-07-28, before the first accepted formal run):** the initial
formal runner inherited `research.klbp_r3_synthetic.py::load_real_coefficients`, which called
`DualCold.panel()` only to obtain component folds. That constructor nevertheless materialized the full
Metz registry, including development affinity, before the runner selected TRAIN. Its formal artifact
was invalidated before acceptance. During the audit, a concurrent replacement attempt also ran after
the source was repaired but before this amendment made the replacement protocol explicit; it too is
invalidated and retained at `reports/active/rmaon_g0_invalidated_pre_a1.json`. A1 rebuilds the identical
fold map directly from reader-filtered TRAIN component IDs, requires the source audit to prove a
Parquet-boundary TRAIN filter, replaces Cholesky jitter with an exact PSD eigensquare-root, and records
input hashes. Seed, replicate count, multiplier count, regimes, effect size, estimators, thresholds, and
verdict rules are unchanged. Only a run executed after A1 and the amended tests may be accepted.

## 1. Scope and claim boundary

G0 asks two independent cheapest-first questions.

1. **G0-A topology:** does a label-blind prospective roster exist that can support the registered
   mechanism and predictive estimands after every firewall?
2. **G0-B estimator:** under the frozen R3 empirical noise regimes, does a direct interaction operator
   have a calibrated regular-null score and recover a known operator magnitude?

G0-B is allowed to run even if G0-A fails because it uses synthetic outcomes and can falsify the
estimation module without commissioning an assay. A G0-B pass is **module evidence only**. Real-label
training, A0 acquisition, and any predictive claim require a separate authorization; a strict
dual-cold predictive stage additionally requires G0-A to pass.

No development, Davis, confirmation, or sealed label may be read. Metz TRAIN labels may enter only
through the existing R3 construction of per-target coefficients and empirical sandwich covariances,
after the Parquet reader has projected named columns and applied `dual_cold_split == train`. The loader
must not call `DualCold.panel()`, whose constructor materializes every affinity before selecting a split.
The historical K-LBP R3 result is immutable and is not reinterpreted.

## 2. G0-A label-blind topology audit

The runner first looks for a frozen prospective manifest at:

```text
manifests/rmaon_prospective_panel.v1.json
```

The manifest, if present, must contain target family, immutable homology/profile/provenance block IDs,
preassigned ligand scaffold IDs, random inclusion probabilities, site allocation, endpoint, and planned
query status, but no affinity outcome. The gate requires all of:

* 70--155 independent multi-family blocks for a mechanism pilot and a documented expansion to at least
  423 blocks for prediction;
* at least 40 preassigned scaffold-diverse high-fidelity query ligands per target after homology,
  binding-profile, scaffold, chemical-neighbour, provenance-family, endpoint, construct, and duplicate
  firewalls;
* `PA2 >= 0.5 pK`;
* at least 80% planned power at the frozen 0.03 paired component-macro ranking gain;
* randomized, inactive-retaining sampling with known nonzero inclusion probabilities and an independent
  provenance lineage.

An absent or incomplete manifest returns:

```text
RMAON_G0_TOPOLOGY_OR_POWER_STOP
```

The existing Metz TRAIN substrate may be summarized using non-label metadata, but its 101 kinase-only
components cannot satisfy the multi-family predictive gate.

## 3. G0-B inputs and deterministic preprocessing

Inputs are fixed to:

* the valid R1 `det_proxy_card` coordinate in
  `reports/active/klbp_r1_coordinates.npz`, with its SHA-256 bound by
  `reports/active/klbp_r1.json`;
* Metz TRAIN per-target coefficient dimensions and empirical sandwich covariance matrices produced by
  the boundary-filtered `research.klbp_r3_synthetic.py::load_real_coefficients`;
* the frozen homology-component assignment, empirical component covariance-scale pool, and exact 22.5%
  S3 stress-target rule from K-LBP R3.

Let component `c` contain `n_c` targets and let there be `C` components. Target weights are:

```text
q_t = 1 / (C n_c(t)).
```

The coordinate is component-balanced centered, then whitened on its non-null eigenspace:

```text
K0 = K - sum_t q_t K_t
G  = K0^T diag(q) K0
X  = K0 U_r Lambda_r^(-1/2)
```

where eigenvalues below `1e-8 * lambda_max(G)` are removed. Thus
`X^T diag(q) X = I` up to numerical tolerance. This transform uses no synthetic outcome.

Each empirical `V_t` is symmetric positive semidefinite and may be rank deficient when a target has no
more rows than the 64-dimensional coefficient. G0-B never inverts `V_t`. Synthetic draws use the exact
eigensquare-root of `sym(V_t)`: eigenvalues in `[-1e-8 * max(lambda_max, 1), 0)` are clipped to numerical
zero, while a more negative eigenvalue stops execution. Adding Cholesky jitter or replacing `V_t` by an
isotropic covariance is prohibited.

## 4. Frozen synthetic regimes

Every regime has 200 replicates. Regime `j`, zero-indexed below, starts from
`seed = 1729 + 1000*j`; bootstrap seeds are derived deterministically from the replicate seed.

| index | regime | truth | noise / fitted coordinate |
| ---: | --- | --- | --- |
| 0 | `S1_null` | `Theta*=0` | empirical `V_t`, unmodified |
| 1 | `S1_active` | `Theta*=0.5 a c^T` | empirical `V_t`, unmodified |
| 2 | `S2_heteroscedastic` | `Theta*=0.5 a c^T` | one empirical covariance-scale draw per component |
| 3 | `S3_degenerate_signal` | `Theta*=0.5 a c^T` | S2 plus the frozen 22.5% non-positive observed-signal stress rule |
| 4 | `S5_null_coordinate` | generated with `Theta*=0.5 a c^T` and real `X` | fitted with an independently drawn, component-balanced whitened coordinate |

`a` and `c` are seeded unit-sphere draws. In each replicate, the shared coefficient is drawn coordinatewise
as `Normal(0, s_bar^2)`, where
`s_bar = sqrt(mean_t ||beta_t - mean_u beta_u||_2^2 / 64)`, exactly the R3 empirical coefficient-scale
rule. Observed coefficients are:

```text
B_t = w_bar + X_t Theta* + L_t epsilon_t,
epsilon_t ~ Normal(0, I),  L_t L_t^T = effective V_t.
```

S3 uses the existing bounded conditional redraw implementation unchanged, with the exact PSD
eigensquare-root substituted for the old jittered factor. The target mask is the frozen R3
`select_degenerate_targets(n_targets, 0.225, seed=1729+55)` mask. For S5, every replicate draws an
independent `Normal(0,1)` target matrix with the same retained dimension as `X`, then applies the identical
component-balanced centering and whitening; the active outcome continues to be generated with real `X`.

## 5. Regular-null score

Under the exact shared-global null, remove only the component-balanced intercept:

```text
R_t = B_t - sum_u q_u B_u
A_c = sum_{t in c} q_t X_t R_t^T
T_obs = || sum_c A_c ||_op.
```

Use 999 homology-component Rademacher multiplier draws:

```text
T_b = || sum_c xi_bc A_c ||_op,  xi_bc in {-1,+1}.
p = (1 + sum_b 1[T_b >= T_obs]) / 1000.
```

Reject at `p <= 0.05`. There is no optimizer and no convergence gate at the null.

The R3 coefficient object has already integrated ligand observations and does not retain a
scaffold-indexed score contribution. G0-B therefore uses the honest one-way component multiplier.
Inventing scaffold IDs for coefficient axes is prohibited. A future edge-level A0/M1 must use the
registered two-way homology-component by ligand-scaffold multiplier and cannot cite G0-B as evidence
that scaffold dependence was calibrated.

## 6. Operator recovery

Because `X` is component-balanced whitened, the unregularized direct estimate is:

```text
Theta_hat = X^T diag(q) R.
```

Let `H_:t = q_t X_t` and let `V_t^eff` be the covariance actually used for the synthetic draw. The exact
second-moment noise correction under S1/S2 is:

```text
Q_noise = sum_t trace(V_t^eff) H_:t H_:t^T
M_hat   = Theta_hat Theta_hat^T - Q_noise
gamma_hat = sqrt(max(lambda_max(sym(M_hat)), 0)).
```

The same frozen estimator is applied in S3 without a special rescue. Recovery is not evaluated for
S1-null or S5.

## 7. Frozen gates and verdicts

| gate | requirement |
| --- | --- |
| `G0B_NULL` | S1-null rejection rate `<= 0.10` |
| `G0B_WRONG_COORDINATE` | S5 rejection rate `<= 0.10` |
| `G0B_POWER_S1` | S1-active rejection rate `>= 0.80` |
| `G0B_POWER_S2` | S2 rejection rate `>= 0.80` |
| `G0B_GRACEFUL_S3` | S3 power is at least S1 power minus 0.15 |
| `G0B_RECOVERY` | in each of S1-active, S2, and S3, median `gamma_hat / 0.5` is in `[0.80, 1.25]` |

All six gates pass:

```text
RMAON_G0_NULL_SCORE_AND_RECOVERY_PASS__MODULE_ONLY
```

Any G0-B gate fails:

```text
RMAON_G0_NULL_SCORE_OR_RECOVERY_STOP
```

An overall real-stage authorization requires both G0-A and G0-B. G0-A failure therefore keeps A0,
M1, and strict dual-cold prediction blocked even if the estimator passes.

## 8. Matched controls and implementation invariants

The mandatory controls are the shared-global null, S5 independently randomized target coordinate, and
the exact `Theta=0` neural fallback. Unit tests must show:

* the component-balanced coordinate has zero weighted mean and identity weighted Gram;
* every affinity-bearing Parquet scan is projected and reader-filtered to TRAIN before materialization,
  and no `DualCold.panel()` call occurs in the G0 input path;
* PSD covariance factors reproduce their input covariances and add no jitter;
* the direct score is finite at `Theta=0` and needs no alternating direction;
* `Theta=0` is bit-identical to the shared ligand-only prediction;
* gradients reach both standard encoders and the direct operator when attached, and disappear from an
  encoder under the explicit detach control;
* the recovery correction matches a Monte Carlo second moment on a small known Gaussian problem;
* changing `--replicates` or `--seed` writes smoke artifacts and cannot certify.

The compact sequence CNN and ligand MLP are standard encoders and are not separate innovations. G0-B
freezes features and tests innovation I2, the regular-at-null direct operator. Innovation I1,
assay-monotone multi-fidelity supervision, cannot be tested without a compliant randomized,
inactive-retaining assay substrate.

## 9. Prohibited rescues and artifacts

No threshold, regime, effect size, replicate count, bootstrap count, whitening tolerance, or seed may
change after a result. No failed regime may be dropped. No `V_t` may be replaced by `sigma^2 I`. No
development, Davis, confirmation, or sealed label may be read.

```text
model/rmaon.py
research/rmaon_g0.py
tests/test_rmaon.py
reports/active/rmaon_g0.json
reports/active/rmaon_g0_decision.md
```
