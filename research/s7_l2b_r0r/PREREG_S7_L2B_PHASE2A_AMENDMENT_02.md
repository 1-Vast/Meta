# Phase 2A preregistration — computational amendment 02

Parent: `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md`
(SHA-256 `4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e`).
Prior amendment: `PREREG_S7_L2B_PHASE2A_AMENDMENT_01.md`.

Written 2026-08-10, **after** Phase 2 (teacher conditionality) completed and
**before** any Phase 3/4 decomposition, coupling or rewiring metric was
computed. Neither the parent nor amendment 01 is edited, so both hashes remain
valid.

Both items are execution-budget and numerical declarations. Neither is
result-dependent.

## B1 — closed-form evaluation of the tied-AP expectation

Amendment A1 fixed the estimator: the exact expectation of average precision
under a uniformly random ordering inside each tied block. A1 evaluated it with a
Python loop over blocks. Phases 3 and 4 evaluate the *same* quantity in closed
form, so the rewiring null over 52,062,975 cells x 20 rewires x 3 arms is
tractable.

Using the cumulative harmonic table `H(x) = sum_{u=1..x} 1/u`, the per-block
term becomes

```text
S1  = sum_{j=1..n} 1/(a+j) = H(a+n) - H(a)
sum_{j=1..n} (j-1)/(a+j)   = n - (a+1) * S1
block contribution = (k/n) * [ (b+1) * S1 + (k-1)/(n-1) * ( n - (a+1) * S1 ) ]
```

with the `(k-1)/(n-1)` factor set to zero when `n = 1`. This is algebraically
identical to A1, not an approximation. It is validated in-run on two fronts,
both reported:

1. with no ties it must reproduce ordinary average precision exactly
   (tolerance `1e-9`);
2. with a binary tied score it must reproduce the Monte-Carlo expectation to
   within Monte-Carlo error (tolerance `0.03` at 400 replicates).

Failure of either check is a fail-closed estimator-contract error.

## B2 — ridge-IRLS budget for the label-fitted Rasch oracle

The parent registers the teacher additive null
`logit P(Y_ra = 1) = mu + alpha_r + beta_a` with ridge `1e-6`. On a matrix that
is 0.07% positive, most residues and most atoms carry no positive at all, so
their coefficients are separated and diverge. The fit is therefore run with a
fixed, declared budget rather than to convergence:

- at most **6** IRLS/Newton steps, each solving a weighted additive least
  squares by at most **12** ALS iterations;
- the linear predictor is clipped to `[-30, 30]`;
- the achieved final step size `max |delta eta|` is reported for every complex
  (maximum and median), so the degree of non-convergence is visible rather than
  hidden.

This object is an **ORACLE ceiling only**. It is label-fitted, it is never a
deployable model arm, and it never enters a Gate. Its budget therefore affects
only the tightness of a reported ceiling, never a verdict.
