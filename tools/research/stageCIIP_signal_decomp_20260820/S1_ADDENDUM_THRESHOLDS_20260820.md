# S1 ADDENDUM — ADJUDICATION THRESHOLDS (frozen 2026-08-20, before any S1 fitting)

Replaces the PREREGISTRATION.md B.8 placeholders. Derived ONLY from the S0
power table (S0_AUDIT.json item 8) and S0 diagnostics (item 7). No S1 arm
result has been computed or observed at freeze time.

## Frozen numbers (from S0)
- sigma_R2 = 0.599, sigma_Spearman = 0.218 (per-pair metric dispersion of a
  legal LOPO profile predictor, train+val 40 pairs).
- MDE(80%): delta-R2 = 0.566; delta-Spearman = 0.208 (parent-cluster
  bootstrap lo2.5 > 0, 2000 draws, clusters [2,2,1,1,2,1]).

## Adjudication rules (frozen)

For each proposition, on its PRIMARY contrast (prereg B.8 list):

- S1-PASS-X := (bootstrap lo2.5 > 0 on the paired contrast)
  AND (LOPO sign stability: the aggregate contrast sign is unchanged when
  excluding any single test parent)
  AND (point estimate >= tau_min, where tau_min is the minimum-effect bar
  below which a claimed effect is practically empty):
    tau_min(delta centered R2) = 0.05   [all T1/T2/T0 contrasts]
    tau_min(delta Spearman)    = 0.10   [T0m contrast]
- If lo2.5 > 0 and LOPO stable but the point estimate < MDE(80%), the PASS
  is labeled "PASS (power-labeled: below MDE at n=9)". This label is
  informational and does not weaken the directional rule.
- S1-NULL-ALL := every trained arm cell is at or below its F7f/F8 floor cell
  + 0.03 R2 (tolerance) on every estimand, AND the within-pair ligand-label
  permutation control shows no systematic loss when destroying labels
  (permutation CI crossing 0), AND no contrast satisfies PASS.
- S1-UNRESOLVED-* := CI crosses 0 (lo2.5 <= 0) while NULL-ALL conditions are
  not all met. Legal terminal state; no further model additions permitted to
  rescue it (prereg B.9 rule 4 and plan 9.3).
- Secondary contrasts (T0 centered R2 for A; residual R2; F2w variants) are
  reported, never gating.

## Pre-registered power expectations (recorded at freeze, non-binding labels)

- The largest legally reachable ceiling gap on T1 (F9 sibling-profile vs F7f
  ligand pattern) is expected in the 0.10-0.20 R2 range (S0: LOPO per-pair
  median 0.326, ligand-global 0.060), i.e. BELOW MDE(80%) = 0.566. The
  likely terminal grade for most contrasts is therefore UNRESOLVED (power)
  unless an arm concentrates an unexpectedly large effect.
- This expectation must not be used to weaken or strengthen any rule after
  results are seen.

## Bootstraps
- All contrasts: parent-cluster bootstrap, 2000 draws, keyed
  S1.boot.<contrast>.<seed>; 2.5/97.5 percentiles; per-parent deltas from
  per-pair metric differences; bootstrap mean never a point estimate.
