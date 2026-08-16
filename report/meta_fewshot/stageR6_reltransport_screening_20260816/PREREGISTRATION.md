# Stage 2 (R6) preregistration: four-arm short training

Frozen before any result existed. Screening runs are for **elimination only**;
no performance claim is made from them. The formal development gates live in
Stage 3 (R7).

## Population and data contract

- Governed double-cold split
  `dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1`
  (assignment sha256 frozen by the builder).
- Development population: `meta_val` — 41 targets, 19 components, 1,411 cells.
- Evaluation bank: `evaluation_seed=73101`, `query_size=16`, one draw, every
  eligible target, nested k in {0,1,2,3,5}. Identical for every arm.
- Training population: `meta_train` only (5,643 cells / 346 targets).
- **`meta_test` is sealed physically** (`QPSMPData include_meta_test=False`);
  no arm can read it. Confirmation happens only in Stage 5 (R9), once, after
  every Stage 3 gate passes.
- Wrong-protein donors at evaluation: most similar target from a *different*
  component **within `meta_val`**, whitening mean/covariance fitted on
  `meta_train` only. Training-time counterfactual donors: `meta_train` pool.

## Arms (fixed seed, budget, episode bank)

| arm | architecture | training |
|---|---|---|
| A0 | incumbent `similarity_only` grammar | Stage R3/R4 A0 recipe (train_qpsmp), half budget |
| A1 | relative-transport trunk | **ordinary**: MSE (level+shape) + 0.5 pairwise ranking on the full prediction; no relative supervision, no routing, no counterfactuals |
| A2 | relative-transport trunk | **full method**: ranking-primary shape + relative supervision + cliff weighting + routing + counterfactuals + delta gate |
| A3 | relative-transport trunk | A2 with the delta gate disabled (`gate == 1` in the transport) |
| A4 | relative-transport trunk | A2 without the counterfactual contrasts |

Shared: 3 seeds (20260815/16/17), **300 steps** (budget amended
2026-08-16 before any arm completed; the first A0 run was killed at step 180
and no result was read), 3 episodes/step, lr 6e-4 cosine,
grad clip 1.0, AdamW, amp bf16. Validation cadence: A0 every 50 steps,
relative-transport arms every 100 steps. Checkpoint selection:
component-target-mean `full_mse_pk` averaged over k, on `meta_val`.

Two amendments recorded before any completed screening result:

1. **Identifiability pin (identify_weight=0.3).** The first screening seed
   reproduced the Stage R3/R4 predrift defect: under routing the shape
   branch's per-target mean is a null direction of every objective and the
   anchor set drifts while `target_level` compensates (measured: ligand_only
   MSE > 99 while full stayed near 3). The label-free pin
   `0.3 * mean_q(shape)^2` closes the freedom. The pre-pin seed-1 artifacts
   are invalidated and deleted; this amendment was made on the strength of
   the drift measurement, not on any arm comparison.
2. **Screening decision rule.** 300-step runs cannot resolve fine
   differences; the screening exists to eliminate clear losers. S1 kills the
   design only if A2 is **more than 10% worse** than A0 at k=0 (3-seed mean);
   S2 only if A2 is more than 15% worse than A1 at k=0. Anything closer
   advances to Stage 3, where the formal gates (Z1-Z7, F1-F5, G1-G5, T1-T4)
   decide at the full 1200-step budget.

## R6a outcome (failed screening, recorded 2026-08-16)

The first complete screening of the multiplicative-gate design **failed and
was eliminated under its own gates**: S1 (A2 k=0 2.552 vs A0 2.174, +17.4%)
and S3 (the gate inert: `nogate` gap 0.000, A3 >= A2 at every k). The
mechanism analysis named one root cause and one design change, both recorded
before the R6b runs:

* the multiplicative form `(1 + tanh(delta)) * r_k` cannot express the
  optimal per-query correction — the optimal k=1 scaling of a residual is a
  ratio, which a saturating gate mis-shapes, and with the nonlinear pair
  function the endpoint shape was invisible to the relative supervision;
* **R6b single-variable change**: the transport now implements the exact
  residual identity `r_q = r_k + (y_q - y_k) - (f0(q) - f0(k))` as
  `t(q) = shrink * sum_k a(q,k) * [r_k + delta_hat(q,k) - (f0(q) - f0(k))]`
  with the bilinear antisymmetric `delta_hat` (A3 becomes "relative
  correction off"). The endpoint shape keeps the anchor centering; the
  bilinear form makes `delta_hat` the shared, directly supervised relative
  estimate. All 23 Stage 1 gates were re-passed before R6b; every R6a
  artifact is retained for audit.

The R6b screening below reuses the identical population, seeds, budget and
gates; its decision rule is the R6a rule unchanged.

## R6b outcome (failed screening, recorded 2026-08-16)

The additive-correction design was **eliminated under the same gates**: S1
(A2 k=0 2.595 vs A0 2.174, +19.4%); the correction again measured inert at
eval time (`nogate` gap 0.001) and A2's k=1 CI (0.510) fell below its own
k=0 (0.544). Two measured facts direct the next hypothesis:

* the binding constraint is **k=0 calibration**: A2 calib 1.662 vs A0 1.289
  at 300 steps, while every new arm's shape term (~0.93) already matches
  A0's (0.885) — the routed level readout (a two-vector MLP) cannot carry
  the calibration the incumbent's whole trunk provides at this budget;
* the ranking objective *does* teach ordering: A3 reaches CI 0.575 vs A0's
  0.554 even though its endpoint spread is ~0.003 pK — a tiny but correctly
  ordered signal, which is exactly the shape signal the design was built to
  produce.

**R6c single-variable change**: the level readout only — `target_level` now
attention-pools over the residue slots with the protein summary as query
(one changed module; routing, shape, transport, and every weight unchanged).
Hypothesis: A2's calibration approaches A0's at the same 300-step budget.
All 23 Stage 1 gates re-passed before R6c; every R6b artifact is retained.

## R6c screening

Population, seeds, budget, arms and gates: identical to R6a/R6b. A0
unchanged (its three runs are reused).

## Screening gates (point estimates, 3-seed mean; elimination only)

- **S1 (primary)** A2 vs A0: eliminate only if A2's k=0 full MSE is more
  than 10% worse than A0's (3-seed mean). CI and sign accuracy reported.
- **S2** A2 vs A1: eliminate only if A2's k=0 MSE is more than 15% worse
  than A1's.
- **S3** A2 vs A3: the delta gate must not be inert — A2 must beat A3 in
  k=1 full MSE or in the mean over k in {1,2,3,5} of (A3 MSE - A2 MSE).
- **S4** A2 vs A4: mean over k in {1,2,3,5} of (A4 MSE - A2 MSE) >= 0
  (counterfactuals must not hurt).
- **S5** A2 internal controls: `full` beats `level_only` and `permuted` at
  k>=2; `full` beats `wrong_protein` at k>=2 (point estimates).

## Decision rule

A2 is eliminated only by S1/S2 above their thresholds, or by S3 showing the
gate inert. A failure of S4/S5 does not by itself disqualify the design but
records a diagnostic that directs the next single-variable hypothesis. All
formal claims wait for Stage 3. No gate is moved after the fact.

## Reporting

`COMPARE_R6_meta_val.json` (+ rows) via `stageR6_compare_arms.py`; per-arm
RESULT.json, PREDICTIONS_meta_val.jsonl, progress.jsonl, checkpoint sha256,
gradient diagnostics, peak memory, wall time. Tanimoto<0.4 tier and activity
cliffs (Tanimoto>=0.6, gap>=1.0 pK) reported separately.
