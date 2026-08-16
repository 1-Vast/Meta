# AnchorDelta P0 preregistration

## Objective

Test whether the existing frozen protein-ligand representation contains a
transferable relative-affinity signal under a fresh homology-component holdout
carved from TRAIN. The experiment is a kill test, not a final benchmark.

## Model boundary

Only `model/anchordelta.AnchorDelta.head` is trainable. The protein encoder,
ligand projection, interaction feature, and train-only `B0` are frozen. The
relative operator is structurally defined as

`Delta(p, q, i) = (h(p, q, i) - h(p, i, q)) / 2`.

The implementation therefore guarantees `Delta(q, i) = -Delta(i, q)` and
`Delta(i, i) = 0` without an auxiliary symmetry penalty.

## Training protocol

- Endpoint: pKi.
- Split: deterministic 20% homology-component holdout from TRAIN; the frozen
  strict development roster is not used for model selection.
- Pair labels: within-target `y_j - y_i`, sampled with a fixed cap per target.
- Orientation: random pair order during training.
- Objective: Huber difference loss for P0. Order loss remains available in the
  module but is not used to make a first-pass gate decision.
- Aggregation: uniform mean of `y_i + Delta(q, i)` across the five support
  anchors. No label-dependent weighting or gating is allowed.

## Required controls

The same episode/query rows must be used for all arms:

- calibration baseline;
- wrong-target support labels (the uniform mean is intentionally invariant to a
  pure column permutation, which is retained only as a label-permutation
  invariance check);
- wrong protein;
- protein-free comparator and parameter-matched absolute head in the full run;
- pair-order swap and diagonal tests;
- similarity-bin and same-assay versus cross-assay reporting.

Bootstrap units are homology components, never individual pairs or rows.

## PASS gate

All of the following must hold on the held-out components:

1. AnchorDelta versus calibration, ligand-only, and gradient baselines has a
   positive 95% component-bootstrap lower bound for RMSE, within-target
   Spearman, and pairwise accuracy.
2. Correct support labels outperform deterministic permutation controls.
3. Correct protein outperforms wrong-protein and protein-free controls.
4. The gain is present outside the highest ligand-similarity bin.
5. A second seed has the same sign and the effect exceeds seed-to-seed noise.

Any failed condition is `P0_STOP`; do not add graph encoders, context
transformers, Bayesian uncertainty heads, or task schedulers after a stop.

## Current smoke status

`scripts.anchor_delta` ran with two fit targets, two gate targets, two pairs per
target, and one epoch. The chain completed and wrote
`reports/active/anchordelta_p0_smoke.json`, but it contains only one independent
homology component and is explicitly non-decisional.
