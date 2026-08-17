# Stage Q preregistration — decoupled frozen-feature level head (single stage)

Frozen before any arm trained. Date: 2026-08-18 (night). Single seed,
development evidence; meta_test sealed and never constructed in this stage.

## Measured rationale

- Q0 (Q0_JOINT_FROZEN_IDENTIFIABILITY.json): a joint frozen-feature probe
  (journal/publisher bag + handcrafted panel statistics + frozen ESM-150M
  pooled) reaches level MSE 1.3416 vs the 2.1547 constant — the best frozen
  predictor on record, above the <=1.45 preregistered threshold.
- Stage L's failure mode was measured precisely: the level head consumed
  TRUNK-DERIVED features (protein summary, ligand encodings), so training it
  reshaped the shared trunk and degraded k>=1 ordering with resolved
  intervals. Stage Q removes that coupling by construction: the head's
  inputs are the frozen ESM bank vector, handcrafted panel statistics and a
  journal-embedding table — none of them carry trunk gradients.

## The two innovations (maximum two)

- I1q (framework): a decoupled level head over frozen assay covariates
  (journal/publisher embeddings + panel statistics + frozen ESM pooled),
  added as a scalar offset to the incumbent zero-shot endpoint, gated by
  support size (active at k=0 only).
- I2q (training): single-stage k=0-specialized supervision — the head trains
  only on k=0 episodes inside the one optimization run; the trunk's own
  protein_head absorbs the residual level through the ordinary full loss,
  so the two level mechanisms compose without any inference-time staging.

## Arms (one code path, matched budget)

- T2 (frozen baseline, already trained/evaluated).
- Q: the gated decoupled head (both innovations).
- Q-UNGATED: the same head active at every k (isolates the gate: with no
  trunk coupling, does an always-on frozen-feature head still compete with
  the transport at k>=1?).

Budget: 1,200 steps, 3 episodes/step, seed 20260815, Stage B recipe, AdamW
(3e-4 backbone 0.25x, transport 3e-4), grad clip 1.0, amp off, float32,
leak-free internal checkpoint selection (Stage B partition), GPU
verification before every arm.

## Gates (single-seed screen, frozen meta_val banks, vs frozen T2)

G1. k=0 MSE improved with a resolved paired interval.
G2. No k in {1,2,3,5} degraded with a resolved interval.
G3. Spearman/CI not degraded at any k by a resolved interval; k=0 ranking
    within the preregistered slack (0.02).
G4. Controls: permuted/matched-wrong above correct, wrong-protein above
    correct at every k (no inversion).
G5. Cost: trainable parameters <= 1.05x T2; VRAM and wall time recorded.

Stop rules: S1 G1 fails; S2 G2 fails; S3 any control inverts.

Promotion path: gates pass -> 3 fixed seeds (nested k, component
bootstrap) -> freeze -> meta_test exactly once with written authorization.

meta_val figures remain single-seed development evidence.
