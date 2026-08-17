# Stage L preregistration — support-gated assay-aware level head

Frozen before any arm trained. Date: 2026-08-18 (evening). Single seed,
development evidence; meta_test sealed and never constructed in this stage.

## Hypothesis (from two measured failures, not analogy)

Stage E and Stage J both measured the same interaction: a learned zero-shot
level head improves k=0 calibration (Stage J: level^2 1.7314 -> 1.30, the
best on record) but degrades k>=1, because at k>=1 the TRUE support labels
calibrate level near-optimally through the transport and any learned
zero-shot level estimate the head supplies replaces part of that calibration
while shrinking the transport's shape residuals (Stage J: resolved k=2/3
ranking degradation). The new hypothesis is that zero-shot level and
few-shot level are DIFFERENT mechanisms and must be trained and routed
separately: the level head should specialize in k=0 and be structurally
absent at k>=1.

## The two innovations (maximum two)

- I1l (framework): the Stage J assay-aware level head (journal/publisher
  embedding + panel set context + protein summary), gated by support size:
  its output enters the endpoint at k=0 only; at k>=1 the endpoint is the
  incumbent trunk + Tanimoto transport, byte-identical to T2.
- I2l (training): k=0-specialized level supervision — the level head receives
  gradients only from k=0 episodes; at k>=1 its output is multiplied by zero
  and its gradient is blocked, so it never competes with the transport.

## Arms

- T2 (frozen baseline, already trained and evaluated).
- L: AssayLevelModel with the support gate (both innovations).
- J (already trained, frozen): the same model ungated — the Stage J arm is
  the preregistered no-gate control; no new training.

Budget: 1,200 steps, 3 episodes/step, seed 20260815, Stage B recipe,
leak-free internal checkpoint selection, GPU verification before the arm.

## Gates (single-seed screen, frozen meta_val banks, vs frozen T2)

G1. k=0 MSE improved with a resolved paired interval.
G2. No k in {1,2,3,5} degraded with a resolved interval (the gate should
    keep L statistically indistinguishable from T2 there).
G3. Spearman/CI not degraded at any k by a resolved interval; k=0 ranking
    not degraded by more than the preregistered slack (0.02).
G4. Controls: permuted/matched-wrong above correct and wrong-protein above
    correct at every k (no inversion).
G5. Cost: parameters equal to J (<=1.15x T2); VRAM and wall time recorded.

Stop rules: S1 G1 fails; S2 G2 fails; S3 any control inverts.

Promotion path: gates pass -> 3 fixed seeds (nested k, component bootstrap)
-> freeze -> meta_test exactly once with written authorization.

meta_val figures remain single-seed development evidence.
