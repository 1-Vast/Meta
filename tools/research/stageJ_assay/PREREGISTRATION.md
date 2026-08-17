# Stage J preregistration — assay-aware level head + paired level alignment

Frozen before any arm trained. Date: 2026-08-18. Single seed, development
evidence; meta_test sealed and never constructed in this stage.

## Measured rationale (not analogy)

- D0c (D0c_JOURNAL_IDENTIFIABILITY.json): journal/publisher codes parsed from
  the legal panel_ids metadata carry level signal — linear probe level MSE
  1.619 vs the meta_train-constant 2.155, shuffled control 2.522, and 100%
  of meta_val episodes share a journal code with meta_train.
- D0_LEVEL_IDENTIFIABILITY: panel composition 1.887 (MLP); ESM-650M linear
  1.688; the three covariate families have non-overlapping measured shares.
- Stage E's panel head failed at k>=1 by double-fitting the level that the
  transport already carries. Stage J therefore keeps the full-prediction
  smooth_l1 term, so the level head is supervised in its RESIDUAL role by
  construction — the identified failure mode is removed, not patched.

## The two innovations (maximum two)

- I1j (framework): an assay-aware level head consuming the protein summary,
  panel composition (mean/max over query ligand encodings) and a learned
  embedding of the panel's journal/publisher codes. Journal embeddings are
  ordinary trainable parameters in the single stage.
- I2j (training): paired cross-target level alignment — within each step,
  predicted level gaps between episode pairs are regressed against true
  (transport-residual) level gaps (smooth_l1, weight 0.25), a
  BatchDTA-inspired implicit alignment that pins a consistent global scale.

## Arms (one code path, matched budget)

| arm | model | loss |
|---|---|---|
| T2 (frozen baseline) | similarity_only | Stage B recipe |
| J-NOJRNL | AssayLevelModel, journal vocab disabled | recipe + paired term |
| J-NOPAIR | AssayLevelModel | recipe only (I2j ablation) |
| J | AssayLevelModel | recipe + paired term (both innovations) |

Recipe = smooth_l1(post) + 1.0 smooth_l1(pre) + 0.5 ranknet + 0.5 centered
+ 0.05 dictionary (Stage B). Budget: 1,200 steps, 3 episodes/step, seed
20260815, AdamW (3e-4 backbone 0.25x, transport 3e-4), grad clip 1.0, amp
off, float32, leak-free internal checkpoint selection (Stage B partition),
GPU verification before every arm.

## Gates (single-seed screen, frozen meta_val banks, vs frozen T2)

G1. k=0 MSE not degraded by a resolved positive interval (expected:
    improvement from the level head).
G2. At least two of k in {2,3,5} improved in MSE with resolved intervals,
    OR k=0 and k=1 both improved with resolved intervals.
G3. Spearman/CI not degraded at any k by a resolved interval.
G4. Correct-support dependence preserved (permuted/matched-wrong above
    correct, resolved); wrong-protein above correct at every k.
G5. Cost: trainable parameters <= 1.15x T2; peak VRAM and wall time <= 1.5x
    T2; gradient coverage recorded.

Stop rules: S1 G1/G2 fail; S2 G3 fails; S3 any control inverts.

Promotion path: gates pass -> 3 fixed seeds -> freeze -> meta_test once.

meta_val figures remain single-seed development evidence.
