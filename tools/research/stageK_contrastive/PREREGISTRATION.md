# Stage K preregistration — contrastive protein-ligand coembedding

Frozen before any arm trained. Date: 2026-08-18. Single seed, development
evidence; meta_test sealed and never constructed in this stage.

## Rationale

This is the last untested framework family in the falsification ledger.
ConPLex (PNAS 2023) showed contrastive coembedding in protein language space
improves zero-shot DTI; the method_ladder family 1 (representation collapse
+ basis reallocation) was never run. The incumbent trunk's ligand
activations collapse within a target (cosine ~0.997, Stage B), and an
episodic contrastive objective is the direct training intervention against
that collapse. It operates on within/between-target ORDERING structure and
cannot create cross-component level information - the gates below reflect
that: this family is tested for its ability to move ranking and shape, not
to reach the k=0 level budget.

## The two innovations (maximum two)

- I1k (framework): a contrastive coembedding branch on the shared trunk —
  the protein summary and each ligand encoding are projected to a 128-dim
  space; the incumbent zero-shot path, transport and level branch stay
  byte-identical.
- I2k (training): episodic InfoNCE on the coembedding — within each
  episode, the (target, ligand) pairs form positives on the diagonal of the
  query grid and the other query ligands form the negatives (temperature
  0.1, weight 0.5), added to the incumbent Stage B recipe. Query labels are
  loss-only; no cross-dataset information.

## Arms

- T2 (frozen baseline): similarity_only, Stage B recipe.
- K-REG: coembedding branch trained by positive/negative regression
  (mse(zp.zl_pos, 1) + mse(zp.zl_neg, 0), weight 0.5) — framework without
  InfoNCE (I2k ablation).
- K: coembedding branch trained by episodic InfoNCE (both innovations).

Budget: 1,200 steps, 3 episodes/step, seed 20260815, AdamW (3e-4 backbone
0.25x, transport 3e-4), grad clip 1.0, amp off, float32, leak-free internal
checkpoint selection (Stage B partition), GPU verification before every arm.

## Gates (single-seed screen, frozen meta_val banks, vs frozen T2)

G1. k=0 MSE not degraded by a resolved positive interval.
G2. At least two of k in {2,3,5} improved in MSE with resolved intervals,
    OR k=0 and k=1 both improved with resolved intervals.
G3. Spearman/CI not degraded at any k by a resolved interval.
G4. Correct-support dependence preserved (permuted/matched-wrong above
    correct, resolved); wrong-protein above correct at every k.
G5. Cost: trainable parameters <= 1.05x T2; peak VRAM and wall time <= 1.5x
    T2; the within-target cosine of the coembedded ligand vectors is
    recorded as the collapse diagnostic.

Stop rules: S1 G1/G2 fail; S2 G3 fails; S3 any control inverts.

Promotion path: gates pass -> 3 fixed seeds -> freeze -> meta_test once.

meta_val figures remain single-seed development evidence.
