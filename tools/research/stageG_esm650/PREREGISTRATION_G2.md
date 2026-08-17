# Stage G2 preregistration — multi-seed confirmation of the ESM-650M lane

Frozen before any G2 arm trained. Date: 2026-08-17 (night).

## Why this is a new stage, not a rescue

Stage G's single-seed screen met its informative content (resolved k=0
centered gain; MSE and ranking improved at every k; controls clean) but
failed gate G2 (no resolved MSE gain at k in {2,3,5} on 41 targets). The
failure mode is statistical power on a new INPUT lane, not mechanism
inertness. Per the research loop, the candidate may not be fine-tuned; a
multi-seed confirmation of the frozen recipe is the protocol's next step
for a promising lane and requires its own preregistration. This document is
it.

## Design (frozen from Stage G)

Arm G: similarity_only trunk with the local ESM-2 650M protein bank
(1280-dim pooled + 128-slot residues), Stage B loss recipe, 1,200 steps,
3 episodes/step, AdamW (3e-4 backbone 0.25x, transport 3e-4), grad clip 1.0,
amp off, float32, leak-free internal checkpoint selection (Stage B
partition, seed 20260818). Arm T2: identical with the governed 150M bank.
Seeds: 20260815, 20260816, 20260817 (the 20260815 pair already exists from
Stage D/G and is reused — its artifacts are frozen).

Evaluation: the same frozen nested meta_val banks (0,1,2,3,5), 16 queries,
2 draws, seed 73101, per seed; per-target paired differences; component
bootstrap pooling across seeds (9999 draws, seed 20260816).

## Gates (multi-seed, meta_val; meta_val stays development evidence)

G2-1. k=0 centered MSE lower for G in all 3 seeds AND the pooled interval
      resolves (hi < 0).
G2-2. MSE not degraded at any k by a pooled resolved interval; at least one
      of k=0 or k=5 improved with a pooled resolved interval.
G2-3. Spearman and CI not degraded at any k by a pooled resolved interval;
      activity-cliff sign not degraded by more than 0.03 at any k.
G2-4. Correct-support dependence preserved in all 3 seeds: permuted and
      matched-wrong MSE above correct (pooled resolved) and no inversion in
      any seed.
G2-5. Cost per seed: wall time and peak VRAM <= 1.5x T2; parameter count
      <= 1.15x T2 (correcting the mis-stated Stage G criterion).

Stop rules: any gate fails -> the lane is recorded as unconfirmed and no
meta_test is opened. All pass -> freeze architecture/hyperparameters/
checkpoints, record the frozen set, then meta_test is opened exactly once
with a written authorization.

## Davis/KIBA

If the lane is confirmed, Davis and KIBA must be trained independently from
scratch in separate experiments before any production promotion — they are
not part of this stage.
