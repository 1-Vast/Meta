# Stage I preregistration — live ESM-2 150M LoRA lane (LM conditioning)

Frozen before any arm trained. Date: 2026-08-17 (night). Single seed,
development evidence; meta_test sealed and never constructed in this stage.

## Rationale (measured)

D0_LEVEL_ANATOMY measured the protein-sequence share of level variance at
+11.9% with FROZEN embeddings; the frozen ESM-150M/650M input lanes (Stage
G/G2) were not confirmed. The remaining sequence-side degree of freedom is
adapting the LM itself. This stage tests whether end-to-end (adapter)
conditioning of the protein LM by the single DTA objective extracts more
level/shape signal than any frozen bank.

## The two innovations (maximum two)

- I1 (framework): the frozen 150M protein bank is replaced by a live ESM-2
  150M encoder (chunked <=1022 residues, mean/slot pooled, identical pooling
  policy), carrying LoRA adapters (rank 8, alpha 16) on every attention
  projection. The DTA trunk (similarity_only) is byte-identical.
- I2 (training): single-stage end-to-end LM conditioning - the adapters are
  trained by the same episodic DTA loss in the one optimization run; no
  pretrain/finetune staging, no closed forms.

External data disclosure: facebook/esm2_t30_150M_UR50D, revision
a695f6045e2e32885fa60af20c13cb35398ce30c, local snapshot.

## Arms

- I: live encoder, LoRA trainable (both innovations).
- I-FROZEN: identical live path, adapters frozen, encode under no_grad
  (isolates the trainability of the LM lane; also measures the live-encoder
  cost versus the frozen-bank T2 reference).
- Baseline: frozen T2 (Stage D/G2, seed 20260815).

Budget: 1,200 steps, 3 episodes/step, seed 20260815, AdamW (3e-4 backbone
0.25x, transport 3e-4, LoRA 3e-4), grad clip 1.0, amp off for the trunk,
bf16 autocast for the LM encode, leak-free internal checkpoint selection.
Declared bank change: the internal bank uses up to 1 target per component
and 1 draw (live-encoder cost); the fit/internal component partition is
unchanged. Declared gradient bound: during training, LoRA gradients flow
through the first 1022-residue chunk only (later chunks are feature-only,
encoded without gradient) so long proteins cannot exhaust GPU memory; the
evaluation encoder always sees the full sequence. GPU verification before
every arm.

## Gates (single-seed screen, frozen meta_val banks, vs frozen T2)

G1. k=0 MSE not degraded by more than +2% and its interval not resolved
    positive.
G2. At least two of k in {2,3,5} improved in MSE with resolved intervals.
G3. Spearman/CI not degraded at any k by a resolved interval.
G4. Correct-support dependence preserved (permuted/matched-wrong above
    correct, resolved).
G5. Cost: trainable parameters <= 2.0x T2 (LoRA adds 1.23M), peak VRAM
    <= 1.5x T2, wall time recorded.

Stop rules: S1 G1/G2 fail; S2 G3 fails; S3 any control inverts.

Promotion path: gates pass -> 3 fixed seeds -> freeze -> meta_test once.

meta_val figures remain single-seed development evidence.
