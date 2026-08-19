# Stage P1 bake-off backbone spec (2026-08-19, frozen)

Frozen BEFORE any learned-arm run, per P1_BAKEOFF_PREREGISTRATION.md §3.
All learned arms (3-8) share this trunk unless the arm's frozen addendum
mandates a delta (recorded as parameter-count delta + ablation).

## Inputs (label-free, per cell)

- protein: 640-dim ESM2-t30-150M global-slot pooled feature from
  dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank
  (frozen external representation; ablated separately as an explicit arm
  variant protein_feature=none in the final report; no ESM fine-tuning).
- ligand: ECFP4, radius 2, 2048 bits (RDKit 2023.09.6, frozen) from
  corpus ligands.jsonl SMILES.
- No target ID, no panel ID, no scaffold string, no assay covariates in
  the trunk (assay context enters only the assay-aware arm variant, if
  any, as an explicit ablated variable).

## Trunk (shared across learned arms)

- p_enc: Linear(640 -> 64, ReLU)
- l_enc: Linear(2048 -> 64, ReLU)
- p_head: Linear(64 -> 1, bias=False)
- l_head: Linear(64 -> 1, bias=False)
- mu: scalar bias
- interaction: low-rank bilinear, rank R=16:
  A: Linear(64 -> 16, bias=False); B: Linear(64 -> 16, bias=False);
  inter_scale scalar (init 1.0); inter_bias scalar (init 0.0);
  yhat = mu + p_head(p_enc(xp)) + l_head(l_enc(xl))
       + inter_scale * sum(A(p_enc(xp)) * B(l_enc(xl)), -1) + inter_bias
- Param count and the exact init scheme are pinned in the implementation
  artifact; init = xavier_uniform everywhere, bias zero, inter_scale 1.

## Training (episodic, on p_train targets only)

- Optimizer: AdamW, lr 3e-4, wd 1e-4; steps 6,000; batch 256 (cells from
  the p_train split, keyed minibatch order shared across arms:
  stageP|order|{seed}); no gradient clipping.
- Episode loss (training phase): for each p_train target drawn by the
  shared episode sampler, support = first k of its ligand-unique order
  (k sampled uniformly from {5,10,20} per episode), query = next 8;
  MSE on support + query (ordinary arms train on both; MAML/CNP/
  FS-CAP arms modify this per their addenda).
- Level/shape ownership: NOT trained here (that is the Phase-6
  innovation, ablated separately with the trunk frozen); the trunk's
  additive heads and interaction head are trained jointly with plain
  MSE.

## Test-time protocol (P1, frozen)

- For each bank record: adapt on support only; checkpoint = best SUPPORT
  loss among the adaptation steps (query labels never enter adaptation
  or selection — asserted by tests); predict query.
- Ordinary fine-tuning arm: 50 support steps lr 1e-3 (same optimizer
  class), full trunk updated (delta-parameter count recorded).
- k=0 path: no adaptation; the pretrained trunk predicts directly; the
  few-shot modules of the other arms are bypassed with zero support
  (frozen assertion test: k=0 predictions identical with the few-shot
  module present or removed).

## Arm-specific deltas (each in its own SHA-frozen addendum before first
## run; AdaMBind addendum requires FULL-text inspection first)

- 3 ordinary fine-tuning: trunk as above.
- 4 first-order MAML: inner loop lr 0.01, 1 inner step on support, outer
  on query during training; test-time = same support adaptation as (3).
- 5 CNP: support encoder = shared MLP over (xp, xl, y) -> mean 64-dim
  context, concatenated to the trunk's p_enc/l_enc outputs; decoder
  unchanged trunk; adaptation happens only through the context.
- 6 FS-CAP-style: ligand-only support encoder (xp dropped from the
  trunk; p path removed -> parameter delta recorded).
- 7 ActFound-style: pairwise within-target supervision: predict
  yhat(q) - yhat(s) from (xp_q, xl_q, xp_s, xl_s); identity-zero and
  exchange-antisymmetric; eval = adapted level prediction.
- 8 AdaMBind-style: per published paper (addendum required before run).
