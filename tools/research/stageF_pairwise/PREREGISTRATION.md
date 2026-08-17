# Stage F preregistration — pairwise learned interaction transport

Frozen before any Stage F arm trained. Date: 2026-08-17 (evening).
Development evidence only; meta_test sealed and never constructed.

## Why this family (evidence, not analogy)

- Stage L measured a pairwise signed-gap direction in embed-space
  (r +0.270 [+0.128, +0.418]) that is orthogonal to Tanimoto (correlation
  +0.026) and therefore invisible to the incumbent fixed kernel.
- The A2 moment-form operator failed because it averaged the support into one
  vector; a pairwise (query, support) operator is the smallest mechanism that
  can consume pairwise signal.
- Every learned support kernel in this project collapsed to uniform; the
  fixed Tanimoto kernel is the only weighting that ever worked. The candidate
  therefore keeps Tanimoto as an additive anchor and learns only the residual
  edge logits.

## The two innovations (maximum two)

- I1f (framework): a learned pairwise transport kernel. For each (query,
  support) pair an edge MLP consumes (query embed, support embed, support
  residual value) and emits a logit; logits = learned_scale * Tanimoto +
  edge, softmax over supports; transport = shrink(n) * sum_k w_qk * r_k.
  Label-locked: only support residual values enter; query labels never do.
- I2f (training): pairwise signed-gap supervision. Within each episode the
  predicted gaps p(q) - f0(k) are regressed against the signed label gaps
  y(q) - y(k) over the query x support grid (smooth_l1, weight 0.5). Query
  labels remain loss-only; inference consumes only support labels. This
  supervision expresses the ranking/ordering objective for transport that the
  absolute MSE cannot (a weighted average of support residuals has no
  pairwise target under pure MSE).

The zero-shot trunk is byte-identical to the incumbent similarity_only model.

## Arms

| arm | model | loss |
|---|---|---|
| T2 (frozen baseline, already evaluated) | similarity_only | Stage B recipe |
| F-ABS | PairwiseTransportModel | Stage B recipe only (framework-only ablation) |
| F | PairwiseTransportModel | recipe + 0.5 pairwise signed-gap term (both innovations) |

Budget: 1,200 steps, 3 episodes/step, seed 20260815, AdamW (3e-4 backbone
0.25x, transport 3e-4), grad clip 1.0, amp off, float32, leak-free internal
checkpoint selection (Stage B partition). GPU verification before every arm.

## Gates (single-seed screen, frozen meta_val banks, vs frozen T2)

G1. k=0 MSE does not degrade by more than +2% and its interval is not
    resolved positive.
G2. At least two of k in {2,3,5} improve MSE with resolved intervals, and no
    k degrades with a resolved interval.
G3. Ranking never trades: Spearman/CI not degraded at any k by a resolved
    interval.
G4. Correct-support dependence: permuted and matched-wrong support MSE above
    correct with resolved intervals, and not below T2's incremental
    dependence minus preregistered slack.
G5. Cost: trainable parameters <= 1.05x T2; wall time and peak VRAM <= 1.5x
    T2; gradient coverage recorded.

Stop rules: S1 G1/G2 fail; S2 G3 fails; S3 any control inverts; S4 the
learned edge logits are within 1e-2 of the Tanimoto-only baseline on average
(the framework is then declared to have collapsed to the fixed kernel).

Promotion path: gates pass -> >=3 fixed seeds -> freeze -> meta_test once.

meta_val figures remain single-seed development evidence.
