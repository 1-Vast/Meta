# Stage P1 bake-off — arms 6/7 implementation addendum AD3 (2026-08-19)

Frozen BEFORE the first run of arms 6 (FS-CAP-style ligand-only support
encoder) and 7 (ActFound-style pairwise). Everything not specified here
is inherited verbatim from P1_BAKEOFF_PREREGISTRATION.md (SHA
59a90ef2...), P1_BACKBONE_SPEC.md (SHA 2d733684...) and the arms-4/5
addendum (SHA 84fea478...) / AD2 (SHA f8909ede...): same split, episode
bank, seeds {1,2,3} (single-seed screening first), minibatch order stream
stable_rng("stageP","porder",seed,"step",step), checkpoint rule (best
p_val draw-0 monitor), same data protocol (support-only adaptation,
query labels never enter adaptation or selection).

## Arm 6 — FS-CAP-style ligand-only support encoder

- Trunk: ligand-only variant of PTrunk — l_enc Linear(2048->64, ReLU),
  l_head Linear(64->1, bias=False), mu scalar. No protein path, no
  bilinear interaction (the interaction term needs a protein side).
  yhat = mu + l_head(l_enc(xl)). Protein-path removal is -43,138 params;
  the support encoder adds +139,584; TOTAL parameter delta vs arm 3 is
  +96,446 (recorded in the artifact).
- Support encoder (same Deep-Sets pattern as AD2): per-item
  phi_l = Linear(2049 -> 64, ReLU) -> Linear(64 -> 64, ReLU) ->
  Linear(64 -> 64); item = concat(ECFP 2048, y 1) = 2049; r = mean of
  item encodings; EMPTY support -> r = 0 (fixed, no trainable prior).
- Decoder: yhat = trunk(xl) + off_head(r); off_head = Linear(64 -> 1,
  bias=False) so the k=0 context correction is EXACTLY zero and k=0
  equals the ligand-only trunk output bitwise.
- Training: per-task query MSE only. AdamW lr 3e-4, wd 1e-4, 6000
  steps, shared sampler/monitor/checkpoint; non-finite-task guard.
- Eval: context encoding only (no gradient steps); support-only by
  construction; query labels never enter.
- Invariants (tests): k=0 zero correction + trunk equality bitwise;
  support permutation invariance (1e-6); query permutation
  equivariance (1e-6); query-label isolation; input dim 2049 matches
  code.
- Role: measures how much SUPPORT-ENCODED LIGAND-CHEMISTRY context alone
  explains query affinity (no protein features anywhere). Not
  protein-conditioned by construction.

## Arm 7 — ActFound-style pairwise within-target supervision

- Pair difference model: D(xq, xs) = h([xp_q|xl_q|xp_s|xl_s]) -
  h([xp_s|xl_s|xp_q|xl_q]) with h = Linear(5376 -> 64, ReLU) ->
  Linear(64 -> 64, ReLU) -> Linear(64 -> 1). Antisymmetric and
  identity-zero BY CONSTRUCTION (D(x,x)=0; D(s,q)=-D(q,s)).
- Training: same task sampler as arms 3-6 (identical rng consumption);
  within each sampled task, ALL ordered (i,j) pairs from the task's
  support+query cells, capped at 256 pairs per outer step (first 256 in
  fixed order); loss = MSE(D(x_i, x_j), y_i - y_j) averaged over pairs.
  AdamW lr 3e-4, wd 1e-4, 6000 steps; monitor/checkpoint rule unchanged.
- Eval = adapted level prediction: yhat(q) = mean_{s in support}
  [y_s + D(x_q, x_s)] (labels enter ONLY through the support anchors;
  query labels never enter). k=0 (no support) -> yhat = p_train label
  mean (frozen constant, computed from p_train cells only; identical to
  the k=0 ligand-only baseline prediction).
- Invariants (tests): identity-zero and antisymmetry on random inputs
  (bitwise); anchor-mean permutation invariance of the support;
  query-label isolation (behavioral); no train/val/test label is ever
  consumed by the pair module at eval except support anchors.
- Role: pairwise-difference supervision as a training mechanism; the
  eval output is an absolute-level prediction anchored to support
  labels, so its MSE/CI are directly comparable with the other arms.

## Consistency

- Monitor cadence (every 600 outer steps), monitor protocol (p_val
  draw-0, k in {0,1,2,3,5,10}), per-record metric set, and artifact
  schema follow the arms-4/5 addendum. Artifacts: P1_ARM6_FSCAP.json,
  P1_ARM7_ACTFOUND.json.
