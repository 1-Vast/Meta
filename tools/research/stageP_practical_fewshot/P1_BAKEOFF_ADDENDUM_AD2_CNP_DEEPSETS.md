# Stage P1 bake-off — arm 5 re-adjudication addendum AD2 (2026-08-19)

Amends P1_BAKEOFF_ADDENDUM_IMPL_ARMS45.md (SHA 84fea478...) BEFORE any
arm-5 run (no arm-5 artifact exists). AD1 (bank addendum) and the rest of
the arms-4/5 addendum are unchanged.

## Adjudication: deterministic Deep-Sets CNP (A), not latent NP (B)

The original arm-5 spec described a latent mu/log-var encoder with an ELBO
term, but the log-var head never participated in prediction or sampling —
decorative latent machinery. Adjudication between:

- A. Deterministic CNP / Deep Sets: per-support-item encoding,
  permutation-invariant aggregation, query decoder reads the context.
  No sampling, no ELBO.
- B. Latent NP: explicit prior/posterior, sampling, KL; strictly more
  machinery and an extra stochastic mechanism that a BASELINE arm should
  not carry (its marginal value is a separate research question).

**Decision: A.** The arm's job is to represent the "support-encoded
context" adaptation class at minimal cost with clean attribution;
determinism also keeps seed-to-seed variance down. Latent-NP language is
removed from the code and the artifact schema.

## Frozen arm-5 spec (replaces the arm-5 section of the addendum)

- Item encoder phi: Linear(2689 -> 64, ReLU) -> Linear(64 -> 64, ReLU)
  -> Linear(64 -> 64). Input per support item = concat(protein 640,
  ECFP 2048, y 1) = 2689 (dims documented in code and artifact).
- Aggregation: r = mean_i phi(item_i); EMPTY support -> r = 0 (fixed
  zero vector, no trainable prior).
- Decoder: yhat = PTrunk(xp, xl) + off_head(r); off_head =
  Linear(64 -> 1, bias=False). Consequently the k=0 context correction is
  EXACTLY zero and k=0 predictions equal the shared PTrunk output
  value-by-value (frozen test).
- Training: per-task query MSE only (no KL, no sampling). AdamW lr 3e-4,
  wd 1e-4, 6000 steps, shared sampler / monitor / checkpoint rule
  (addendum 84fea478...). Non-finite-task guard identical to arm 4.
- Eval: context encoding only (no gradient fine-tuning). Support-only by
  construction; query labels never enter.
- Parameter delta vs arm 3: 180,544 = enc1 172,160 (2689x64+bias)
  + enc2 4,160 (64x64+bias) + enc3 4,160 (64x64+bias) + off 64
  (64x1, no bias); recorded in the artifact.
- Mandatory report wording: ordinary FT / MAML use 50 support gradient
  steps; CNP uses context encoding — DIFFERENT adaptation mechanisms
  compared under the SAME data protocol (support-only, query never
  adapted on); inference procedures are not claimed identical.

## Frozen invariants (tests)

1. k=0: yhat == shared-trunk yhat (bitwise); context correction exactly 0.
2. Support permutation invariance (any support order -> same output,
   1e-6).
3. Query permutation equivariance (query row permutation permutes
   outputs identically).
4. Query-label isolation: prediction reads labels ONLY through the
   support ids (structural + behavioral test with corrupted query-label
   lookup).
5. Input dimension 640 + 2048 + 1 = 2689 matches code.
