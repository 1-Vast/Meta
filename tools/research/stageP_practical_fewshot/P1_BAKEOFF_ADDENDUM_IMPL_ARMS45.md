# Stage P1 bake-off — arms 4/5 implementation addendum (2026-08-19)

Frozen BEFORE the first run of arms 4 (first-order MAML) and 5 (CNP).
Everything not specified here is inherited verbatim from
P1_BAKEOFF_PREREGISTRATION.md (SHA 59a90ef2...) and P1_BACKBONE_SPEC.md
(SHA 2d733684...): same split, episode bank, seeds {1,2,3} (single-seed
screening on p_val first), minibatch order stream
stable_rng("stageP","porder",seed,"step",step), checkpoint rule
(best p_val draw-0 monitor metric), and the frozen test-time protocol
(50 support steps, lr 1e-3, AdamW, best SUPPORT loss; query labels never
enter adaptation or selection).

## Arm 4 — first-order MAML

- Backbone: PTrunk (unchanged). Training episodes from the shared
  sampler (p_train targets with >= 28 unique ligands, k sampled from
  {5,10,20}, query 8; accumulate tasks until >= 256 cells, keep first
  256).
- Inner loop (training only): 5 steps, SGD, inner lr 1e-2, no weight
  decay, on support cells of each episode; outer loss = MSE on the same
  episode's query cells; outer optimizer AdamW lr 3e-4, wd 1e-4, 6000
  outer steps. First-order only: no second derivatives.
- Eval: identical to arm 3 (frozen test-time protocol; the inner loop is
  not used at eval — adaptation is the same 50-step support fit).

## Arm 5 — CNP (context-encoder style)

- Extra module (parameter delta vs arm 3 recorded in the artifact):
  support encoder g: mean over support of concat(protein_feat(640),
  ECFP(2048), y) -> MLP 2888 -> 64 (ReLU) -> latent mu (64-dim);
  decoder offset h: 64 -> 1. Prediction yhat = PTrunk(xp, xl) + h(mu).
  Diagonal log-var head exists (64-dim) but only mu is used at
  inference (MAP-style context); the log-var head is trained with the
  standard CNP ELBO at train time (beta=1.0).
- k=0: support empty -> mu = learned prior (trainable 64-dim vector);
  the PTrunk part is unchanged (frozen assertion test: k=0 predictions
  identical with the encoder module present vs removed-and-replaced-by-
  prior, when the same PTrunk weights are used).
- Eval: identical to arm 3/4 (frozen test-time protocol; the encoder
  consumes only support rows and y, never query).

## Consistency and reporting

- All three learned arms share monitor cadence (every 600 outer steps),
  monitor protocol (p_val draw-0 records, k in {0,1,2,3,5,10}, 10
  adaptation steps), and the per-record eval metric set.
- Artifacts: P1_ARM4_MAML.json, P1_ARM5_CNP.json (same schema as
  P1_ARM3_ORDINARYFT.json + arm-specific protocol block). Comparison via
  p_report.py (paired target-level bootstrap vs ligand_only/Tanimoto and
  across learned arms).
