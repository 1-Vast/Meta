# Stage 3 (R7) preregistration: three-seed formal development validation

DRAFT — frozen only after Stage 2 (R6) screening completes and before any
Stage 3 run starts. The arm set may be reduced by the screening decision
rule; the gates below never change after freezing. If screening eliminates
the design, this document is replaced by a failure report and a new
single-variable hypothesis.

## Recorded screening history (fixed before the R7 runs)

R6a (multiplicative saturating gate) and R6b (additive correction) were each
eliminated under the R6 gates; R6c (attention-pooled level readout) failed
S1 at 300 steps (A2 k=0 2.677 vs A0 2.174) while reaching CI parity (0.550
vs 0.554). The measured mechanism facts behind the final design, all
recorded before R7: the additive correction self-cancels by construction
(`delta_hat` converges to the endpoint's implied relative; measured
nogate gaps 0.001/0.000 in R6b/R6c), so query-specific k=1 requires a
residual-scaling gate; a saturating gate cannot express the optimal
per-query scaling; zero-initialising both gate factors dead-starts the gate.
The R7 transport is therefore `t(q) = shrink * sum_k a(q,k) * rho(q,k) * r_k`
with the linear zero-mean gate `rho = 1 + u_g^T[g(e_q)h(e_k) + g(e_k)h(e_q)]`
(small random init, starts at the level shift), the bilinear antisymmetric
`delta` kept for the zero-shot shape with relative supervision. The R3R4
ladder's measured fact — the routed training's calibration advantage over
ordinary training materializes only at the full budget (B1 1.49 vs B3 1.13
at 1200 steps) — is why the decisive test runs at 1200 steps rather than the
screening budget. All 23 Stage 1 gates pass on this exact design.

## Population and data contract

- Double-cold `meta_val` (41 targets, 19 components); evaluation bank
  `evaluation_seed=73101`, `query_size=16`, one draw, all eligible targets.
- `meta_test` sealed (physical). Confirmation only in Stage 5, once, and
  only if every gate here passes.
- Wrong-protein donors: same-split most-similar cross-component,
  meta_train-fitted whitening. Training counterfactual donors: meta_train.

## Arms and budget

- **A0** incumbent: the frozen Stage R3/R4 `similarity_only` checkpoints
  (seeds 20260815/16/17, 1200 steps) — the retained comparator.
- **A2** candidate: relative-transport trunk + full training method,
  1200 steps, 3 episodes/step, lr 6e-4 cosine, seeds 20260815/16/17.
- **A1** same architecture, ordinary training (same-arch ablation for
  Innovation B), 1200 steps, 3 seeds.
- **A3/A4** key ablations (no-gate / no-counterfactual) if they survive
  screening, 3 seeds each.
- Selection: component-target-mean `full_mse_pk` averaged over k on
  `meta_val`, identical for every arm (A0 uses its retained checkpoints).

## Aggregation and inference

- equal component -> equal target -> equal draw; seeds averaged inside a
  target before components are resampled (intervals conditional on the
  trained seeds, stated as such); 9999 paired component bootstrap draws;
- per-seed direction tables; Tanimoto<0.4 tier and activity cliffs
  (Tanimoto>=0.6, gap>=1.0 pK) reported separately.

## Gates

Zero-shot (k=0):
- **Z1** A2 full MSE >= 10% below A0 (3-seed mean);
- **Z2** paired component bootstrap lower bound > 0;
- **Z3** all three seeds improve;
- **Z4** the improvement holds on the Tanimoto<0.4 tier;
- **Z5** CI and Spearman do not degrade (point estimate), and pairwise sign
  accuracy does not degrade;
- **Z6** correct protein beats matched wrong protein with lower bound > 0;
- **Z7** interaction-cut (`ligand_only`) is clearly worse than full.

k=1:
- **F1** correction is query-specific (evaluated gate, not trained);
- **F2** full beats level-only;
- **F3** correct support beats the magnitude-matched wrong label;
- **F4** with the support label and residual fixed, replacing the support
  ligand changes the prediction sensibly (population-level mean absolute
  change above a preregistered floor of 0.01 pK);
- **F5** MSE, CI and sign accuracy move in the same direction relative to k=0.

k>=2:
- **G1** full beats level-only at every k;
- **G2** full beats the fixed-Tanimoto-only baseline (`nogate` arm);
- **G3** correct support beats permuted labels and matched-wrong support;
- **G4** MSE, CI, Spearman and sign accuracy improve in the same direction;
- **G5** the gain is not explained by a better zero-shot endpoint alone
  (full-minus-zero gains exceed the A0 full-minus-zero gains at each k).

Innovation B:
- **T1** A2 beats A1 on k=0 MSE with a positive bootstrap lower bound;
- **T2** A2's improvement over A1 is not calibration-only: A2's shape_pk or
  CI improves over A1;
- **T3** the training method contributes >= 50% of the total gain:
  (A2 - A1) / (A2 - A0) >= 0.5 on k=0 MSE (point estimate);
- **T4** T1-T3 hold in the three-seed and component-bootstrap aggregates.

## Decision rule

All of Z1-Z7, F1-F5, G1-G5, T1-T4 must pass for admission. A failure is
recorded as a failure; no gate moves, no alternative metric rescues it. The
failure analysis names the mechanism cause and proposes the next
single-variable hypothesis before any further run.

## Resources

Reported per arm and seed: parameters, peak GPU memory, wall time,
throughput (episodes/s), checkpoint sha256, gradient coverage census,
gradient conflict summary, activation statistics, per-target predictions
(jsonl), and the exact command lines.
