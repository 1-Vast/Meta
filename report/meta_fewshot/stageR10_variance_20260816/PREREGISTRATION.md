# Stage R10 preregistration: reduce the shape variance term on the C1 base

Frozen before any R10 run. Population, bank, seeds, budget, selection and the
wrong-protein contract are identical to R7/R8/R9. `meta_test` remains sealed.

## Hypothesis

The R9 audit localized the remaining CI gap to margin compression: C1's mean
absolute predicted margin is 0.097 vs A0's 0.121, and the small-gap/mid-gap
strata remain below A0 while no stratum is resolved anymore. The leading
suspect is the shape variance term (`shape_variance_weight=1.5`), which pulls
the shape branch toward the centered labels and flattens it on uncertain
pairs — variance-optimal predictions are regression-to-mean, exactly the
margin compression observed. Literature basis: ranking losses are invariant
to per-target constants and must not be diluted by a regression term on the
same branch (ActFound pairwise formulation; RankNet/LambdaLoss family).

## The single variable

C1's configuration exactly, except `shape_variance_weight 1.5 -> 0.5`:

`--no-gate --shape-variance-weight 0.5 --relative-loss-weight 1.0
 --cliff-pair-weight 1.0`, 1200 steps, 3 episodes/step, lr 6e-4 cosine,
seeds 20260815/16/17. The transport stays the retained Tanimoto+key
baseline (no query-specific gate — the gate family remains closed).

## Arms

- **A0** frozen R3R4 incumbent checkpoints;
- **C1** the R9 checkpoints (variance 1.5, cliff weight 1.0) — the
  single-variable control, already trained;
- **D1** the R10 arm (variance 0.5).

## Gates

- **G1 (primary)** D1 k=0 CI improves over C1's 0.562 and is no more than
  0.02 below A0's 0.580;
- **G2** D1 k=0 MSE does not regress beyond C1's 2.235;
- **G3** D1 k=5 activity-cliff sign accuracy stays >= 0.70 (the R8/R9 cliff
  gain is not lost);
- **G4** D1's k=0 shape term stays <= A0's 0.913;
- **G5** the improvement holds in all three seeds (direction) and the
  D1-vs-C1 CI contrast has a positive bootstrap lower bound.

## Failure conditions

If G1 fails (CI does not improve over C1), the variance term is recorded as
not the binding constraint and the next single-variable hypothesis targets
the shape parameterization itself (the anchor-mean-of-delta readout vs a
direct interaction-head readout), preregistered separately. If G1 passes but
G2/G3 fail, the tradeoff is recorded and the dose (0.5) is the formal
candidate for a three-seed re-evaluation against the standing Z gates.

## Resources

Three 1200-step runs on the RTX 4060 Laptop (~5 min each), executed by
`scripts/run_stage.py` with the smoke-first discipline; commands recorded in
the stage's `commands.jsonl`.
