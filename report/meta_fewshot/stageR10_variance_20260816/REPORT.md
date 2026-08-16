# Stage R10: shape-variance reduction — hypothesis falsified

Numerical authority: `COMPARE_R10_meta_val.json` (+`.rows.jsonl`), per-arm
`RESULT.json`. Gates fixed in `PREREGISTRATION.md` before any run. Three
seeds, matched 1200-step budget, C1 base (cliff weight 1.0), one changed
variable (`shape_variance_weight 1.5 -> 0.5`). Executed by
`scripts/run_stage.py` (smoke first; commands in `stage_spec.commands.jsonl`).
`meta_test` not opened.

## Result (3-seed mean)

| arm | k=0 MSE | CI | Spearman | calib | shape | k=5 MSE | k=5 cliff sign |
|---|---:|---:|---:|---:|---:|---:|---:|
| A0 | 2.149 | **0.580** | **0.223** | **1.236** | 0.913 | 0.915 | 0.675 |
| C1 (var 1.5) | 2.235 | **0.562** | 0.159 | 1.332 | 0.903 | 1.096 | **0.782** |
| D1 (var 0.5) | 2.285 | 0.552 | 0.148 | 1.358 | 0.927 | 1.111 | 0.694 |

## Gates

- **G1 fail** — D1 CI 0.552 does not improve over C1's 0.562 (seed 1's
  0.605 was a seed outlier);
- **G2 fail** — k=0 MSE 2.285 regresses past C1's 2.235;
- **G3 fail** — k=5 cliff sign 0.694 falls below 0.70 (the R8/R9 cliff gain
  was lost);
- **G4 fail** — shape 0.927 regresses past A0's 0.913;
- G5 not evaluated.

## Verdict

**The variance-term hypothesis is falsified.** The shape variance term is
not the cause of the margin compression: halving it degrades CI, MSE, shape
and the cliff ordering together. It is recorded as a clean negative under
the preregistered failure condition, which directs the next single-variable
hypothesis to the shape parameterization itself (the anchor-mean-of-delta
readout vs a direct interaction-head readout) — implemented as R11, the
incumbent grammar trunk trained with the shape-first method
(`scripts/train_grammar_shape.py`, zero architecture change), which tests
the parameterization question on the architecture with the best-known
calibration.
