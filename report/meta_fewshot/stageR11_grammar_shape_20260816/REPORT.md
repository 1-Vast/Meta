# Stage R11: shape-first training on the incumbent grammar trunk — falsified

Numerical authority: `COMPARE_R11_meta_val.json` (+`.rows.jsonl`), per-arm
`RESULT.json`. Gates fixed in `PREREGISTRATION.md` before any run. Three
seeds, matched 1200-step budget, zero architecture change (the exact
`similarity_only` trunk, Tanimoto transport, trained by the shape-first
method with cliff weight 1.0). Executed by `scripts/run_stage.py` (GPU
smoke first; the smoke caught two real device/shape bugs that the CPU smoke
had missed — commands and exit codes in `stage_spec.commands.jsonl`).
`meta_test` not opened.

## Result (3-seed mean)

| arm | k=0 MSE | CI | Spearman | calib | shape | k=5 MSE | k=5 cliff sign |
|---|---:|---:|---:|---:|---:|---:|---:|
| A0 (ordinary) | **2.149** | **0.580** | **0.223** | **1.236** | 0.913 | 0.915 | **0.675** |
| G1 (shape-first) | 2.405 | 0.525 | 0.073 | 1.488 | 0.917 | 0.920 | 0.650 |

G1 vs A0 paired bootstrap: k=0 -0.256 [-0.677, +0.130]; k=5 -0.005
[-0.048, +0.038].

## Gates

- **H1 fail** — k=0 MSE regresses 11.9%;
- **H2 fail** — CI 0.525 vs 0.580;
- **H3 fail** — cliff sign 0.650 vs 0.675;
- H4/H5 not evaluated (primary gates failed).

## What this falsifies

The incumbent trunk's calibration **lives in the interaction branch itself**:
routing the level term away from it (the shape-first method's core device)
strips the trunk of exactly the capacity that makes A0 the incumbent
(calibration 1.236 -> 1.488). Combined with R7/R10, the level/shape routing
trades calibration for shape **whichever architecture carries it** — on the
factorized trunk the routed level recovers only at the full budget and never
exceeds the incumbent; on the incumbent trunk the routing destroys the
calibration outright. The shape-first gains (R8/R9 cliff and shape terms)
are real, but no tested configuration converts them into a joint k=0 MSE +
CI improvement over A0.

## Status of the preregistered variable ladder (R9-R11)

1. cliff pair weight: resolved — x4 is a net negative; w=1 maximizes CI
   (0.562) and cliffs (0.606); w=2 minimizes k=0 MSE (2.119, calib 1.218).
2. shape variance term: falsified as the margin-compression cause.
3. shape parameterization (incumbent trunk + shape-first): falsified —
   routing breaks the incumbent's calibration.

The remaining recorded lever is budget scaling, which the protocol permits
only with a matched-budget A0 retrain and a preregistered learning-curve
condition. That experiment is left preregistered, not launched, pending the
next round's decision.
