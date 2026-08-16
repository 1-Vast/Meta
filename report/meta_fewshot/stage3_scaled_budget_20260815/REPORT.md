# Stage 3 result: budget, schedule and diagnosed capacity

Seed 20260812, 2000 steps, 4 episodes/step (8,000 episodes), cosine schedule
(lr 6e-4, 5% warmup, 10% floor, `backbone_lr_scale` 1.0), validation bank
2 targets/component, frozen protocol test bank unchanged.

Numerical authority: `arm*/RESULT.json`, `WIDE_arm*.json`, `DIAG_armC.json`,
`RECHECK_armC_frozen.json`.

## Arms

| arm | architecture | trainable parameters | peak CUDA MB | best step |
|---|---|---:|---:|---:|
| A | `bpsf` (retained) | 3,788,937 | 9,110 | 1500 |
| B | `grammar` | 1,821,970 | — | 2000 |
| C | `grammar` | 7,294,171 | 6,507 | 2000 |

Arm A's peak allocation exceeds the 8,188 MB device; the retained architecture
no longer fits this budget without spilling.

## Wide bank, 42 episodes over all eligible meta-test targets (primary)

| configuration | k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---:|---:|---:|---:|---:|
| meta-train global mean | 3.441 | 3.441 | 3.441 | 3.441 | 3.441 |
| meta-train ligand prior | 3.119 | 3.119 | 3.119 | 3.119 | 3.119 |
| raw support mean | — | 2.346 | 2.180 | 1.918 | 1.523 |
| retained baseline (100 steps) | 3.589 | 2.386 | 1.993 | 1.831 | 1.581 |
| **A** `bpsf`, 2000 steps | 3.720 | 2.407 | 2.023 | 1.839 | 1.597 |
| **B** `grammar` 1.8M | 3.080 | **1.875** | **1.583** | 1.469 | 1.244 |
| **C** `grammar` 7.3M | **3.021** | 1.955 | 1.641 | **1.461** | **1.241** |
| oracle target mean (level ceiling) | 1.100 | 1.100 | 1.100 | 1.100 | 1.100 |

## Frozen protocol bank, 6 episodes (retained comparator)

| configuration | k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---:|---:|---:|---:|---:|
| retained baseline, three-seed mean | 2.12 | 1.70 | 1.34 | 1.27 | 1.25 |
| A `bpsf`, 2000 steps | 2.402 | 1.843 | 1.391 | 1.279 | 1.250 |
| B `grammar` 1.8M | 1.656 | **1.433** | **1.148** | 1.092 | 1.019 |
| C `grammar` 7.3M | **1.640** | 1.502 | 1.163 | **1.055** | **0.993** |

## What each change bought

**The budget and schedule alone bought nothing on the retained architecture.**
Arm A received the identical 2000-step cosine schedule and the identical
data-path fix, and matched its own 100-step version to within noise
(wide k=0 3.720 against 3.589; k=5 1.597 against 1.581). Its zero-shot spread
across queries *fell* to 0.0065 pK — it collapsed harder onto a constant. This
repeats the Stage 0 400-step probe and confirms the retained recipe does not
convert updates into cold-target accuracy.

**The architecture bought the improvement.** At matched seed, budget, schedule
and banks, arm C beats arm A by 0.70 / 0.45 / 0.38 / 0.38 / 0.36 MSE at
k=0/1/2/3/5 on the wide bank.

**Capacity 1.8M -> 7.3M bought no consistent MSE gain.** B wins k=1 and k=2, C
wins k=0, k=3 and k=5, and the differences are 0.02-0.08. Its one reproducible
benefit is control sign: arm B's magnitude-matched permutation control is
**inverted** at k=2, 3 and 5 on the wide bank (permuted 1.563/1.417/1.229
against correct 1.583/1.469/1.244), whereas arm C's is correctly signed at
every k. Reported as a negative result for width.

**The data-path fix bought the wall clock.** Episode materialization fell from
1,158 ms to 15.2 ms per episode; the training step fell from 5.57 s to 0.75 s
(`bpsf`) and 0.205 s (`grammar`, 1.8M). That is what made 8,000-episode runs
affordable at all.

## Controls, arm C, wide bank

| k | full | sar_cut (adaptation cut) | permuted (matched) | matched wrong support | wrong-protein gap |
|---|---:|---:|---:|---:|---:|
| 1 | 1.955 | 2.048 | 4.936 | 3.871 | +0.246 |
| 2 | 1.641 | 1.776 | 1.670 | 4.048 | +0.049 |
| 3 | 1.461 | 1.612 | 1.473 | 5.276 | -0.011 |
| 5 | 1.241 | 1.361 | 1.251 | 5.115 | -0.010 |

* The query-specific adaptation channel contributes 0.093 / 0.135 / 0.151 /
  0.119 MSE at k=1/2/3/5. **At k=1 this channel does not exist in the retained
  model**, where `sar_cut` equals `full` identically.
* Correct support beats magnitude-matched permuted support at every k, but by
  only 0.010-0.029 at k>=2. The adaptation is query-specific; the
  *support-identity*-specific component is small.
* Wrong protein costs 0.246 MSE at k=1 and 0.049 at k=2; at k=3 and k=5 level
  calibration dominates and the gap is within noise.

## Zero-shot degeneracy: improved, not solved

| quantity | retained baseline | arm A (2000 steps) | arm C |
|---|---:|---:|---:|
| zero-shot spread across queries in an episode | 0.065 pK | 0.0065 pK | **0.186 pK** (wide) / 0.196 pK (48 random draws) |
| zero-shot change on cross-component protein swap | 0.0093 pK | — | **0.438 pK** |
| wrong-protein zero-shot MSE gap | +0.012 | +0.005 | **+0.674** |

Carried-forward **R2 passes** (gap 0.674 against a 0.05 threshold): the
zero-shot endpoint is now genuinely protein-conditioned. Carried-forward
**R1-secondary does not pass**: 0.186-0.196 pK against a 0.20 pK threshold and
a 0.93 pK label spread. The endpoint still explains only a small share of
within-target ligand variation.

## Gradient coverage, arm C

| k | 0 | 1 | 2 | 3 | 5 |
|---|---:|---:|---:|---:|---:|
| zero-gradient trainable tensors | 7 | 2 | **0** | **0** | **0** |
| retained baseline, same measurement | 23 | 22 | 17 | 17 | 17 |

The k=0 zeros are the seven `transport.*` tensors, which have nothing to act on
without support; the k=1 zeros are `transport.key.weight` and
`transport.log_temperature`, because a softmax over one support is 1 regardless
of the key. Both are semantic, not unreachable code. The query-specific gate
**is** live at k=1, which is the defect this stage set out to fix.

## Reproduction

`RECHECK_armC_frozen.json` reloads `armC_grammar_7M/checkpoint.pt` and
reproduces 1.639625 / 1.502116 / 1.163350 / 1.055060 / 0.993290 exactly.

## Gate outcome

All five Stage 3 gates pass for arm C. Budget was **not** saturated: arm C's
best validation step was the last step, 2000.
