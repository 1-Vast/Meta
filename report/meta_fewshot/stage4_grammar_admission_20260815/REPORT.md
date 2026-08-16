# Stage 4 result: three-seed governed admission run

Configuration exactly as preregistered (`InteractionGrammarModel`, 7,294,171
trainable parameters, 2000 steps, 4 episodes/step, cosine schedule, lr 6e-4,
seeds 20260812 / 20260813 / 20260814). Each seed saw 8,000 episodes, all 399
meta-train targets, and ~10,700 of 12,633 meta-train query cells. Peak CUDA
memory 6,507-6,863 MB on an 8,188 MB device. Best validation steps: 2000, 1750,
750.

Numerical authority: `nested/RESULT.json`, `nested/PREDICTIONS.jsonl`,
`ADMISSION_TABLE.md`, `FROZEN_seed*.json`, `WIDE_seed*.json`.

## Verdict

**The governed admission gate does not pass.** Gates 4, 5 (permutation half) and
6 fail. The improvement over the retained baseline is large, reproducible across
three seeds, and holds on every bank — but the evidence attributes it to the
**zero-shot interaction trunk plus target-level calibration**, not to
query-specific SAR transfer.

## Absolute performance

### Frozen protocol bank, `evaluation_seed=73101`, three-seed mean

| configuration | k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---:|---:|---:|---:|---:|
| retained baseline | 2.115 | 1.704 | 1.341 | 1.270 | 1.253 |
| **this stage** | **1.786** | **1.493** | **1.155** | **1.077** | **1.026** |
| reduction | -0.329 | -0.211 | -0.186 | -0.193 | -0.227 |
| reduction | -15.6% | -12.4% | -13.9% | -15.2% | -18.1% |

Per seed (frozen bank):

| seed | k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---:|---:|---:|---:|---:|
| 20260812 | 1.640 | 1.502 | 1.163 | 1.055 | 0.993 |
| 20260813 | 1.815 | 1.445 | 1.092 | 1.033 | 0.983 |
| 20260814 | 1.902 | 1.531 | 1.209 | 1.143 | 1.100 |

Every seed is below the **best** retained-baseline seed at every k. The retained
baseline's per-seed range was 2.047-2.150 (k=0) and 1.217-1.283 (k=5).

### Wide bank, all 42 eligible meta-test targets, three-seed mean

| configuration | k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---:|---:|---:|---:|---:|
| meta-train global mean | 3.441 | 3.441 | 3.441 | 3.441 | 3.441 |
| meta-train ligand prior | 3.119 | 3.119 | 3.119 | 3.119 | 3.119 |
| raw support mean | — | 2.346 | 2.180 | 1.918 | 1.523 |
| retained baseline (seed 20260812 only) | 3.589 | 2.386 | 1.993 | 1.831 | 1.581 |
| **this stage** | **3.246** | **1.977** | **1.642** | **1.493** | **1.271** |
| oracle target mean (level ceiling) | 1.100 | 1.100 | 1.100 | 1.100 | 1.100 |

The retained baseline is available on the wide bank for one seed only; that
comparison is therefore point-to-mean, not seed-matched.

## Governed nested manifest, 42 targets, 6 components, three seeds pooled

| k | arm | MSE | CI | Spearman |
|---|---|---:|---:|---:|
| 0 | full = zero-shot | 2.847 | 0.647 | 0.372 |
| 1 | full | **1.882** | 0.610 | 0.257 |
| 1 | level / adaptation-cut | 1.993 | 0.647 | 0.372 |
| 1 | permuted (magnitude matched) | 4.846 | 0.672 | 0.416 |
| 1 | matched wrong support | 3.417 | 0.678 | 0.414 |
| 1 | wrong protein | 2.044 | 0.602 | 0.239 |
| 2 | full | **1.694** | 0.603 | 0.255 |
| 2 | level / adaptation-cut | 1.746 | 0.647 | 0.372 |
| 2 | permuted | 1.632 | 0.637 | 0.347 |
| 3 | full | **1.491** | 0.571 | 0.169 |
| 3 | level / adaptation-cut | 1.513 | 0.647 | 0.372 |
| 3 | permuted | 1.416 | 0.622 | 0.307 |
| 5 | full | 1.329 | 0.608 | 0.257 |
| 5 | level / adaptation-cut | **1.318** | 0.647 | 0.372 |
| 5 | permuted | 1.328 | 0.610 | 0.302 |

## Gate outcome

| # | requirement | outcome |
|---|---|---|
| 1 | k=0 does not materially regress, preferably improves | **pass** (2.115 -> 1.786 frozen; 3.589 -> 3.246 wide) |
| 2 | k=1 gain is query-specific (`full` beats `sar_cut`) | point estimate passes (1.882 against 1.993); no positive lower bound |
| 3 | k=2,3,5 improve against the retained baseline in every seed | **pass** (frozen bank, all three seeds) |
| 4 | `full` beats `sar_cut` with a positive bootstrap lower bound | **FAIL** at every k: +0.110 [-0.069, +0.413], +0.052 [-0.182, +0.375], +0.022 [-0.156, +0.286], -0.011 [-0.115, +0.126] |
| 5 | correct support beats permuted and matched wrong support | matched wrong support **pass** with positive lower bounds at every k; permutation **FAIL** at k=2,3,5 (permuted is *better* by 0.062 / 0.075 / 0.001) |
| 6 | CI and Spearman improve together with MSE | **FAIL**: CI falls from 0.647 to 0.571-0.610 and Spearman from 0.372 to 0.169-0.257 at every k>0 |
| 7 | checkpoint re-evaluation reproduces | **pass**, exactly (`RECHECK_armC_frozen.json`) |
| 8 | no dead trainable branch | **pass**: 0 zero-gradient tensors at k>=2; the k=0 and k=1 zeros are semantic (`transport.*` has nothing to select from) |

## What the failures mean

1. **The level channel does all the statistically supported work.** Adding a
   constant per target cannot change within-target ranking, which is why
   `level_only` and `zero_shot` have identical CI and Spearman. The full model's
   MSE advantage over `level_only` is 0.11 / 0.05 / 0.02 / -0.01 at k=1/2/3/5
   and none of these survives a paired component bootstrap over six components.
2. **The transferability gate degrades ranking.** `full` is the only arm whose
   CI and Spearman fall below the zero-shot values, at every k. The
   query-specific coefficient is buying a small squared-error reduction by
   shrinking predictions toward the target level, at the cost of within-target
   discrimination.
3. **Support identity is not being used.** Permuting the support labels leaves
   `mean(r)` exactly unchanged, so the permutation contrast isolates the
   query-specific channel. It is negative at k=2, 3 and 5. Only the k=1
   magnitude-matched flip (which does change the level) shows a large gap.
4. **Matched wrong support is rejected strongly** (+1.5 to +3.2 MSE with
   positive lower bounds), but that control is dominated by the donor target's
   different absolute level, so it is not evidence of binding specificity.

## What is established

* A protein-conditioned zero-shot endpoint now exists. Swapping in a
  cross-component donor protein moves the zero-shot output by 0.438 pK
  (retained baseline: 0.0093 pK) and costs 0.21-0.67 MSE (retained: 0.012).
* k=0 improved by 15.6% on the frozen protocol bank and now beats the
  meta-train global-mean and ligand-prior references on the wide bank, which the
  retained baseline did not.
* Every k improved by 12-18% on the frozen bank, in all three seeds.
* The zero-shot endpoint still explains little within-target ligand variation:
  spread across queries is 0.087-0.186 pK against a 0.93 pK label spread, so the
  carried-forward R1-secondary gate (> 0.20 pK) still does not pass.

## Standing caveats

* The meta-test split has been consumed by prior architecture search
  (`test_used_for_tuning: true`); there is no untouched confirmation cohort.
* Six independent homology components bound every bootstrap.
* No SOTA, biological-mechanism, or confirmatory Cold Target claim is
  authorized by this stage.
