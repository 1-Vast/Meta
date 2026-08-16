# Stage 1: frozen-trunk transport swap — the cross-arm contradiction is resolved

Numerical authority: `FROZEN_meta_val.json` (+`.rows.jsonl`),
`BOOTSTRAP_meta_val.json`, `FROZEN_meta_test.json` (+`.rows.jsonl`),
`BOOTSTRAP_meta_test.json`. Gates fixed in `PREREGISTRATION.md` before running.

## Hypothesis and decision

**Hypothesis.** The incumbent `grammar` trunk and the validated fixed chemical
kernel are complementary; the Stage 6 cross-arm contradiction was trunk
co-adaptation, not a mechanism conflict.

**Decision: ACCEPT.** On three frozen incumbent checkpoints, with `f0`,
shrinkage and episodes held identical and only the transport rule swapped, the
Tanimoto kernel beats the incumbent's own level baseline at k=2, 3 and 5 in
**MSE, concordance and Spearman simultaneously, with every component-level
bootstrap lower bound above zero**. No retraining was involved, so trunk
co-adaptation is excluded by construction.

**`nearest` is REJECTED** under the preregistered rule.

## Exact code change

`scripts/stage1_frozen_trunk_transport.py`. No model parameter is modified. Each
checkpoint supplies `f0` via `model.encode` for support and query ligands, plus
its own `s(n)`. Only `w_qk` varies:

| transport | rule |
|---|---|
| `mean` | `f0(q) + s(n) * mean_k r_k` |
| `tanimoto` | `f0(q) + s(n) * sum_k softmax_k(8 * Tanimoto1024) r_k` |
| `nearest` | `f0(q) + s(n) * r_argmax(Tanimoto)` |
| `incumbent` | the checkpoint's own trained transport |

## Results, `meta_val`, complete 44-episode bank, three seeds pooled

| k | mean | tanimoto | nearest | incumbent |
|---|---:|---:|---:|---:|
| 1 | 1.418 | 1.418 | 1.418 | **1.365** |
| 2 | 1.260 | 1.130 | **1.011** | 1.179 |
| 3 | 1.231 | 1.073 | **0.896** | 1.168 |
| 5 | 1.178 | 0.938 | **0.693** | 1.149 |

Concordance / Spearman at k=5: mean 0.560/0.165, tanimoto **0.671/0.438**,
nearest 0.642/0.379, incumbent 0.577/0.217.

### Paired component-level bootstrap, 9,999 draws

`tanimoto` against `mean` — the preregistered Stage 1 gate:

| k | MSE | CI | Spearman |
|---|---|---|---|
| 2 | +0.130 [+0.061, +0.208] ✓ | +0.037 [+0.021, +0.054] ✓ | +0.089 [+0.051, +0.126] ✓ |
| 3 | +0.158 [+0.080, +0.244] ✓ | +0.049 [+0.002, +0.085] ✓ | +0.143 [+0.036, +0.232] ✓ |
| 5 | +0.240 [+0.133, +0.349] ✓ | +0.110 [+0.050, +0.182] ✓ | +0.273 [+0.116, +0.444] ✓ |

**9 of 9 lower bounds above zero.** This is stronger than the Stage 6
within-checkpoint result, which had two marginal cells, and it is obtained on a
trunk that never trained with the mechanism.

`tanimoto` against the incumbent's own learned transport:

| k | MSE | CI | Spearman |
|---|---|---|---|
| 1 | -0.053 [-0.136, +0.034] | -0.007 | -0.017 |
| 2 | **+0.049 [+0.008, +0.097]** ✓ | +0.025 [-0.012, +0.061] | +0.059 [-0.019, +0.133] |
| 3 | **+0.095 [+0.032, +0.182]** ✓ | +0.038 [+0.011, +0.064] ✓ | +0.107 [+0.043, +0.169] ✓ |
| 5 | **+0.210 [+0.068, +0.379]** ✓ | +0.093 [+0.040, +0.151] ✓ | +0.221 [+0.089, +0.354] ✓ |

At k=1 the incumbent is better, as expected: the Tanimoto kernel is
mathematically degenerate at k=1 while the incumbent's gate is active there.

### `nearest` is rejected

| k | MSE (nearest − tanimoto) | CI | Spearman |
|---|---|---|---|
| 2 | +0.119 [-0.056, +0.409] | **-0.018 [-0.042, -0.001]** | -0.040 [-0.102, +0.001] |
| 3 | +0.177 [-0.109, +0.658] | +0.002 | -0.004 |
| 5 | +0.246 [-0.105, +0.865] | **-0.028 [-0.041, -0.017]** | **-0.059 [-0.086, -0.035]** |

Its MSE advantage is a point estimate with intervals that always cross zero, and
it **significantly degrades ranking** at k=2 and k=5. The preregistration
required a paired interval excluding zero before adopting it; that requirement
is not met. Hard selection is high variance — the correct reading of the Stage 6
`nearest_residual` observation is that it was an artefact of point estimates.

## `meta_test`, descriptive only

| k | mean | tanimoto | nearest | incumbent |
|---|---:|---:|---:|---:|
| 2 | 1.987 | 1.906 | 1.882 | **1.731** |
| 3 | 1.759 | 1.641 | 1.749 | **1.525** |
| 5 | 1.510 | 1.418 | 1.589 | **1.314** |

`tanimoto` against `mean` replicates partially: MSE lower bound above zero at
k=3 and k=5, Spearman at k=3 and k=5, CI at k=5; k=2 crosses zero.

`tanimoto` against the incumbent transport: MSE point estimates favour the
**incumbent** at every k, but no interval excludes zero; CI point estimates
favour `tanimoto` at every k, with a lower bound above zero only at k=2.

So on `meta_test` neither transport dominates: the incumbent has the better
squared error, the chemical kernel the better ranking, and nothing is
significant. This is reporting evidence and selected nothing.

## What is now established, and what is not

**Established.** The fixed Morgan/Tanimoto residual kernel is a genuinely
transferable k>=2 mechanism. It improves MSE *and* ranking over the level
baseline on a frozen trunk that never co-adapted to it, with 9/9 component-level
lower bounds above zero on the development split. The Stage 6 cross-arm
contradiction is explained: it was trunk and training variance, not a conflict
between trunk and mechanism.

**Not established.** That the chemical kernel should *replace* the incumbent's
learned transport. It wins on `meta_val` at k>=2 in MSE, but on `meta_test` the
incumbent has the better MSE point estimate at every k with nothing significant.
The two are separated only by ranking, where the chemical kernel is consistently
ahead.

**Not addressed.** k=0 is untouched (the transport is inactive) and k=1 is
degenerate; the incumbent remains better at k=1. k=0 is still the dominant error
term (2.0-3.7 pK^2 depending on bank).

## Interval semantics

All intervals are conditional on the three trained checkpoints: seeds are
averaged per (component, target) before components are resampled, so
seed-to-seed retraining variance is not included. Eleven `meta_val` and six to
seven `meta_test` homology components bound every interval.

## Resources

Zero training cost. Inference-only over 3 checkpoints x 4 transports x 4 support
sizes x 44 (or 42) episodes. No parameter changed, no checkpoint rewritten, no
historical result modified.

## Next decision

Stage 1's preregistered gate passed, so Stage 2 (Mac-Diff-inspired
locality-aware protein refinement, targeting k=0) is authorised. The Tanimoto
transport is held **fixed** in Stage 2 so protein refinement cannot be confounded
with transport changes.
