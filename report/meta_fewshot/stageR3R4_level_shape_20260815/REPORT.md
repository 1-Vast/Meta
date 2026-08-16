# Stage R3/R4: level-shape factorization and gradient-routed training

Numerical authority: `COMPARE_v2_3seed.json` (+`.rows.jsonl`),
`LADDER_v2_3seed.json`, and the per-arm `RESULT.json` files.
`COMPARE_seed20260815.json` retains the superseded low-resolution ladder. Gates
fixed in `PREREGISTRATION.md` before any result existed.

Population: governed double-cold `meta_val` — 41 targets, 19 components, **zero**
exact-ligand / scaffold / component / document overlap with the training block,
81.6% of ligands below Tanimoto 0.40. Three seeds (20260815/16/17), matched
1,200-step budget, identical evaluation bank. Seeds are averaged inside a target
before components are resampled, so **every interval is conditional on those
three trained seeds**. **`meta_test` was not opened.**

## Verdict

| claim | verdict |
|---|---|
| the level and shape objectives actively conflict inside one trunk | **supported — measured directly** |
| the training method beats conventional training of the same architecture | **supported — resolved at k=0 and on the low-similarity tier** |
| the factorized architecture beats the incumbent trunk | **rejected — significantly worse** |
| routing alone, or counterfactual alone, is separately resolved | **not supported** |
| the method produces protein specificity the incumbent lacks | **supported at k >= 2** |
| removing the level/shape conflict unlocks within-target ordering | **falsified** |
| the zero-shot admission target against the incumbent | **not met** |

## 1. The measurement that motivates Innovation B

Per-objective, per-module gradients at step 1,200, arm B1 (factorized
architecture, ordinary shared-gradient training):

| module | level-objective norm | shape-objective norm | **cosine** |
|---|---:|---:|---:|
| `interaction` | 6.76e-2 | 2.20e-2 | **-0.334** |
| `ligand_encoder` | 2.35e-1 | 1.55e-1 | **-0.532** |
| `channels` | 1.32e-1 | 3.86e-2 | -0.146 |
| `level_head` | 9.53e-1 | 4.73e-8 | -1.000 |

The two halves of the squared error pull the interaction trunk and the ligand
encoder in **opposing directions**. This is the mechanism behind every failed
transport in this project, measured rather than inferred. `level_head`'s shape
gradient is 4.7e-8 — numerically zero, as it must be, since `var(p-y)` is
analytically independent of a constant.

Under routing (B2, B3) the same table reads `interaction` level-norm
**0.000e+00**, shape-norm 0.186, and `level_head` shape-norm **0.000e+00**. The
routing does exactly what it claims, verified at the gradient level.

## 2. Three-seed results, double-cold `meta_val`

| k=0 | MSE | < 0.40 tier | wrong-protein | CI | Spearman | ligand-only |
|---|---:|---:|---:|---:|---:|---:|
| **A0** incumbent | 2.1488 | 2.1033 | 2.2292 | **0.5800** | **0.2228** | 2.2742 |
| **B1** factorized, ordinary | 2.4171 | 2.3844 | 2.5120 | 0.5190 | 0.0745 | 2.5956 |
| **B2** + routing | 2.3228 | 2.2641 | 2.5687 | 0.5243 | 0.0686 | 2.3344 |
| **B3** + counterfactual = full | **2.0554** | **1.9670** | 2.2654 | 0.5312 | 0.0855 | 2.0312 |

| k=5 | MSE | < 0.40 tier | wrong-protein | CI | Spearman |
|---|---:|---:|---:|---:|---:|
| A0 | 0.9154 | 0.9237 | 1.0451 | 0.6293 | 0.3557 |
| B1 | 0.9568 | 0.9760 | 1.1986 | 0.6188 | 0.3289 |
| B2 | 0.9385 | 0.9495 | **1.4552** | **0.6336** | **0.3606** |
| B3 | **0.9023** | **0.9058** | 1.2728 | 0.6163 | 0.3182 |

## 3. The training innovation is real and is a major performance source

Paired component bootstraps along the ablation ladder, **identical architecture**
throughout, so the only variable is the training method (`LADDER_v2_3seed.json`):

| contrast | k=0 MSE | k=0, `< 0.40` tier |
|---|---|---|
| **B3 vs B1** (full method vs ordinary) | **+0.3618 [+0.0275, +0.7002]** | **+0.4174 [+0.1248, +0.7308]** |
| B2 vs B1 (routing alone) | +0.0944 [-0.1556, +0.3597] | +0.1203 [-0.1075, +0.3703] |
| B3 vs B2 (counterfactual alone) | +0.2674 [-0.0695, +0.6535] | +0.2971 [-0.0077, +0.6698] |

**The complete method reduces k=0 MSE by 15.0% and low-similarity-tier MSE by
17.5% against conventional training of the same model, with component-level
lower bounds above zero.** It also clears the bootstrap at k=5 (+0.0545
[+0.0194,+0.0917]) and on the k=3 tier.

Gate **T4 passes decisively**: the training method contributes **+0.362** at k=0
while switching architecture contributes **-0.268**. It is not a cosmetic
regularizer — it is the only thing making the new architecture viable at all.

Gates **T1 and T2 fail**: neither half is separately resolved at k=0. Routing
alone and counterfactual alone each move the point estimate in the right
direction, and only their combination clears a component bootstrap. The honest
claim is about the method as a whole, not about either half.

## 4. The factorized architecture is rejected on performance

**B1 vs A0 at k=0: -0.2683 [-0.4823, -0.0592]** — the factorized architecture
under ordinary training is *significantly worse* than the incumbent, with a
resolved interval. It is also worse at k=2 (-0.0939 [-0.156,-0.038]) and k=5
(-0.0414 [-0.068,-0.016]).

B3 recovers to +0.0935 [-0.1418, +0.3138] against the incumbent — a 4.3% point
gain that does **not** clear the bootstrap. Innovation A is therefore **not
admitted as a performance contribution**. Its role is structural: it is what
makes the routing expressible at all, and it delivers the protein specificity
below. Its 19 structural gates all pass; that is a correctness result, not a
performance result.

## 5. The first resolved protein specificity in this project

Wrong-protein donor = the **most similar** target from a different homology
component **within the evaluation split**, chosen by Stage R2's `esm_whitened`
metric. Both proteins are equally unseen, so the contrast varies identity alone.

| arm | k=2 gap | k=5 gap |
|---|---|---|
| A0 incumbent | +0.1016 [-0.0060, +0.2234] | +0.1297 [+0.0406, +0.2347] |
| B1 | +0.1614 [-0.0243, +0.3740] | +0.2419 [+0.0780, +0.4432] |
| **B2** | **+0.4216 [+0.1180, +0.7436]** | **+0.5167 [+0.2443, +0.8230]** |
| **B3** | **+0.2697 [+0.1227, +0.4495]** | **+0.3705 [+0.2051, +0.5774]** |

B2's k=2 gap is **4.2x** the incumbent's, with a lower bound above zero where the
incumbent's crosses it. At k=0 every arm's gap crosses zero, so **Z5 fails at
k=0** and the specificity claim is limited to k >= 2.

An earlier version of this control drew donors from `meta_train` and produced
*inverted* gaps for every arm. That was a defect in the control, not in the
models: a `meta_train` donor is a protein the model has fitted, so it confounds
wrong identity with seen-versus-unseen. Both versions are retained.

## 6. What this falsifies

The cycle's central hypothesis was that the level/shape gradient conflict was
what prevented the trunk from learning within-target ordering. The conflict was
real (cosine -0.334) and removing it produced a resolved MSE gain. **It did not
unlock ordering.**

| k=0 | MSE | = calibration | + shape | query spread |
|---|---:|---:|---:|---:|
| A0 | 2.1488 | 1.2358 | **0.9130** | 0.1317 |
| B1 | 2.4171 | 1.4897 | 0.9274 | 0.1154 |
| B2 | 2.3228 | 1.3894 | 0.9333 | 0.1204 |
| B3 | 2.0554 | **1.1309** | 0.9244 | 0.1435 |

Every point of B3's advantage is **calibration** (1.4897 -> 1.1309); shape is
flat and marginally worse than the incumbent's. Removing the interaction branch
entirely changes k=0 MSE by **-0.0241 [-0.0587, +0.0024]** for B3 — the branch
contributes nothing, and its sign is negative.

So the protein specificity in section 5 lives in `target_level`, not in
`centered_interaction`: swapping the protein moves the level, and the level is
what the routing improved. **The gradient conflict was a real and fixable
problem, and it was not the binding constraint on within-target shape.** That is
the sharpest result of this cycle and it closes the hypothesis.

## 7. Gate outcomes

**Zero-shot** — Z1 fail (4.3% against a 10% target), Z2 fail (interval crosses
zero), **Z3 fail** (CI regresses by 0.049 [0.010, 0.090] — resolved), Z4 fail
(2 of 3 seeds), Z5 fail at k=0 / pass at k >= 2, Z6 fail (tier gain 6.5%,
unresolved).

**Few-shot** — F1 **fails by construction**: a softmax over one support is
identically 1, so the k=1 correction is a scalar level shift and its concordance
equals the k=0 concordance to four decimals. F2 **passes** for every arm (full
beats the level-only transport, all lower bounds above zero, e.g. B3 k=5 +0.2582
[+0.155,+0.367]). F3 **passes** for every arm (permutation gaps +0.32 to +0.51,
all lower bounds above zero). F4 fails against the incumbent. F5 **passes** —
the resolved level-only gap shows the transport contributes beyond the endpoint.

**Training innovation** — T1 fail, T2 fail, **T3 pass for the combined method**,
**T4 pass**.

Since the mandate authorises extending the few-shot core only *after* zero-shot
admission, and zero-shot admission did not occur, the few-shot core is left
unchanged and F1 is recorded as failed by construction rather than by
experiment.

## 8. Two defects the gates caught, both retained

**Structural unidentifiability.** The centering makes any constant introduced
after the last nonlinearity cancel in `s(P,L) - mean_m s(P,anchor_m)`. The
gradient-coverage gate flagged `mix.3.bias` and `readout.bias` with exactly zero
gradient; both were removed. Removing them is required by the centering, not a
tuning choice.

**An identifiability defect the routing created.** The first routed ladder
(`A1/A2/A3_*_predrift`, retained) showed k=0 falling 2.3704 -> 1.9482 -> 1.9086
together with a `ligand_only` MSE of **175.8**. Routing leaves one direction
free: a constant added to `centered` within a target is detached out of the level
term and invisible to the shape term, which is a variance. It drifted and
`target_level` silently compensated. The fix is a **label-free** identifiability
constraint `mean_q(centered)^2`, active exactly when routing is; it contains no
label so it cannot buy fit. The gradient table shows it is *required*, not
optional: under routing the anchor parameter receives **0.000e+00** from both
routed terms — analytically so, since the centering is constant across queries
and therefore invisible to a variance — and without this constraint the anchors
are a dead parameter. With it, `ligand_only` returns to within 0.003 of the
endpoint and the reported gains fall to their honest size.

**A capacity defect.** The first `TypedLigandChannels` pooled each ligand to five
vectors before any protein contact, against an incumbent that attends per atom;
arm A1 reached k=0 concordance 0.4634 with Spearman -0.0948, worse than chance.
Three learned query slots per pharmacophore type (16 tokens, anchors in the same
space) recovered it. One changed variable; both ladders are retained.

## 9. Resources

1,915,880 trainable parameters. Peak CUDA memory **401 MB** against an 8,188 MB
device. 120-153 s per 1,200-step run on an RTX 4060 Laptop. The factorization
makes the wrong-protein counterfactual cheap: `encode_ligand` is protein-blind,
so a counterfactual protein re-runs only the protein encoder, the level head and
the interaction — which is why it fits in every training step.
