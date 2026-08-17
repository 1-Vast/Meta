# Stage B: residual-complementary partial meta-adaptation — rejected, with a cause

**Decision: the AdaMBind-inspired target-task meta-adaptation framework is NOT
admitted to the production MetaSieve model.** No code was promoted to `model/`
or `scripts/`.

Framework inspiration only. This is **not** a reproduction of AdaMBind, MAML,
ANIL, Meta-SGD or FS-CAP, and nothing here evaluates those methods.

Single seed, `meta_val` read once, **development evidence only**. `meta_test`
was not read; its audited status is *logical exclusion after parsing with an
open process-isolation incident* — not "untouched".

Authority: `RESULT.json`, `STAGE1_meta_val.json`, `STAGE2_meta_val.json`,
`LEAKAGE_meta_val.json`, `CONDITIONING_STAGEB.json`. Preregistration frozen
before any arm trained; `CORRECTION_AUDIT.md` records the eight Stage A defects
this stage corrects. Stage A's artifacts are preserved unmodified.

---

## 1. The largest measured effect in this cycle is a leak, not a mechanism

`Tleak` is `T` with one difference: checkpoints selected on `meta_val` — the
rule Stage A and the incumbent trainer use — on the **same** 227 fit components.

| k | `Tleak` − `T` MSE | resolved |
|---|---|---|
| 0 | −0.6180 [−1.6622, +0.1625] | no |
| 2 | −0.1741 [−0.3875, −0.0092] | **yes** |
| 3 | −0.1086 [−0.2398, −0.0016] | **yes** |
| 5 | −0.0839 [−0.1638, −0.0133] | **yes** |

k=0 MSE 2.7425 → **2.1246**; k=0 Pearson 0.0561 → **0.1557**.

Against Stage A's `A0` (2.0753, full `meta_train`, `meta_val`-selected), the gap
to leak-free `T` (2.7425) is 0.6672 pK², and it decomposes:

| source | pK² | share |
|---|---:|---:|
| **`meta_val` checkpoint selection** | **0.6180** | **93%** |
| 12% smaller training set (227 vs 258 components) | 0.0492 | 7% |

**The leak is 5.6× the largest mechanism effect measured in this cycle (0.111)
and 10.7× the recorded same-config retraining spread (0.058).** Every `meta_val`
figure produced by the standard trainer — including the recorded incumbent band
— is optimistic by roughly this margin. This is a measured fact about the
evaluation protocol, not about any method.

## 2. Root cause: the ligand representation collapses within a target

Mean pairwise cosine between the readout hidden vectors of a target's query
ligands, `meta_val`, k=5:

| arm | cos(hᵢ, hⱼ) | ‖h − h̄‖ / ‖h̄‖ | same for `embed` |
|---|---:|---:|---:|
| `T` | **0.99859** | 0.03408 | 0.03292 |
| `M` | 0.99784 | 0.04080 | 0.03292 |
| `H` | 0.96309 | 0.17903 | 0.03172 |
| `C` | 0.99659 | 0.05221 | 0.03383 |

Within one target, every ligand produces almost the same readout activation.
That is not created by the readout MLP — `embed`, its input, is already just as
collapsed (0.033 vs 0.034).

The algebra follows immediately. A weight update `dw = −lr·Σᵢ cᵢ hᵢ` moves query
`q` by `⟨dw, h_q⟩`. With `hᵢ ≈ hⱼ` that is **the same number for every query — a
level shift**. And if the target `c` is mean-zero, `Σᵢ cᵢ hᵢ ≈ (Σᵢ cᵢ)·h̄ = 0`
and **the adapter produces nothing at all**.

This one measurement explains every arm's behaviour:

| arm | mean \|meta\| | mean \|centered shape\| | share that is level |
|---|---:|---:|---:|
| `M` (k=5) | 0.6542 | 0.0021 | **99.7%** |
| `H` (k=5) | 0.1482 | 0.0004 | **99.7%** |
| `C` (k=5) | 0.0389 | 0.0002 | inert |

*(normalized label units)*

## 3. `H` — the naive hybrid triggers the stop pattern

| k | MSE | Spearman | CI |
|---|---|---|---|
| 0 | −0.2252 [−0.6415, +0.1538] | −0.0277 | −0.0124 |
| 1 | −0.2802 [−0.7669, +0.0742] | −0.0263 | −0.0118 |
| 2 | **−0.2077 [−0.4706, −0.0104]** | −0.0377 | −0.0159 |
| 3 | −0.1049 [−0.2788, +0.0338] | −0.0203 | −0.0052 |
| 5 | −0.0593 [−0.1857, +0.0417] | −0.0235 | −0.0054 |

Mean MSE gain +11.9% across k>0, resolved at k=2, surviving the low-recall
stratum. **And ranking degrades at every single k.** That is the preregistered
stop condition, and the mechanism explains it exactly: a 99.7%-level correction
lowers squared error by recalibrating the target mean and cannot improve
ordering.

The gain also *shrinks* with support (0.280 → 0.059), which is the
residual-duplication the whole stage was built to test — confirmed as a
diagnosis.

## 4. `C` — the candidate improves both metrics, and it is still a rejection

`C` is the only arm in this project's record to improve MSE **and** ranking
together with resolved intervals:

| k | MSE | Pearson | Spearman | CI |
|---|---|---|---|---|
| 0 | −0.2220 [−0.5689, +0.0505] | **+0.0529 R** | **+0.0920 R** | **+0.0363 R** |
| 1 | **−0.1364 R** | +0.0529 R | +0.0920 R | +0.0363 R |
| 2 | **−0.0754 R** | +0.0403 | **+0.0618 R** | +0.0228 |
| 3 | −0.0304 | +0.0273 | **+0.0627 R** | **+0.0267 R** |
| 5 | −0.0256 | +0.0217 | **+0.0495 R** | **+0.0228 R** |

**But none of it comes from the meta-adapter.** Three independent facts:

1. `C`'s correction is inert: mean |meta| = 0.0000 (k=1), **0.0006** (k=2),
   0.0434 (k=3), 0.0389 (k=5).
2. Its no-adaptation control is +0.0004 / −0.0360 / −0.0159 — removing the
   adapter changes nothing, or slightly *helps*.
3. **The `C` − `T` ranking contrast is bitwise identical at k=0 and k=1**
   (+0.052852 Pearson, +0.092001 Spearman, +0.036328 CI). At k=1 `C`'s meta is
   exactly zero and the transport is a pure level shift, which cannot reorder
   anything — so the ranking gain must already exist at k=0.

**`C`'s entire advantage over `T` is a zero-shot trunk difference produced by
*training with* the complementary objective, not by the few-shot mechanism it
deploys.** It is a training-time representation effect. Reporting it as evidence
for meta-adaptation would be wrong.

### The preregistered gates

| gate | verdict | value |
|---|---|---|
| G1 mean MSE gain ≥ 5% | **FAIL** | +4.89% |
| G2 non-decreasing with support | **FAIL** | 0.1364 → 0.0256 |
| G3 k=0 degradation ≤ 1% | PASS | −8.09% (improves) |
| G4 no material ranking loss | PASS | Spearman +0.049…+0.092 |
| G5 beats transport-only **and** meta-only | **FAIL** | loses to `M` at k=1, k=2 |
| G6 depends on correct support more than `T` | **FAIL** | −0.21 / −0.42 / −0.65 / −0.79 |
| G7 survives low-recall stratum | PASS | none resolved |
| G8 no adaptation overshoot | PASS | α = 0.24–0.88, 0% overshoot |

## 5. What the counterfactuals actually show

Stage A reported a large correct-vs-wrong-support gap for `A1` and attributed it
to adaptation. With the corrected pre-adaptation anchor and **every arm
measured**, the baseline shows the same thing:

| k | `T` matched-wrong − correct |
|---|---|
| 1 | +2.4051 [+0.8774, +4.7312] |
| 5 | +5.8874 [+2.5703, +10.6243] |

The incumbent Tanimoto transport is label-driven *by construction*, so this gap
was never evidence for the inner loop. The correct statistic is the incremental
dependence, and for `C` it is **negative at every k** (−0.21 to −0.79): the
candidate depends on correct support *less* than the baseline does.

## 6. Decision and stop rules

Three preregistered stop conditions fired:

- the complementary arm does not beat both constituent mechanisms;
- improvement is only level calibration (`H` and `M`, both 99.7%);
- ranking degrades while MSE improves (`H`, every k).

**Rejected. Nothing promoted. No multi-seed run, no Davis, no KIBA** — the
single-seed gate is the precondition for those and it did not pass.

Per the governing contract, no rescue: no attention blocks, no full-backbone
MAML, no learned task selector, no extra datasets. The measured obstacle is
representational collapse *upstream* of the adapter, and no adaptation rule
downstream of a collapsed representation can repair it.

## 7. What is worth keeping

1. **The leakage measurement.** `meta_val` checkpoint selection is worth ~0.62
   pK² at k=0. Any future comparison must remove it, and prior `meta_val`
   numbers should be read with that margin in mind.
2. **The collapse measurement.** Within-target ligand cosine of 0.997 bounds
   what *any* partial-head adaptation can achieve, and connects Stage R's inert
   operator (query spread 0.0027 pK), Stage P's uninformative protein response
   and Stage L2's weak directional signal to a single upstream cause.
3. **The `C` trunk effect** — training with a complementary-residual objective
   improved zero-shot Spearman by +0.092 (resolved). That is a *training-time*
   observation deserving its own preregistered test, and explicitly not evidence
   for meta-adaptation.

## Commands

```bash
conda run -n drug python -m tools.research.stageB_complementary.train_stageb --mode C --steps 1200 --seed 20260815 --output report/meta_fewshot/stageB_complementary_20260817/C
```

```bash
conda run -n drug python -m tools.research.stageB_complementary.evaluate_stageb --stage report/meta_fewshot/stageB_complementary_20260817 --arms T M H C --output tools/research/stageB_complementary/STAGE2_meta_val.json
```

```bash
conda run -n drug python -m tools.research.stageB_complementary.evaluate_stageb --stage report/meta_fewshot/stageB_complementary_20260817 --arms T Tleak --allow-mixed-selection --output tools/research/stageB_complementary/LEAKAGE_meta_val.json
```

```bash
conda run -n drug python -m pytest tools/research/stageB_complementary/tests -q
```
