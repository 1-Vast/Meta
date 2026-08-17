# Stage A: target-task inner/outer-loop meta-learning, single-seed screen

Verdict: **NOT PROMISING** on the preregistered conjunction (4 of 6 gates pass).

* the **inner loop (`A1`)** is a *weak positive that did not resolve* — better
  than the baseline at every k on every metric, never with an interval
  excluding zero. Retained as a three-seed candidate.
* the **task selector (`A2`)** is **rejected** — worse than `A1` on MSE,
  Pearson, Spearman and CI at every k.

AdaMBind-**inspired** framework test. This is **not** a reproduction of
AdaMBind and nothing here evaluates that method. **Single seed, development
evidence, directional screening only** — not a performance claim. `meta_test`
was not read.

Authority: `STAGE_A_meta_val.json` (+ `.rows.jsonl`), `RESULT.json`.
Preregistration: `PREREGISTRATION.md`, frozen before any arm was trained.
Population: double-cold `meta_val`, 41 targets / 19 components, 82 episodes per
k, k ∈ {0,1,2,3,5}, component-paired bootstrap (9,999 draws).

## What was built

One code path serves all three arms, so `--inner-steps 0` *is* the accepted
recipe and `A0` is matched by construction. Both halves of that are tested:

| identity | result |
|---|---|
| re-implemented readout vs `InteractionGrammarModel.encode` | **bitwise equal** at k ∈ {0,1,3,5} |
| `A0` episode loss vs the production episode loss | equal to **2e-5 on all seven terms**, each verified non-zero |

Adaptable scope: `interaction_head.2.{weight,bias}` — **97 parameters, 0.0054%
of the 1,798,833 trainable**. Functional fast weights: no module is ever
written to, so one task cannot mutate the persistent model or another task's
state. First-order (`create_graph=False`). Inner step size 0.1 at 1 step, chosen
on `meta_train` component folds against the frozen A0 checkpoint;
`meta_val` was not read for it.

## Results

### Per-arm, each at its own operating condition

| k | arm | MSE | Pearson | Spearman | CI |
|---|---|---:|---:|---:|---:|
| 0 | A0 | 2.0753 | 0.1735 | 0.1689 | 0.5484 |
| 0 | **A1** | **2.0579** | **0.2069** | **0.1981** | **0.5673** |
| 0 | A2 | 2.1405 | 0.1038 | 0.0685 | 0.5170 |
| 1 | A0 | 1.5352 | 0.1735 | 0.1689 | 0.5484 |
| 1 | **A1** | **1.4242** | **0.2084** | **0.1991** | **0.5676** |
| 1 | A2 | 1.4512 | 0.1051 | 0.0689 | 0.5171 |
| 2 | A0 | 1.1718 | 0.2459 | 0.2291 | 0.5756 |
| 2 | **A1** | **1.0738** | **0.2597** | **0.2525** | **0.5875** |
| 2 | A2 | 1.0844 | 0.2142 | 0.1812 | 0.5653 |
| 3 | A0 | 1.1245 | 0.2768 | 0.2641 | 0.5906 |
| 3 | **A1** | **1.0803** | **0.2844** | **0.2798** | **0.5997** |
| 3 | A2 | 1.0998 | 0.2529 | 0.2153 | 0.5790 |
| 5 | A0 | 0.9075 | 0.3291 | 0.3500 | 0.6314 |
| 5 | **A1** | **0.8869** | **0.3307** | **0.3549** | 0.6295 |
| 5 | A2 | 0.8983 | 0.3204 | 0.3048 | 0.6133 |

`A0`'s k=0 MSE of 2.0753 sits inside the recorded incumbent band (frozen
checkpoint 2.1488, `A0repro` 2.0911, retraining spread 0.058), so the baseline
reproduced.

### The gates

| gate | verdict | evidence |
|---|---|---|
| **G1** A1 beats A0 at k>0, k=0 undamaged | **pass, unresolved** | MSE −0.1111 / −0.0980 / −0.0443 / −0.0206 at k=1/2/3/5; k=0 −0.0175 (−0.8%). Every interval crosses zero. |
| **G2** A2 beats A1 | **FAIL** | A2 worse at every k on every metric |
| **G3** MSE and ranking agree | pass | A1 improves both; A2 degrades both |
| **G4** correct support beats wrong support | **pass, resolved** | see below |
| **G5** improvement grows with support | **FAIL** | it *shrinks*: 0.1111 → 0.0980 → 0.0443 → 0.0206 |
| **G6** no leakage | pass | query labels structurally excluded and tested; banks identical; no extra data |

Only k=1 (−0.1111) and k=2 (−0.0980) exceed the recorded 0.058 retraining
spread. k=3 and k=5 fall below it and are not claimed as real.

### Support labels genuinely matter (G4, resolved)

Control minus correct MSE; positive means the control is worse. Every interval
below excludes zero.

| k | permuted support | matched-wrong support |
|---|---|---|
| 1 | — (a 1-element permutation is the identity) | **+0.8099 [+0.1709, +1.5868]** |
| 2 | **+0.4296 [+0.1580, +0.7869]** | **+1.4868 [+0.7063, +2.4267]** |
| 3 | **+0.3616 [+0.1225, +0.6964]** | **+1.6088 [+0.7994, +2.5939]** |
| 5 | **+0.3905 [+0.1983, +0.6009]** | **+2.0516 [+1.1117, +3.2183]** |

This is the strongest part of the result: the mechanism is genuinely
label-bound, not a calibration artifact.

### `A1` also shows the first resolved wrong-protein gap in the record

| k | A1 wrong-protein minus correct | A2 |
|---|---|---|
| 1 | +0.0253 [−0.0009, +0.0513] | −0.0178 |
| 2 | **+0.0188 [+0.0052, +0.0327]** | +0.0004 |
| 3 | **+0.0177 [+0.0079, +0.0282]** | +0.0002 |
| 5 | **+0.0085 [+0.0037, +0.0137]** | +0.0005 |

Small — 0.008 to 0.019 pK² — but resolved at k=2, 3 and 5. Every prior
wrong-protein gap in this project crossed zero. `A2` erases it to four
decimals.

## k=1: shape or level? **Level.**

The 97-parameter scope splits exactly, so the same run answers this with no
extra training.

| condition | k=1 MSE | share of the gain |
|---|---:|---:|
| no adaptation | 1.5377 | — |
| bias only (pure level shift) | 1.4456 | **81%** |
| weight only (pure shape) | 1.4364 | 89% |
| full update | 1.4242 | 100% |

A pure scalar bias update recovers **81%** of the k=1 adaptation gain, and the
shape residual is +0.0214 [−0.0913, +0.1446] — unresolved. **The k=1 effect is
not distinguishable from a scalar recalibration.** The shape share does grow
with support (19% → 39% → 36% → 41% at k=1/2/3/5) but never resolves.

## Why `A0` cannot simply be given an inner loop at test time

The inner-step sweep looks pathological until it is measured. At k=1 with one
support point and squared error, `alpha = 2·lr·‖h‖²` and the post-step residual
is `(1 − alpha)` times the pre-step residual.

| arm | alpha | overshooting (alpha>1) | residual multiplier |
|---|---:|---:|---:|
| **A0** | **1.514** | **100%** | **−0.514** |
| A1 | 0.241 | 0% | +0.759 |
| A2 | 0.216 | 0% | +0.784 |

`A0` overshoots on every single episode, so each step flips the residual sign
while shrinking it. That predicts an alternating sweep, and the observed sweep
alternates exactly: **1.5352 → 2.8626 → 1.4349 → 2.2049** at 0/1/2/3 steps —
odd steps bad, even steps good, magnitude decaying. `A1`'s training pushed
alpha **6.3× down** into the stable contraction region.

**Consequence: the inner loop cannot be bolted onto an already-trained model.
Training with the loop is what makes the loop well-conditioned.**

## Why the selector failed

It worked as designed, and that is the problem. It concentrated on candidates
whose support and query gradients almost perfectly agree — mean cosine
**+0.9897** among selected against **+0.6555** across all 8,640 candidates,
with 82.8% of candidates positive and none exactly zero.

High support/query gradient agreement identifies the tasks where the support
*already* predicts the query: the easy, redundant ones. Training on them cut
effective task diversity to **5.29 of 9** candidates per step, cost **2.33×**
the encoder forwards, and destroyed the protein specificity `A1` had. This is
precisely the failure mode the preregistration named — high weights selecting
easy tasks rather than informative ones.

## Cost

| arm | steps | encoder forwards | peak GPU | added params |
|---|---:|---:|---:|---:|
| A0 | 1200 | 6,480 | 1030 MB | 0 |
| A1 | 1200 | **6,480** | 1030 MB | 0 |
| A2 | 1200 | 15,120 (2.33×) | 1020 MB | 0 |

**The inner loop is free in encoder terms.** Because the adapted scope sits
downstream of the encoder, `A1` re-evaluates 97 parameters on cached features
and uses exactly as many encoder passes as `A0`. Only the selector costs more.
(Wall times — 1793 / 1100 / 1138 s — are not comparable: `A0` ran alongside a
test-suite sweep.)

## Diagnosis of what did not work

Per the preregistration, naming the cause rather than adding complexity:

1. **G5 (gain shrinks with k): competition with the incumbent transport.** The
   `SimilarityTransport` already shifts predictions by Tanimoto-weighted support
   residuals. The inner loop fits the *same* support labels, so as adaptation
   absorbs the residual, `locked = support_y − support_prediction` shrinks and
   transport does less. At larger k transport already captures the support
   information, leaving less for adaptation to add. This was preregistered as
   conflict #1 and is now observed.
2. **k=1 is mostly a level shift.** 97 parameters fitted on one point move
   predominantly in the direction that shifts the mean.
3. **The selector's objective is miscalibrated,** not its implementation.
   Gradient agreement measures redundancy, not informativeness.

Not the cause: parameter scope was sufficient to produce shape (the mechanism
demonstrably reorders queries), the optimization was stable in the trained arms
(alpha 0.24, no overshoot), and the task construction is clean (0 contract
violations in 400 draws).

## Disclosure that travels with every number here

**The trainer selects checkpoints on `meta_val` labels.** The training gradient
never touches `meta_val` and the label scale is fitted on `meta_train`
(verified numerically), but `train()` keeps the state with the best `meta_val`
admission score. It is identical across all three arms, so it cannot
manufacture a between-arm difference — the only quantity this screen decides —
but **every `meta_val` figure above is an optimistic development estimate, not
a held-out one.** It was disclosed rather than repaired because repairing it
would stop `A0` from reproducing the accepted baseline.

## Bugs found and repaired

1. **The inner loop could not run under `torch.no_grad()`** — adaptation *is* a
   gradient computation and evaluation calls it inside `no_grad`. Support
   adaptation at inference was structurally impossible until an explicit
   `enable_grad` scope was added. The most consequential bug in this stage.
2. **Three auxiliary loss terms were silently dropped** from the first trainer.
   Had that survived, any `A1` gain would have been "the inner loop recovers
   what we deleted". Restored and pinned by
   `test_baseline_equivalence.py`.
3. **Python's `hash()` reintroduced** for a permutation seed — the
   PYTHONHASHSEED defect repaired in Stage R. Caught before any run.
4. **`paired()` compared a condition against itself**, making every
   counterfactual an exact zero that read as a clean null.
5. **Two trainers wrote to one output directory** after a killed shell wrapper
   walked on to the next arm. A fail-closed guard now refuses to start on an
   existing `progress.jsonl`.

## Decision

Per the preregistered rule, `PROMISING` required all six gates; two failed, so
the framework is **NOT PROMISING as a whole on this single seed**.

Split by mechanism, as the preregistration requires:

* **Reject the task selector.** Uniform sampling is retained. Do not implement
  the learned bi-level variant — the simple rule did not produce credible
  evidence, and its failure is explained rather than noisy.
* **Retain the inner loop as an unresolved weak positive.** It is directionally
  better than the accepted baseline at every k on every metric, it is free in
  encoder cost, its support-label dependence is resolved and large, and it
  produced the first resolved wrong-protein gap in this project's record. One
  seed cannot separate that from retraining noise. The honest next step is
  three matched seeds of `A0` vs `A1` alone — **not** additional architecture.

## Commands

```bash
conda run -n drug python -m tools.research.stageA_innerloop.audit_dataflow
```

```bash
conda run -n drug python -m tools.research.stageA_innerloop.select_inner_lr
```

```bash
conda run -n drug python -m tools.research.stageA_innerloop.train_meta --arm A1 --inner-steps 1 --inner-lr 0.1 --steps 1200 --seed 20260815 --output report/meta_fewshot/stageA_innerloop_20260817/A1
```

```bash
conda run -n drug python -m tools.research.stageA_innerloop.evaluate --stage report/meta_fewshot/stageA_innerloop_20260817 --output tools/research/stageA_innerloop/STAGE_A_meta_val.json
```

```bash
RUN_SLOW=1 conda run -n drug python -m pytest tools/research/stageA_innerloop/tests -q
```
