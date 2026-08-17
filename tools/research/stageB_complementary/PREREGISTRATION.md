# Stage B preregistration: residual-complementary partial meta-adaptation

Frozen 2026-08-17, before any Stage B arm was trained. Stage A's artifacts are
preserved unmodified; `CORRECTION_AUDIT.md` records the eight defects this
stage corrects.

**Framework inspiration only.** This is not a reproduction of AdaMBind, MAML,
ANIL, Meta-SGD or FS-CAP, and no result here evaluates any of them.

## Borrowed and rejected, stated explicitly

| source | principle borrowed | rejected here |
|---|---|---|
| **AdaMBind** | one protein target = one task; disjoint same-target support/query; inner adaptation on support, outer optimization on query loss | its adaptive task selection — **measured and rejected in Stage A**, where it selected candidates with support/query gradient cosine +0.9897 against a +0.6555 population mean, i.e. the redundant tasks; also its label-noise injection, not tested here |
| **MAML** | an initialization trained so that a few gradient steps adapt well | second-order meta-gradients and full-backbone adaptation, both excluded as first choices |
| **ANIL** | adaptation confined to the task-specific head | nothing; this is the core borrowed idea |
| **Meta-SGD** | a learned update magnitude | its full per-parameter learning-rate vector — 97 extra degrees of freedom fitted on as few as two support points. Replaced by **two bounded scalars** |
| **FS-CAP** and few-shot bioactivity models | query-conditioned use of same-assay support compounds | large-scale multi-assay pretraining and any cross-dataset support memory |

**AdaMBind's own novel-target few-shot results are mixed** — on BindingDB and
Davis some correlations improve while MSE/CI degrade. That is precisely why
admission here is decided by local evidence against local gates, and why a
correlation-only improvement is not sufficient for admission.

## Central hypothesis

Stage A's inner loop fitted the **same** support residual the incumbent
Tanimoto transport already explains, so the two mechanisms substituted rather
than composed. That is the measured reason its gain *shrank* as support grew
(0.111 → 0.098 → 0.044 → 0.021 pK² at k=1/2/3/5) instead of growing.

An adapter trained only on the residual that target-level calibration and
Tanimoto transport do **not** explain may add complementary few-shot signal
while leaving zero-shot untouched.

## Arms

| arm | inner loop | inner target | transport | role |
|---|---|---|---|---|
| `T` | none | — | incumbent | corrected `A0` |
| `M` | weight+bias | raw support label | **off** | meta alone |
| `H` | weight+bias | raw support label | incumbent | corrected `A1` (naive hybrid) |
| `C` | **weight only** | zero-shot + **complementary** residual | incumbent | the candidate |

Verified identity: the incumbent transport `shrink * Σ w_qk r_k` equals
`shrink * (mean(r) + Σ w_qk (r_k − mean(r)))` because softmax weights sum to
one. So `C = T + meta_correction` **exactly** — the contrast isolates one
additive term. Pinned by `test_C_equals_T_plus_the_meta_correction`.

`C` adapts the weight only: the level already has an explicit closed term, and
letting the bias adapt too would fit the same level twice.

## The residual decomposition

```
r_i = y_i − zero_shot_i                    support residual
L   = mean(r)                              target level        (raw mean)
t_i = leave-one-out Tanimoto( r − L )_i    neighbourhood transfer
c_i = r_i − L − t_i                        complementary — the adapter's target
```

Leave-one-out is load-bearing: a support item is its own Tanimoto-1.0
neighbour, so an ordinary support-to-support transport would predict `r_i` from
`r_i` and drive `c_i` to numerical noise. The diagonal is masked before the
softmax.

**At k = 1 this yields `c ≡ 0` exactly** — `centered = r − mean(r) = 0` and there
is no other support item. The adapter is therefore *inert at k=1 by
construction*, which is the honest encoding of "one support label cannot
identify within-target shape". `C` makes **no structural SAR claim at k=1**.

At k = 2 the leave-one-out weight is forced onto the single other item, so
`c = 2 × centered` — a rescaling, not new information. Genuine neighbourhood
weighting begins at k ≥ 3.

## Leakage controls, all frozen before training

- **checkpoint selection never reads `meta_val`.** `meta_train`'s homology
  components are partitioned once (seed 20260818) into 227 fit and 31
  internal-validation components. Training episodes come only from fit
  components; selection reads only the internal bank. Partitioning is by
  component, never by target.
- `meta_val` is read **exactly once**, after the candidate is frozen.
- normalization (`training_label_scale`) is fitted on `meta_train` only.
- inner learning rate selected on `meta_train` folds (Stage A's
  `INNER_LR_SELECTION.json`, 0.1 at one step); not re-tuned here.
- query labels are loss-only; structurally excluded from prediction,
  adaptation and inference, and tested by perturbation.
- `stable_seed` everywhere; never Python `hash()`.
- nested support banks with an identical query panel across k (verified).
- `meta_test`: **logical exclusion after parsing, process-isolation incident
  open.** Not "untouched". It may be evaluated exactly once, only after
  architecture, hyperparameters, seed set and thresholds are frozen.

## Metrics and controls

Reported per k ∈ {0,1,2,3,5}: MSE, RMSE, CI, Spearman, Pearson, R², centered
within-target MSE, and activity-cliff sign accuracy where powered.

Controls, **for every arm**: no adaptation; correct support; permuted support;
matched-wrong support anchored to the **pre-adaptation** support prediction;
wrong-protein (a **full-system** perturbation, not an inner-loop test);
level-only and shape-only decompositions of the query correction; zero-shot
identity at k=0.

Strata: ligand novelty by max Tanimoto to `meta_train`, with a low-recall
stratum (< 0.4) so any benefit can be separated from near-duplicate recall.

## Gates

### Stage 1 stop rule

If corrected `H` loses its Stage A k=1/k=2 direction against `T`, the framework
closes and no complexity is added.

### Stage 2 single-seed screen — all required to advance

1. mean MSE across k=1/2/3/5 improves by **≥ 5%** vs `T`;
2. benefit is strongest or non-decreasing as support grows;
3. k=0 MSE degrades by **≤ 1%**;
4. no material CI or Spearman degradation;
5. `C` beats **both** `T` and `M`;
6. `C` depends on correct support binding **more than `T` does**
   (difference-of-differences on the permuted and matched-wrong controls);
7. benefit survives in the low-recall novelty stratum;
8. finite, stable inner gradients; no adaptation overshoot (`alpha < 1`).

### Final admission — multi-seed, ≥ 3 matched seeds

Same-direction improvement across seeds; a hierarchical seed/component interval
excluding zero on the primary MSE contrast; no significant CI/Spearman
degradation; correct-support advantage over permuted and matched-wrong; k=0
non-inferiority; and no benefit attributable to `meta_val` selection, query-label
leakage, changed banks or extra data.

**A component bootstrap is within-checkpoint uncertainty and is never a
substitute for retraining seeds.** The recorded same-config retraining spread is
0.058 pK² in k=0 MSE.

## Promotion

Promotion to `model/` and `scripts/` requires the multi-seed BindingDB gate, not
the single-seed screen. It does not become the default until independent Davis
or KIBA evidence reproduces the direction. Each dataset is trained from scratch;
labels and examples are never pooled.

## Data

Exactly one dataset: governed BindingDB-Ki `bindingdb_ki_double_cold_v1`. No
Davis, KIBA, ChEMBL, structural corpus, cross-dataset support, retrieval,
normalization statistic, label or checkpoint.
