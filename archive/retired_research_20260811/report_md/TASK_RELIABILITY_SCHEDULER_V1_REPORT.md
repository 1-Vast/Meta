# AdaMBind Task Reliability/Transferability Scheduler Audit

**Date:** 2026-08-11
**Status:** formal CUDA development Gate completed; scheduler integration
rejected; no `meta_test` or fresh-confirmation result was used for fitting,
scoring, selection or stopping.

## Decision

```text
GATE0_STATIC_AVAILABILITY_PASS
GATE1_TASK_SCHEDULER_IDENTIFIABILITY_FAIL
REJECT_TASK_SCHEDULER_GATE1_FAIL_CLOSED
SHORT_CUDA_TRAINING_AUTHORIZED=false
TRAIN_MAIN_V1_INTEGRATION_AUTHORIZED=false
BIOLOGY_ADMISSION_AUTHORIZED=false
```

AdaMBind remains a core methodological reference for target-as-task episodic
learning and for the hypothesis that source tasks should not be sampled
uniformly. Its public scheduler/noise implementation is not transplanted. The
current MetaSieve loss/alignment statistics do not identify partner-specific
transfer utility after dependency-aware cross-fitting and destructive
controls.

## Primary-source and official-code audit

Sources were the [Nature Communications article](https://www.nature.com/articles/s41467-026-70554-5),
its complete 19-page Supplementary Information (84 paragraphs, 22 tables), and
official repository commit
[`01a169a6`](https://github.com/Moohyun-w/AdaMBind/tree/01a169a6d62fba0d6c003f47bfba539e55f5b344).

The paper-level adaptive module takes query loss, layer-wise support/query
gradient agreement and progress, then describes a separate-validation
approximate bilevel resampling cycle. This supplies optimization information,
not label reliability or a partner-specific biological observable. The
reported uniform label perturbation grid (`0.1` to `0.6`) is not based on
replicate disagreement or assay uncertainty.

The official implementation differs materially from that description:

- `train.py` concatenates `train_idxs + val_idxs` and samples scheduler-train
  and scheduler-validation tasks from the same pool, so their roles can
  overlap.
- The adaptor update is a no-baseline REINFORCE loss, not the differentiable
  one-step scheduler update depicted in Supplementary Algorithm 1.
- `sample_task(..., replace=True)` is the default, allowing duplicate tasks in
  one selected meta-batch.
- Both LSTMs receive sequence length one (`reshape(1, tasks, features)`), so
  bidirectionality does not implement a task curriculum sequence.
- The advertised set context sums feature coordinates within each task
  (`sum(dim=1)`) rather than aggregating across tasks; it is not a leave-one-task
  set context.
- `adaptive_tasks` is parsed but never gates the adaptive path. With
  `noise=0`, local training target `y` is undefined; query perturbations are
  drawn in one path but the reported query MSE still uses clean labels.
- `Trainer.gradients` is persistent and is not reset between `train` calls,
  so the meta-gradient accumulator carries prior calls forward.
- Support/query construction is a shuffled per-target prefix with no scaffold,
  assay/document or recipient closure.

These discrepancies make the official implementation `REFERENCE_ONLY`; direct
reuse is `REJECT`.

## MetaSieve adaptation

The frozen preregistration is
`research/meta_fewshot/PREREG_TASK_RELIABILITY_SCHEDULER_V1.md`. The retained
V1 `uniform_clean` checkpoints generate source and meta-validation episode
gradients on CUDA. A two-input ridge scorer uses only:

1. `log1p(clean query MSE)` as difficulty;
2. support/query gradient cosine as local directional compatibility.

Whole CD-HIT40 components define five outer cross-fitting folds and nested
regularization folds. The target is mean first-order gradient utility to
correct, isolated `meta_val` components. The equal-capacity null has the same
two inputs, folds, ridge grid and compute, but permutes the two statistic rows
within task-size strata.

Task size, label spread, scaffold overlap and provenance density are nuisance
audit columns only; they never enter scorer fitting. Empty Murcko strings are
missing values, not a shared scaffold. They are excluded from overlap; an
all-missing task receives mean imputation only for nuisance residualization and
an explicit missing-fraction indicator. `replicate_count` and `panel_count`
are provenance density, not empirical reliability. Replicate disagreement and
continuous protein familiarity are `NA`.

The already-fitted out-of-fold scorer is re-evaluated, never refit, under
protein shuffle, wrong-target support, cross-task label permutation,
ligand-only and intercept-only destruction. These outcomes are not scheduler
supervision.

## Gate results

Gate 0 passed with 285 `k=5` source targets in 207 components and 37
meta-validation targets in nine components. The largest validation component
share was `0.243`, below the frozen `0.35` ceiling.

| k | clean score-utility rho | 95% LCB | correct-minus-null rho | 95% LCB | protein-shuffle utility delta | 95% LCB |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | -0.0597 | -0.1574 | +0.0522 | -0.0627 | +0.000329 | -0.001210 |
| 2 | -0.1635 | -0.2642 | -0.0876 | -0.1909 | +0.000196 | -0.001526 |
| 3 | +0.0102 | -0.0981 | -0.0694 | -0.2000 | +0.001407 | -0.000509 |
| 5 | -0.0885 | -0.2057 | -0.2091 | -0.3556 | +0.000922 | -0.001290 |

The absolute utility criterion (`rho >= 0.20`, LCB > 0) failed at every k.
Matched-null superiority failed at every k. Protein-shuffle necessity failed
at every k; score correlation often increased rather than decreased after
protein or wrong-support destruction. Label and shortcut controls were also
not all-k stable. No ordinary MSE result can override these failures.

## Interpretation

The nuisance-aware result is more negative than the earlier online ATS result.
Loss/alignment can track properties shared by source and validation episodes,
but that relation disappears out of component and does not require the correct
protein representation. This is a curriculum/calibration shortcut, not a
transferable partner-specific scheduler signal.

The code adaptation is retained as a reusable falsification harness and test
suite. Because Gate 1 failed, the main trainer was deliberately not changed,
no short three-seed scheduler training was run, and the fixed label-noise lane
remains separately rejected. Improving Cold Target performance still requires
a new affinity-directed partner observable; task scheduling cannot repair the
F-132 identifiability failure.

The authoritative artifact is
`report/meta_fewshot/task_reliability_scheduler_v1/RESULT.json`. It binds the
preregistration, runner, scorer, model and three checkpoints by SHA256 and
records PyTorch `2.6.0+cu124`, the RTX 4060 Laptop GPU, and 131.05 seconds of
formal audit compute.
