# PSEP operator gate — pre-registration

Written **before** `research/psep_operator.py` results were read. Seed 20260802,
reported on `meta_test` components inside `discover`; `validate` and `confirm`
sealed.

## The bar, derived from measurement rather than chosen

All gains are within-document pair concordance **over the support-free base**,
component bootstrap. From `psep_m2_records` (379 components), closed-form ridge
in the fixed 266-d basis:

| k | 1 | 3 | 5 | 10 | 20 |
|---|---:|---:|---:|---:|---:|
| ridge − base | −0.0036 | +0.0007 [−0.0007,+0.0021] | **+0.0019** [−0.0006,+0.0046] | **+0.0056** [+0.0020,+0.0096] | +0.0096 [+0.0044,+0.0153] |

Full-information ceiling (own head fitted on the whole support pool, ~140 labels)
is **+0.0154** over base.

**k=1 is structurally −0.0036**: a single centred support label determines only an
intercept, and within-document concordance is invariant to intercepts, so any
transport-form operator must reproduce the target-agnostic head exactly there.

## Pre-registered success criteria

1. **Primary.** A learned operator at **k=5** must reach a component-bootstrap
   lower bound **> +0.005** over the support-free base — i.e. it must deliver at
   k=5 what closed-form ridge needs **k≈10–20** to deliver. This is the
   "move the label-budget curve left" claim, stated numerically.
2. **Necessary.** Correct support must beat the best of four wrong-support
   controls (random target, protein-hard, chemistry-matched, label-permuted) by
   **> 0.005**.
3. **Necessary.** Must beat intercept calibration, Tanimoto-KRR and fixed-basis
   ridge — each granted its own grid-selected ridge and global transport scale on
   `meta_val` — by **> 0.005**.

## Structural validity checks (not criteria — harness correctness)

- `intercept` must score **exactly 0.000** on within-document concordance at every
  k. If it does not, the metric is not intercept-invariant and the run is void.
- `attention` (pure transport, `Delta = sum_i alpha_i(x_q) e_i`) must score
  **exactly 0.000** at k=1 for the same reason. `cnp`, whose decoder reads `x_q`
  alongside the pooled support, is not so bound — the k=1 gap between them is the
  clean read on whether anything beyond transport is being learned.

## Declared in advance

If no operator clears criterion 1, the verdict is
`NO_OPERATOR_PASSES_ADMISSION_GATE`, and it will **not** be rescued by more
parameters, longer training, extra losses or additional modules. A negative here
is a statement about the substrate's identifiability, already quantified in
`PSEP_CORE_MECHANISM_DECISION_2026-08-02.md`, not about tuning.
