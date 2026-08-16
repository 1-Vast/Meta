# Stage 6 audit addendum

Recorded after the multi-seed jobs completed. `PREREGISTRATION.md` is left
unchanged; this file records what the audit found, what was corrected, and how
the preregistered gates were resolved. No historical outcome is rewritten.

## Findings, verified independently

| # | claim | verified? | invalidates results? |
|---|---|---|---|
| 1 | seed 20260813 shows no F gain on the small automatic meta_test bank | **yes** — k=1 1.642/1.656, k=2 1.264/1.263, k=3 1.159/1.198, k=5 1.118/1.137; F `full` even worse than F `level_only` there | no — 6 episodes, underpowered; the 42-episode re-evaluation reverses it (F beats its own level in 9/9) |
| 2 | selection used ~17 meta_val episodes, automatic evaluation ~6 meta_test episodes | **yes** | no, but all decisions were re-derived on the complete 44-episode meta_val bank |
| 3 | audit used 2048-bit fingerprints, production uses 1024 | **yes** | no — production rerun gives 0.976/0.886/0.747 against 0.976/0.887/0.747 and Pearson -0.349 against -0.347 |
| 4 | the raw-label audit does not match the residual the model transports | **yes** | no — a frozen-checkpoint residual audit was added and confirms the weighting helps residuals on a checkpoint that never saw the mechanism |
| 5 | attribution must be within-checkpoint | **accepted** | reporting restructured; k=0/k=1 explicitly non-attributable |
| 6 | wrong-protein is confounded for this mechanism | **accepted** — Tanimoto consumes no protein information | reframed as a full-system perturbation |
| 7 | "similarity only"/"zero-learning" inaccurate; gamma ~ 7.99 against init 8.0 | **yes** — 7.982 / 7.985 / 7.990 | wording corrected |
| 8 | `use_learned_key=False` leaves dead parameters | **yes** — `transport.key.weight` and `transport.log_temperature` return `grad=None` | **no** — AdamW skips `grad=None` parameters, so completed runs are numerically unaffected |
| 9 | zero-fingerprint wording wrong | **yes** — a zero fingerprint still receives a finite softmax weight | no — 0 of 9,880 ligands are affected; docstring corrected |
| 10 | data-integrity facts | **yes** — 9,880 rows, 9,880 unique keys, 0 parse failures, 9,880 unique canonical molecules, 0 duplicate groups, 17,717 cells | n/a |

## Fixes made, all surgical

* `model/similarity_grammar.py`: freeze (not remove) the two unused tensors when
  `use_learned_key=False`, preserving state-dict keys and initialisation order;
  corrected the module docstring's description of the mechanism.
* `scripts/qpsmp_data.py`: corrected the zero-fingerprint docstring.
* `scripts/audit_transferable_signal.py`: `--bits` / `--radius` flags and a
  `matches_production_pipeline` field.
* New: `scripts/audit_residual_transport.py`, `scripts/stage6_paired_analysis.py`.
* New tests: `test_similarity_only_has_no_dead_trainable_parameter`,
  `test_similarity_only_checkpoints_remain_loadable`.
* All six existing checkpoints verified to load with `strict=True` after the fix.

**No training run was interrupted, altered or discarded. All outputs preserved.**

## Preregistered gates, resolved on the complete meta_val bank

| id | requirement | outcome |
|---|---|---|
| G1 | `full` beats `level_only` by more than arm A's margin at k in {2,3,5} | **pass**, 9/9 seed-k cells, margins 2-14x arm A's |
| G2 | permutation gap exceeds arm A's | **pass**, +0.40 to +0.51 against +0.06 to +0.20 |
| G3 | CI and Spearman at or above `level_only` | **pass**, 9/9 on meta_val and 9/9 on meta_test |
| G4 | k=0 within 0.05 of arm A | **fails as written**; correctly reclassified as non-attributable, since the mechanism is inactive at k=0 |
| G5 | `full` MSE not worse than arm A | **pass on meta_val, fails on meta_test** — the cross-arm contradiction |

## Decision against the stated policy

> Accept H1 only if the multi-seed, full-meta_val paired analysis shows that
> k=2/3/5 full consistently beats the same checkpoint's level-only baseline,
> with MSE improvement and non-degrading CI/Spearman, while label-permutation
> gaps remain positive, with aggregate CI support rather than one seed.

**ACCEPT** — satisfied on `meta_val` (the deciding split) and independently
replicated on `meta_test`: 18/18 positive point estimates, 16/18 component-level
lower bounds above zero, permutation gaps +0.40 to +0.51 in every cell.

Two qualifications on that sentence, both required corrections:

* the component-level bounds are **not** all positive — `meta_val` k=3 CI
  ([-0.015, +0.106]) and `meta_val` k=3 Spearman ([-0.005, +0.289]) cross zero,
  so the ranking claim at that one cell is not established;
* all intervals are **conditional on the three trained seeds**, because
  `scripts/stage6_paired_analysis.py` averages seeds per (component, target)
  before resampling components. Seed variance is not resampled, so these are
  component-level intervals given these checkpoints, not retraining intervals.

> If similarity_only improves label sensitivity but not predictive performance
> across seeds, retain it as a strong classical baseline or diagnostic.

Not the case: predictive performance improves within-checkpoint across seeds on
both splits. However, **superiority over the incumbent `grammar` transport is
NOT established** — cross-arm results contradict across splits with no lower
bound excluding zero. H1 is therefore admitted as a validated mechanism at
k>=2, not as a replacement for the incumbent.

> If it passes robustly, freeze Stage 6 before exploring k=0/k=1 improvements.

Stage 6 is frozen. No k=0/k=1 change and no second architectural change is
bundled here.

## Claims explicitly NOT made

No performance breakthrough, no SOTA, no protein-specific transport, no core
meta-learning advance, no k=0 or k=1 improvement claim, and no claim that this
mechanism is better than the incumbent transport.
