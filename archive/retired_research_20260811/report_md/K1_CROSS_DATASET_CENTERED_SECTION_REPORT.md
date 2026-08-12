# K1 cross-dataset centered-section migration report

Date: 2026-08-11  
Preregistered verdict: `CENTERED_SECTION_CROSS_DATASET_FAIL`  
Migration: `NOT AUTHORIZED`

## Experiment

The K1 contract was frozen before affinity metrics were computed. It used the
public AdaMBind snapshots of BindingDB, Davis and KIBA as three separate
protein-as-task engineering benchmarks. Target split was outcome-blind and
SHA-derived. All datasets used the same Morgan-256 representation, `k=5`, ten
support draws, at most 50 queries and ridge `lambda=1.0`; no dataset-specific
tuning was allowed.

The public benchmark labels are now consumed for this architecture decision.
This experiment is not an independent biological confirmation.

## Results

| Dataset | Eval targets | Population | Intercept | Uncentered | Centered |
|---|---:|---:|---:|---:|---:|
| BindingDB | 94 | 1.7621 | 1.2219 | **1.0746** | 1.1449 |
| Davis | 75 | 0.6095 | 0.7008 | **0.5998** | 0.6881 |
| KIBA | 42 | 0.8218 | 0.8246 | **0.7285** | 0.7999 |

Positive reduction is better for the centered candidate.

| Dataset | Intercept - centered | one-sided LCB | Uncentered - centered | one-sided LCB |
|---|---:|---:|---:|---:|
| BindingDB | +0.0770 | +0.0473 | -0.0703 | -0.1115 |
| Davis | +0.0127 | +0.0086 | -0.0882 | -0.1104 |
| KIBA | +0.0247 | +0.0153 | -0.0713 | -0.0920 |

Equal-dataset standardized pooled reductions:

```text
intercept - centered    +0.03157   LCB +0.02419
uncentered - centered   -0.09011   LCB -0.10627
```

Numerical failure fraction was `0.0` in all three datasets.

## Interpretation

The centered residual branch contains reproducible ligand-specific signal: it
beats a pure support intercept in all three datasets, with positive individual
and pooled lower bounds. The proposed replacement nevertheless fails because
the retained uncentered ridge is better in every dataset, by about 0.07--0.09
target-macro MSE.

This changes the architecture diagnosis:

1. the generic closed-form ridge solver is not the failed module;
2. forcing an unpenalized support mean into the predictor increases 5-shot
   variance, especially where the source population is already strong;
3. the intercept should remain a standing calibration-only control, not replace
   the working predictor;
4. MetaSieve's calibration-dominant v0 result is therefore localized to its
   biological coordinates/T-BASIS, rather than a universal property of the
   ridge section.

The negative result is useful: it prevents a theoretically tidy but empirically
worse solver from being migrated.

## Decision

```text
KEEP_UNCENTERED_POSITIVE_RIDGE
KEEP_SUPPORT_INTERCEPT_AS_REQUIRED_BASELINE
DO_NOT_MIGRATE_CENTERED_KERNEL_SECTION
NEXT_FAILED_MODULE_TO_TEST = BIOLOGICAL_PAIR_REPRESENTATION
```

No file was copied to `model/` or `scripts/`. The next model-changing stage is
R0/R1 admission of a dataset-independent pair representation. It must pass
structural partner controls and measured crossed-affinity tests on more than one
dataset/domain before being connected to the retained uncentered ridge.

Machine-readable evidence:
`report/meta_fewshot/k1_cross_dataset_centered_section.json`.

