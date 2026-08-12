# Preregistration: Task Reliability/Transferability Scheduler V1

**Frozen:** 2026-08-11, before computing any new scheduler score or opening a
new performance result.

## Scientific axis

This experiment changes one axis only: source-task scheduling. The frozen
MetaSieve V1 representation, target-as-task construction, support/query
isolation, support-only Meta-Section, uncentered positive ridge solve, loss,
optimizer and label-noise policy are not changed. Fixed uniform label noise is
a prior negative control and is not part of this experiment.

The AdaMBind/ATS hypothesis is narrowed to the following falsifiable claim:
a task score learned from dependency-isolated development utility can rank
source tasks for transfer better than an equal-capacity row-permuted scorer,
without reducing to task difficulty, ligand familiarity, task size, label
spread, a wrong support set, or a wrong protein representation.

## Data firewall

- Allowed: `meta_train` for candidate tasks and scorer fitting; `meta_val` for
  first-order transfer-utility targets and development-only evaluation.
- Forbidden: `meta_test`, every sealed/fresh confirmation cohort, R0-C labels,
  and any selection based on the final Cold Target test result.
- Protein groups are the governed CD-HIT40 components. Folds and bootstrap
  units are whole `protein_group_40` components.
- Destructive controls are audit outcomes only. Correct-versus-shuffled
  protein performance, wrong-support identity, or any artificial partner
  label must never enter the scheduler's training features or targets.

All output must say `development`, never `confirmation`.

## Gate 0: availability and estimand

Gate 0 is evaluated before fitting a new scheduler.

1. At least eight eligible `meta_val` protein components and 100 eligible
   `meta_train` target tasks must exist at `k=5`.
2. The largest `meta_val` component may contribute at most 35% of eligible
   targets.
3. Every candidate task must expose task size and label spread. Ligand
   scaffold familiarity is usable only if exact Bemis-Murcko scaffolds can be
   derived from admitted ligand records without touching forbidden data.
4. `replicate_count`, `panel_count`, and panel/document identifiers are
   provenance-density covariates, not measured label reliability. Because the
   governed cell artifact contains no replicate disagreement, empirical assay
   reliability is `NA`; no reliability-improvement claim is permitted.
5. Continuous protein familiarity is `NA` unless an already-governed,
   label-blind sequence-similarity artifact exists. CD-HIT40 separation alone
   establishes component closure but is not substituted for a continuous
   familiarity score.

Failure of items 1-2 stops score fitting. Missing optional covariates do not
authorize proxies: they narrow the claim and must be reported as `NA`.

## Offline score and cross-fitting

Use the frozen final `uniform_clean` V1 development checkpoint for each of
three registered seeds (`20260831`, `20260832`, `20260833`). For each
`k in {1,2,3,5}`, draw deterministic source candidate episodes and disjoint
meta-validation episodes. The three original ATS statistics remain separated:

1. difficulty: `log1p(clean query MSE)`;
2. local transfer direction: support/query gradient cosine;
3. progress: fixed to the checkpoint endpoint and therefore excluded from an
   offline ranker when it has zero variance.

Task size, label spread, scaffold familiarity, and provenance density are
audited as nuisance/provenance columns. They are not renamed reliability.

The candidate scheduler is a small ridge scorer fitted out-of-fold over whole
source protein components to predict first-order utility to held-out
`meta_val` component gradients. Regularization is selected inside each
training fold only. The matched null has exactly the same design, parameter
count, folds, regularization grid, and compute, but deterministically permutes
the informative statistic rows within task-size strata. Predictions for a
source component are produced only by a scorer that did not fit that
component.

## Frozen destructive controls

Recompute statistics and the already-fitted out-of-fold score on the identical
episodes under:

- correct protein and correct support;
- a precomputed protein-shuffle feature from a different CD-HIT40 component;
- a support episode donated by a different CD-HIT40 component;
- deterministic label permutation across candidate episodes (including query
  labels; this keeps `k=1` falsifiable);
- ligand-only population prediction, with no pair coordinate;
- an intercept-only scalar predictor.

No scorer is refit on a destruction. No destructive-control result is used as
a feature, target, hyperparameter selector, or stopping signal.

## Frozen success criteria

Report target/component-macro estimates over all three seeds and a 9,999-draw
whole-component bootstrap. The scheduler axis passes only if every condition
below passes at every `k`:

1. **Matched-null superiority:** out-of-fold Spearman correlation with clean
   first-order utility exceeds the equal-capacity permuted-statistics null by
   at least `0.05`, and the paired 95% bootstrap lower confidence bound is
   greater than zero.
2. **Absolute utility tracking:** clean out-of-fold Spearman is at least
   `0.20`, with a 95% lower confidence bound greater than zero.
3. **Partner/support necessity:** clean first-order utility exceeds each of
   protein-shuffle and wrong-support utility, with paired 95% lower confidence
   bounds greater than zero. The fitted score must also lose at least `0.05`
   Spearman correlation with clean utility under each destruction.
4. **Label necessity:** the label-permuted score loses at least `0.05`
   Spearman correlation with clean utility, with paired 95% lower confidence
   bound greater than zero.
5. **Shortcut controls:** clean score-utility correlation exceeds ligand-only
   and intercept-only correlations by at least `0.05`, with paired 95% lower
   confidence bounds greater than zero.
6. **Nuisance survival:** matched-null superiority remains positive after both
   score and utility are residualized against every available nuisance column
   (task size, label spread, scaffold familiarity, provenance density), with a
   paired 95% lower confidence bound greater than zero.
7. All three seeds have the same favorable direction for criteria 1-6.

Ordinary MSE is not a Gate. Any all-k failure yields `REJECT` for scheduler
integration and forbids short end-to-end scheduler training. Thresholds may
not be relaxed after seeing results.

## Authorization after Gate 0/1

- Gate 0 failure: stop before fitting or scoring.
- Gate 0 pass but any Gate 1 criterion failure: retain the audit code and
  report; do not modify `train_main_v1.py` scheduling behavior and do not run
  a short training Gate.
- Gate 1 all-k pass: a later, separately recorded change may integrate the
  cross-fitted scorer into source-task selection and run a minimum-cost,
  three-seed CUDA development Gate against uniform and the matched null.
- No outcome in this lane can admit a biological representation or reopen
  R1/V1 exact-pair integration. A scheduler gain is training regularization
  unless an independent fresh biological Gate passes.

## Non-threshold implementation clarification before the final hash-bound run

A preliminary engineering run exposed two static contract ambiguities. That
run is void as scientific evidence. Before the final run, with every Gate
threshold above unchanged, the implementation is fixed as follows:

- the fitted scorer has exactly two inputs: difficulty and support/query
  gradient alignment; task size, label spread, scaffold familiarity and
  provenance density are nuisance-audit columns only and never enter scorer
  fitting;
- empty Murcko scaffold strings are missing values, never a shared scaffold;
  they are excluded from overlap calculations. Tasks with no nonempty
  scaffold use mean imputation only inside nuisance residualization and carry
  an explicit missing-fraction column.

The equal-capacity null therefore also has exactly two inputs. The final
result must bind this amended preregistration and both runner source files by
SHA256.
