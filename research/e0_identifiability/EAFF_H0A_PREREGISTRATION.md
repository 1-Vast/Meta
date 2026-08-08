# E-AFF-H0A Task-Local Radial Headroom Diagnostic

Status: registered after E-AFF-P0 shared-direction failure and before any H0A
feature generation, fit, or metric computation.

## Question

When one population-shared 288D radial direction fails, does the same frozen
T-BASIS basis contain task-local held-out affinity headroom?

This is an oracle/headroom diagnostic. It is not a production model and cannot
establish target-specific biology because each task is one
`target x endpoint x assay x context` unit.

## Selection

Eligible tasks have at least 40 distinct ligand states. Select one task per
closure component by SHA-256 of `EAFF-H0A|task_id`. Select 40 ligand states by
SHA-256 of `EAFF-H0A-LIGAND|task_id|ligand_state_key`; the first 20 are fit and
the remaining 20 are untouched task-local test. No affinity value enters
selection.

## Frozen Inputs

The P1B checkpoint, local states, chemistry classes, six RBF centers, T-BASIS
radial calibrator, D1 folds and closure assignments are identical to E-AFF-P0.
The closure-OOF ligand prior is refit by the same fixed `alpha=10` procedure on
the full governed source corpus.

## Head

For each task independently, fit one deterministic standardized Ridge direction
with `alpha=10` on all residual differences among the 20 fit ligands. No point
loss, deranged pair, coupling null, task ID embedding or test value enters the
fit. Apply the frozen task direction to 20 held-out ligands.

## Arms And Inference

`L`, `C`, `D` and `N` have the same meanings as E-AFF-P0. Each task direction is
shared across its correct, wrong-protein and coupling-null test arms. Report
task CI and component-macro contrasts; bootstrap closure components.

Headroom requires `C-L >= 0.03` with positive 95% LCB. Partner-specific
headroom additionally requires `C-D >= 0.03` with positive LCB. `C-N` is an
attribution diagnostic and must be positive with positive LCB for a coupling
claim.

## Verdicts

- `TASK_LOCAL_RADIAL_HEADROOM_AND_PARTNER_SPECIFICITY_OBSERVED`
- `TASK_LOCAL_RADIAL_HEADROOM_WITHOUT_PARTNER_SPECIFICITY`
- `TASK_LOCAL_RADIAL_HEADROOM_NOT_OBSERVED`
- `H0A_DATA_OR_MAPPING_FAIL_CLOSED`

Only the first verdict may authorize a separately registered H0-B cross-assay
target-transport diagnostic. It still does not admit a statistic to `z`.
