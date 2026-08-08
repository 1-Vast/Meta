# E-AFF-P0 Fixed-Radial Source Feasibility

Status: registered before any E-AFF model fit or metric computation.

## Question

Does the frozen, structure-calibrated T-BASIS `8 x 6 x 6` radial tensor contain
a population-shared direction that improves within-assay affinity ranking beyond
a closure-OOF ligand prior? If not, is there label-blind data support for a
separate task-local or cross-assay target-headroom study?

This is a lightweight source feasibility study, not the terminal E-AFF Gate.
It reads ChEMBL37 Ki/Kd source labels. It must not read DAVIS or recipient
labels and cannot authorize RFSA, production `z`, or P2-P4.

## Frozen Inputs

- D0/D1 ChEMBL37 immutable corpus and five closure-component folds.
- P1B epoch-5 checkpoint with SHA-256 recorded in the output manifest.
- E0 local-state cache and its manifest hashes.
- T-BASIS-R0 six-bin RBF definition and train-only `6 x 6` radial calibrator.
- Ligand chemistry channels `chemistry[32:40]` and six residue chemistry
  channels from the label-blind E0 cache.

The affinity stage cannot update P1B, local states, RBF centers, chemistry
classes, radial calibration, or tensor normalization.

## Score-Blind Panel

For each of the 245 E0-input closure components, select one task using SHA-256
of `EAFF-P0|task_id`, then select 20 distinct ligand states using SHA-256 of
`EAFF-P0-LIGAND|task_id|ligand_state_key`. Duplicate measurements for one
task-ligand state are collapsed by median pAffinity. This produces 245 tasks
and 4,900 task-ligand observations before mapping failures.

Wrong proteins are selected before scoring from the same fixed outer fold,
with a different D1 homology component and exact local identity below 0.40.
The map is score-blind and hashed. Wrong proteins are nuisance controls, not
non-binder labels.

## Ligand Prior

Fit five closure-OOF Ridge ligand priors on all eligible E0 task-ligand median
rows using only the frozen 128D pooled ligand state. `alpha=10` is fixed; each
task has total sample weight one. No mechanism feature enters this prior.

## Shared Direction

For each outer fold, fit one `w in R^288` on the other four folds. The primary
target is the residual difference

```text
(pAffinity_a - pAffinity_b) - (ligand_OOF_a - ligand_OOF_b).
```

All within-task pairs are used, but each task has total loss weight one. Fit a
deterministic standardized Ridge with fixed `alpha=10`, no intercept, no point
loss, and no deranged/null examples. This pilot deliberately isolates ranking
semantics and avoids assay-offset fitting.

## Evaluation Arms

- `L`: closure-OOF ligand prior.
- `C`: `L + w^T phi(correct protein, ligand)`.
- `D`: `L + w^T phi(wrong protein, ligand)` with the same `w`.
- `N`: `L + w^T phi_null(correct protein, ligand)` with the same `w`.

The project-defined coupling null preserves chemistry-pair and radial
marginals and deletes their coupling. It is not a physical nonbinding state.

Primary aggregation is task CI, averaged within closure component, then macro
over components. Confidence intervals resample closure components. Task/row-IID
bootstrap is prohibited.

## Feasibility Criteria

```text
C - L >= 0.03 and 95% LCB > 0
C - D >= 0.03 and 95% LCB > 0
C - N > 0 and 95% LCB > 0
```

Ki and Kd are reported separately. `pAffinity` must pass an executable audit
that stronger binding is numerically larger.

## Conditional H0 Boundary

If either `C-L` or `C-D` fails, this run performs only an H0 data census. It
does not fit a task-specific 288D head, because 20-point tasks are insufficient
and selecting its regularization after the shared result would confound the
triage. A separate H0 stage is supportable only if the census establishes:

- at least 40 distinct ligands for task-local fit/test; and
- multiple assay/document tasks for the same target for biological target
  transport.

Task-local headroom and cross-assay target transport must remain distinct.

## Terminal Pilot Verdicts

- `SHARED_RADIAL_AFFINITY_FEASIBILITY_OBSERVED`
- `SHARED_DIRECTION_NOT_OBSERVED_H0_DATA_SUPPORTED`
- `SHARED_DIRECTION_NOT_OBSERVED_H0_DATA_INSUFFICIENT`
- `PILOT_DATA_OR_MAPPING_FAIL_CLOSED`

No verdict identifies physical free energy or admits a coordinate to `z`.
