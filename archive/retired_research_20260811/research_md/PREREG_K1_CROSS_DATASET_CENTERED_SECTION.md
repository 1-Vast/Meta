# K1 preregistration: cross-dataset centered-section migration Gate

Frozen before any K1 affinity metric is computed: 2026-08-11.

## Question

Does separating a support-derived task intercept from a centered, positive-ridge
ligand section improve unseen-target 5-shot regression across multiple public
DTA datasets, rather than only on MetaSieve's BindingDB corpus?

This is an engineering migration Gate. It consumes the public AdaMBind
BindingDB, Davis and KIBA CSVs for this architecture decision and is not an
independent scientific confirmation of biological specificity.

## Frozen datasets and task definition

- `bindingdb-full-data.csv`
- `davis-full-data.csv`
- `kiba-full-data.csv`
- one exact `target_sequence` is one task;
- rows require finite affinity and a valid RDKit molecule;
- duplicate canonical `(target_sequence, canonical_smiles)` observations are
  aggregated by median;
- eligible evaluation tasks require at least 25 unique ligands.

Target split is deterministic and outcome-blind. Let
`u = int(sha256(dataset + "|" + sequence)[:8],16) / 2**32`:

- source: `u < 0.70`;
- validation: `0.70 <= u < 0.80` (census only; no tuning in K1);
- evaluation: `u >= 0.80`.

The source affinity mean is the sole population predictor. Evaluation labels do
not enter it.

## Frozen representation and episodes

- RDKit Morgan bit fingerprint, radius 2, 256 bits, canonical SMILES;
- binary float64 features, divided by their L2 norm per ligand;
- `k=5` support;
- 10 deterministic support draws per eligible target;
- per draw, shuffle ligands with SHA-derived seed, use first 5 as support and
  at most the next 50 as query;
- the same support/query rows are paired across all arms;
- ridge `lambda=1.0`, no dataset-specific tuning.

## Arms

1. `population`: source-target global affinity mean.
2. `intercept`: population plus mean support residual.
3. `uncentered`: population plus the previous uncentered dual linear ridge.
4. `centered`: population plus explicit support intercept plus centered dual
   linear ridge from `CenteredKernelSection`.

The K0 primal/dual equality is a code precondition, not an empirical arm.

## Metrics and inference

- squared loss per query;
- first average within draw, then target; target is the statistical unit;
- dataset target-macro MSE for each arm;
- primary paired reductions:
  - `delta_uncentered = MSE_uncentered - MSE_centered`;
  - `delta_intercept = MSE_intercept - MSE_centered`;
- 2,000 target bootstrap draws per dataset, one-sided 95% lower confidence
  bound (5th percentile), seed 20260811;
- cross-dataset pooled statistic first standardizes each target reduction by
  the source affinity variance of its dataset, then gives each dataset equal
  weight and bootstraps targets within each dataset.

Rows, draws and queries are never treated as independent inference units.

## Migration Gate

`CENTERED_SECTION_CROSS_DATASET_PASS` requires all of:

1. K0/unit/numerical tests pass before labels are evaluated;
2. both paired point reductions are positive in every dataset;
3. for `delta_uncentered`, at least two of three dataset LCBs are positive;
4. for `delta_intercept`, at least two of three dataset LCBs are positive;
5. both equal-dataset pooled LCBs are positive;
6. no dataset has non-finite predictions or more than 1% rejected eligible
   episodes from numerical failure.

If any condition fails, the module remains under `research/` and is not copied
to `model/` or `scripts/`. Thresholds, fingerprint size, ridge and datasets are
not changed after seeing results.

## Scope limits

A PASS supports a dataset-agnostic centered support solver interface. It does
not admit T-BASIS, PAREF, protein specificity, a rich kernel, calibrated
uncertainty, CSMO `z`, or a production DTA model.

