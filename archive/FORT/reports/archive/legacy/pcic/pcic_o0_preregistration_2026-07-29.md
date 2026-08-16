# PCIC-O0 label-blind operator-identifiability preregistration

Date: 2026-07-29
Status: frozen before implementation and before any new label access

## Purpose

`PCIC-O0` tests the only increment not answered by `P0-Cycle-A`:

> Does the existing ChEMBL TRAIN nuisance-null edge space contain a stable,
> provenance-replicated, strict-dual-cold-estimable direction in one frozen
> target-ligand operator coordinate?

It does not test whether the null space is nonempty, does not enumerate
circuits, and does not train an affinity predictor.

## Protected boundary

O0 may read identifiers and covariates from ChEMBL TRAIN only. It may not
decode, deserialize, aggregate, compare, or log:

- `pK`;
- `value_nM`;
- registry `affinity`;
- development or confirmation outcome fields;
- sealed outcomes.

The safe raw projection must reject protected keys at every nesting depth or
skip their bytes without conversion. The existing `load_train_source`
function is prohibited because it explicitly decodes `record["pK"]`.

Development, confirmation, and sealed query identities may be used only
after a separate firewall review. The primary O0 query-support audit uses
deterministically held TRAIN homology-scaffold-lineage components.

## Frozen inputs

- ChEMBL source version: 37, bound to
  `dataset/public/chembl_37/processed/api_source_manifest.json`.
- Raw endpoint files:
  `chembl37_pKi.jsonl.gz` and `chembl37_pKd.jsonl.gz`.
- Split registry:
  `dataset/public/chembl_37/processed/dualcold/registry.parquet`.
- Protein features:
  `dataset/public/chembl_37/processed/dualcold/target_esm2.npz`.
- Ligand features:
  `dataset/public/chembl_37/processed/dualcold/ligand_features.npz`.
- Seed: `1729`.
- Endpoints: pKi and pKd analyzed and decided separately.

Input sizes and SHA-256 values must be written to the result before any
gate is evaluated.

## O0-P: metadata and topology phase

Construct one row per exact
`(endpoint, target, ligand-connectivity, assay, document)` cell using only
TRAIN targets and non-outcome fields.

The row projection must provide:

- target ChEMBL ID and accession/homology mapping;
- canonical nonisomeric ligand connectivity and Murcko scaffold;
- assay ID;
- document ID;
- exact ChEMBL-37-bound publication/patent lineage closure where available;
- a declared nonnested context inventory.

DOI, PMID, patent family, institution, and site are not assumed present.
Each must be classified as exact, missing, nested in document, or
non-identifying. No unversioned live API result may enter the primary
registry.

Build dependency components by transitive closure over target homology,
ligand scaffold, and exact provenance lineage. Report raw components and the
2-core, largest-component share, endpoint coverage, and missing-lineage
share.

O0-P stops with `STOP_PCIC_O0_PROVENANCE_OR_TOPOLOGY_INADEQUATE` if any of
the following holds:

1. protected outcome keys are decoded or emitted;
2. an input is unversioned or lacks a recorded hash;
3. exact lineage closure is unavailable for more than 5% of retained cells;
4. fewer than 88 joint dependency components remain;
5. the largest joint dependency component exceeds 20% of cells;
6. fewer than five nonempty lineage-disjoint folds can be constructed for
   either endpoint.

An endpoint may stop independently. A pKi pass cannot rescue pKd.

## Frozen operator coordinate

Fit all coordinate transforms without outcomes and on TRAIN entities only.

- Protein: center frozen pooled ESM-2 features and retain the first eight
  deterministic SVD coordinates.
- Ligand: center the existing frozen ligand features and retain the first
  eight deterministic SVD coordinates.
- Fix component signs by requiring the largest-absolute loading to be
  positive.
- Form `x_(t,l) = v_l tensor u_t`, giving a direct `8 x 8`, 64-parameter
  operator coordinate.

No alternative protein model, ligand descriptor, coordinate dimension, or
learned encoder may replace the primary coordinate after seeing O0 results.
A `16 x 16` coordinate may be reported as confirmatory sensitivity only and
cannot change the primary decision.

## Nuisance projection

The primary unweighted design is:

```text
Z = [target, ligand, assay, document, nonnested exact context].
```

Homology, scaffold, and lineage are used for blocked replication and
effective-sample accounting. A metadata field enters `Z` only if it is
demonstrably nonnested in the existing columns.

For `D = W^(1/2)`, the weighted sensitivity analysis uses only pre-outcome
weights derived from record count or prespecified measurement metadata:

```text
Z_tilde = D Z
X_tilde = D X
A = (I - Z_tilde Z_tilde^+) X_tilde
M = A^T A.
```

The unweighted result is primary. A weighted result cannot rescue it.
No ridge, diagonal jitter, or singular-value floor may enter the rank or
query-span calculation.

The projector is implemented implicitly. Explicit minimal-cycle enumeration
and pseudo-outcome row expansion are prohibited.

## O0-I: operator-information phase

Run on CUDA in `D:\anaconda\envs\drug`:

- multiway residualization of all 64 operator columns;
- `M = A^T A`;
- singular/eigen decomposition;
- five leave-lineage-fold-out information matrices;
- common-subspace and query-span calculations;
- 20 target-feature and 20 ligand-feature entity permutations.

CPU work is restricted to parsing, RDKit canonicalization, metadata/graph
construction, and independent sparse LSMR/LSQR/KKT checks.

Report:

1. numerical rank using `sigma_j / sigma_1 >= 1e-8`;
2. condition number within the admitted nonzero subspace;
3. fold-wise rank and minimum canonical correlation of the common subspace;
4. information effective sample size
   `(sum trace_g)^2 / sum(trace_g^2)` by lineage and joint component;
5. maximum information-trace fraction for any lineage, homology family, and
   scaffold family;
6. the fraction of held strict-dual-cold TRAIN queries with
   `r_q <= 0.01`;
7. the same stability and query metrics for all entity-permuted controls;
8. KKT, idempotence, and CPU/GPU Gram disagreement.

## Frozen O0-I gates

An endpoint passes only if all conditions hold:

1. primary numerical rank is nonzero;
2. admitted-subspace condition number is at most `1e6`;
3. a nonzero common direction exists in all five lineage-disjoint folds;
4. the minimum fold-to-common canonical correlation is at least `0.80`;
5. joint-component information ESS is at least `88`;
6. no single lineage carries more than `5%` of total information trace;
7. no homology or scaffold family carries more than `20%`;
8. at least `95%` of held strict-dual-cold TRAIN queries have
   `r_q <= 0.01`;
9. semantic common-subspace stability and query coverage each exceed the
   95th percentile of both matched entity-permutation controls;
10. maximum relative KKT error is at most `1e-8`, idempotence error at most
    `1e-7`, and CPU/GPU Gram disagreement at most `1e-6`;
11. zero protected outcomes are decoded and zero
    development/confirmation/sealed outcomes are accessed.

Failure returns `STOP_PCIC_O0_OPERATOR_NOT_IDENTIFIABLE`. The result must
state every failed gate; no single passing endpoint authorizes a generic
affinity claim.

Pass returns `REQUEST_PCIC_T1_PREREGISTRATION`. It does not authorize outcome
training or development scoring.

## Conditional T1 requirements

A future T1 is permitted only by a superseding preregistration after O0 pass.
It must:

- use the frozen 64-dimensional direct operator;
- compare equal-budget ligand-only, additive, HCRR, rectangle-only, full
  nuisance-null, no-provenance, PB-CEC, and entity-shuffled arms;
- retain wrong-target and wrong-ligand destruction;
- use joint dependency components, not circuits, as independent units;
- preserve the existing `0.0586` macro within-target Spearman MDE gate and
  `1.02 x B0` RMSE ceiling;
- require each target/ligand destruction to remove at least 70% of the gain;
- require provenance-disjoint sign stability and held-component
  risk-coverage;
- treat HCRR as a baseline only. Its prior withdrawal remains in force until
  the superseding preregistration exists.

## Compute estimate

For pKi, `231,090 x 64` float32 operator features occupy about 59 MB and the
Gram step requires approximately `9.5e8` multiply-adds. Expected elapsed
time is dominated by safe metadata projection and blocked controls:

- O0-P: 4-8 working hours if ChEMBL-37 metadata is locally recoverable;
- O0-I primary endpoint: 4-8 GPU hours including validation;
- both endpoints, five folds, and 40 permutations: up to three working days.

If O0 fails, no T1 compute is spent. If it passes, an equal-budget
exploratory T1 and mandatory controls are estimated at four additional
working days.
