# PCIC-O0 execution clarification

Date: 2026-07-29
Status: frozen after synthetic parser tests and before any real O0 projection,
topology metric, or operator metric

This clarification resolves implementation ambiguities in
`pcic_o0_preregistration_2026-07-29.md`. It strengthens the firewall and makes
the O0-I tests non-circular. It does not relax an O0-P or O0-I gate.

## Mixed raw-file identity projection

The two frozen gzip inputs contain all split identities in flat JSON rows.
It is impossible to know whether a row belongs to a TRAIN target without
first reading its `target_chembl_id`.

The allowed operation is therefore:

1. scan every raw row as bytes;
2. decode only the exact fixed metadata whitelist declared in
   `research/pcic_o0.py`;
3. skip `value_nM` and `pK` as byte spans without numeric or textual
   conversion;
4. retain and emit only rows whose target occurs in the Parquet
   `dual_cold_split == train` column projection;
5. discard all other safe identities immediately.

The projector rejects unknown fields, duplicate fields, nested values, and
non-flat rows. Exceptions and logs must not include raw records or value
spans. The full historical `json.loads` loader in
`research/p0_identifiability.py` remains prohibited.

Registry `affinity` and label-derived `replicate_sd` are both protected. O0
may load neither column.

## Version and digest binding

The official ChEMBL status endpoint reported:

```text
chembl_db_version = ChEMBL_37
chembl_release_date = 2026-05-01
status = UP
```

before real projection. A source contract must freeze path, byte size, and
SHA-256 for:

- both raw gzip files;
- the registry Parquet;
- the API source manifest;
- both frozen feature archives;
- the preregistration;
- this clarification;
- the exact implementation.

Preparation verifies the contract before reading, hashes each raw file again
after streaming, and writes only a safe-cell artifact and a release-37
document-metadata cache. A second run contract freezes those two generated
artifacts before topology gates are evaluated.

The ChEMBL document cache is admissible only if a fresh status response still
reports `ChEMBL_37`. It may contain document ID, DOI, PMID, exact patent ID,
document type, source ID, and release metadata. It may not contain activity
outcomes. Unknown lineage remains missing; a document ID alone is not treated
as proof of independent provenance.

## Coordinate definition for O0-I

If and only if O0-P passes:

- fit target PCA on unique TRAIN target entities, not edge-weighted rows;
- fit ligand PCA on unique TRAIN ligand connectivities, not edge-weighted
  rows;
- center and whiten both coordinates;
- retain exactly eight target and eight ligand coordinates;
- order `x_(t,l) = v_l tensor u_t`;
- stop for coordinate instability if the relative eigengap between the
  eighth and ninth component is at most `1e-6`;
- use `allow_pickle=False`, record exact NPZ keys/dtypes, and bind row maps by
  digest.

Whitening is part of the primary parameterization. Unwhitened coordinates
are not a rescue sensitivity.

## CUDA nuisance projection

The primary projector operates in CUDA float64. For unweighted O0 it cycles
exact group-mean subtractions over target, ligand, assay, document, and any
certified nonnested context.

It converges only when:

- maximum relative KKT error is at most `1e-10`;
- relative full-sweep change is at most `1e-12`;
- repeat-projection idempotence satisfies the original `1e-7` gate.

Maximum sweeps are `20,000`. Failure to converge is
`PCIC_O0_NUMERICAL_NO_DECISION`, not a biological failure. Ridge, jitter,
damping, eigenvalue clipping, and singular-value flooring remain forbidden.

Rank is obtained from a direct float64 tall QR/SVD of the residual operator,
not from thresholding float32 Gram eigenvalues. The original
`sigma_j / sigma_1 >= 1e-8` rule applies to singular values of `A`.

Independent CPU validation uses four individual columns and four
seed-1729 Rademacher combinations with undamped LSMR and LSQR. Full
`64 x 64` CPU/GPU Gram comparison checks arithmetic on the already projected
matrix; it does not imply 64 independent CPU sparse solves.

## Non-circular replication and query tests

Joint homology-scaffold-lineage components, not rows or circuits, are
assigned to five deterministic folds.

For every fold `f`:

- `A_f`, built from that fold alone, measures source replication;
- `A_-f`, built from the other four folds, defines the estimable row space
  for held queries;
- held target-ligand pairs are counted once regardless of the number of
  assays/documents.

The within-fold replication subspace is the leading eight right-singular
directions of `A_f`. All five folds must have rank at least eight. For a held
fold, form an eight-dimensional consensus from the normalized leading-eight
projectors of the other four folds. The frozen replication statistic is:

```text
T_rep = minimum over held folds of the smallest canonical correlation
        between the held leading-eight subspace and the four-fold consensus.
```

The original `0.80` threshold applies to `T_rep`. This is not the exact
intersection, so the correlation is not identically one.

Query residual `r_q` is computed against the admitted right-singular basis of
`A_-f`:

```text
r_q = sqrt(max(0, ||x_q||^2 - ||V_-f^T x_q||^2)) / ||x_q||.
```

Zero-norm queries are unsupported rather than assigned residual zero. Full
data may not define a held query's span.

## Permutation correction

Generic 64-dimensional row rank and row-span coverage can be algebraically
saturated for semantic and entity-permuted coordinates. Therefore query
coverage is retained as an absolute estimability gate but is not used as a
semantic permutation statistic.

The two semantic negative-control statistics are:

- `T_rep` above;
- joint-component information ESS from row-energy concentration.

Target and ligand mappings are independently permuted 20 times within
endpoint, fold, and frozen log2 degree bin. PCA loadings, topology, row
multiplicity, and folds are not refit. Every map digest is saved.

For either family and either statistic:

```text
p = (1 + number of permuted statistics >= semantic statistic) / 21.
```

A semantic pass requires strict superiority to all 20 controls (`p = 1/21`)
for both statistics and both permutation families. Interpolated percentiles
are prohibited.

If row rank or query span is saturated across all controls, report
`CONTROL_METRIC_ALGEBRAICALLY_SATURATED`; do not interpret saturation as a
biological semantic failure.

## Conditional T1 dimension

An O0 pass does not authorize all 64 ambient coordinates. A later T1 may fit
only the frozen eight-dimensional cross-fold consensus subspace that passed
O0. The other 56 ambient directions remain unidentified and must not be
restored by ridge or a factorized model.

