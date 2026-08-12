# Preregistration: BindingDB development interaction existence

Stage: `E-AFF-CQ-R1_BINDINGDB_DEVELOPMENT_INTERACTION`

Status: frozen after CQ-R0 metadata census and before numeric BindingDB values
are parsed by the new CQ pipeline.

## Inputs

- BindingDB-curated Articles release 202608 only;
- the sealed CQ-R0 metadata projection with SHA-256
  `8ac11beb8e058f373b03e7f392c8e993b49d64010012d887aaa43d312ef2ce00`;
- the primary `document + endpoint + normalized protocol` panel definition;
- single-chain targets with an explicit sequence and stereo-aware ligand key.

ChEMBL-derived, patent, PDSP, PubChem, DAVIS, recipient and mutation values are
not read. Ki and Kd are separate strata and are never pooled numerically.

## Label contract

Only finite positive uncensored numeric Ki/Kd values in nM are accepted. The
canonical direction is:

```text
pK = 9 - log10(value_nM)
higher pK = stronger affinity
```

Rows with `<`, `>`, `~`, ranges or non-numeric text are excluded and counted.
Repeated measurements of the same `(panel, protein, ligand, endpoint)` cell are
averaged before quotient construction. Replicate sample variance and variance
of the cell mean are reported when at least two measurements exist.

## Quotient statistic

Within each connected panel graph, fit the declared additive nuisance

```text
y_e = mu + alpha_protein[e] + beta_ligand[e] + residual_e
```

by deterministic float64 least squares. The retained rank is
`E - rank(X_full)`; for the protein/ligand incidence design this equals the
cycle rank. Report the rank-normalized residual mean square and numerical
orthogonality `max(abs(X.T @ residual))`.

Panel statistics are averaged within source document and then document-macro.
Bootstrap units are documents, never cells, edges or quotient coordinates.

## Development signal criterion

Ki is primary and Kd is a separately reported replication stratum. A stratum is
development-trainable only if all conditions hold:

```text
cycle_positive_panels >= 50
total_retained_rank >= 1000
document_units >= 30
construct_mapping_coverage >= 0.95
quotient_RMS >= 0.25 pK units
95% document-bootstrap LCB(quotient_RMS) > 0
max_additive_orthogonality_error <= 1e-7
```

The 0.25 pK floor is fixed before value parsing and denotes a modest
approximately 1.8-fold affinity ratio. Replicate noise is diagnostic in R1;
noise-corrected population claims require a separate covariance/power stage.

## Verdicts

Exactly one:

```text
CQ_R1_LABEL_OR_NUMERICAL_CONTRACT_FAIL_CLOSED
CQ_R1_DEVELOPMENT_INTERACTION_NOT_OBSERVED
CQ_R1_DEVELOPMENT_INTERACTION_OBSERVED
```

The last verdict authorizes construction and smoke-testing of one research-only
linear response on the frozen 288D T-BASIS. It does not authorize a biological
claim, model/ migration, affinity `z`, few-shot adaptation, DAVIS, recipient or
production use.
