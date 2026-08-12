# Amendment 01: executable open-data quotient contract

Stage: `E-AFF-CQ-R0_BINDINGDB_PANEL_IDENTIFIABILITY`

Status: frozen before the new BindingDB 202608 census. This amendment corrects
four implementation-level defects in the original preregistration. It does not
change any measured result because no CQ-R0 metric has yet been computed.

## Source chronology and firewall

The BindingDB 202608 Articles and Assays archives were acquired on 2026-08-08,
before the CQ-R0 preregistration. Their hashes and acquisition manifest already
exist locally. The preregistration is therefore prior to the **new projection,
census and adjudication**, not prior to acquisition.

The article TSV is a monolithic 640-column file. A trusted extractor necessarily
traverses rows containing affinity text. It must never deserialize, log, return
or persist numeric affinity values during R0. The machine report must distinguish:

```text
affinity_bytes_traversed_by_trusted_extractor = true
numeric_affinity_values_parsed = 0
numeric_affinity_values_exposed = 0
numeric_affinity_values_used = 0
```

Only the sealed metadata projection may be consumed by the census.

## Correct panel keys

`BindingDB Reactant_set_id` joins `ENTRYID` in the assay table. The pair
`(ENTRYID, ASSAYID)` identifies an edge-to-assay mapping and is not a panel;
using it as a panel key would force zero cycle rank.

Materialize these label-blind candidates:

1. `DOI/PMID + endpoint + normalized assay protocol family` (primary);
2. `DOI/PMID + endpoint` (broader sensitivity);
3. exact `ENTRYID + ASSAYID + endpoint` (edge-level negative control).

The primary protocol family is a deterministic hash of normalized assay name
and description tokens. It is not selected using affinity values. Because
target-specific assays may occur inside one publication matrix, any positive
result supports only panel-compatible selectivity until replicated by an
independent publication or protocol.

## Canonical quotient

For a panel with response vector `y`, prediction `q`, declared covariance
`Sigma`, and full nuisance design `X_full`, use:

```text
y_star = Sigma^(-1/2) y
q_star = Sigma^(-1/2) q
X_star = Sigma^(-1/2) X_full
P_perp = I - X_star pinv(X_star)
loss_panel = ||P_perp (y_star - q_star)||^2 / retained_rank
retained_rank = E - rank(X_star)
```

This avoids a basis-dependent cycle truncation. When only diagonal measurement
variances are available, the claim is limited to heteroscedastic whitening.
Panel losses are averaged within dependency component and then component-macro;
raw quotient coordinates are never treated as IID.

## Two distinct authorization gates

### Development training readiness

This gate may authorize research-only optimization on open source data:

- at least one release-pinned source has a non-zero retained quotient/ranking
  design after endpoint and panel separation;
- source, panel, target sequence/construct and ligand connectivity are traceable;
- train/dev partitions enforce document, protein-homology and scaffold conflicts,
  with all discarded conflicts reported;
- training and evaluation masks are frozen and disjoint;
- modalities are not numerically mixed: Ki, Kd, Kdapp and single/two-dose
  profiling have separate loss strata.

This gate does not require 60 independent confirmation components and cannot
support a population-level biological claim.

### Biological claim/admission readiness

The original thresholds remain unchanged:

```text
components >= 60
largest_component_share <= 0.25
construct_mapping_coverage >= 0.95
ChEMBL_overlap_in_primary == 0
```

The old `effective_quotient_rank >= 245` proxy is not itself an independence or
power certificate. It is retained as a descriptive continuity statistic. A
claim additionally requires component-level null simulation with type-I error
at most 0.05 and power at least 0.80 at a preregistered biological effect.

## R0 verdicts

R0 emits both booleans and one source verdict:

```text
development_training_ready: true | false
biological_claim_ready: true | false

CQ_R0_SOURCE_OR_LICENCE_FAIL_CLOSED
CQ_R0_PANEL_DATA_NOT_IDENTIFIABLE
CQ_R0_DEVELOPMENT_SOURCE_IDENTIFIABLE
CQ_R0_CLAIM_SOURCE_IDENTIFIABLE
```

Only a source-identifiable verdict authorizes a separate numeric-label stage.
No affinity value, GPU training, model change, few-shot section or biological
`z` is authorized by this amendment itself.
