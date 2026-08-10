# Preregistration: BindingDB cycle-quotient source census

Stage: `E-AFF-CQ-R0_BINDINGDB_PANEL_IDENTIFIABILITY`

Status: registered before BindingDB acquisition or parsing. Audit-only.

## Frozen source order

Primary source is the release-pinned BindingDB curated-articles TSV plus assay
mapping and target FASTA. ChEMBL-, PubChem-, PDSP- and CSAR-derived BindingDB
rows are excluded. Patent data are a separate future stratum and cannot rescue
the article census after its result is known.

Before acquisition, record URL, release date, published checksum, local SHA-256,
licence/terms and redistribution decision. Failure to verify these fields is
`CQ_R0_SOURCE_OR_LICENCE_FAIL_CLOSED`.

## Label-blind fields

R0 may parse only source/entry/assay identifiers, DOI/PMID/patent, publication
and curation dates, target sequence and construct annotation, organism, ligand
connectivity/stereochemistry identifiers, endpoint availability/relation flags,
and assay description. Numeric affinity values are not deserialized.

## Panel candidates

Materialize all three definitions before any value read:

1. strict BindingDB EntryID-AssayID;
2. DOI plus normalized experimental-method signature plus endpoint;
3. deposited matrix identifier plus endpoint.

No definition may be selected by its eventual affinity statistic. Definitions
that mix organisms, incompatible constructs, endpoints, units or measurement
types fail closed.

## Dependency closure

Union dependencies by source document/patent, exact measurement reaction,
target sequence homology, ligand connectivity and Bemis-Murcko scaffold.
ChEMBL overlap is closed by source flag, DOI and reaction identity. Report the
strict union closure and a preregistered DataSAIL-style two-dimensional
sensitivity; neither may be chosen after values are read.

## Census outputs

For every panel report `E,P,L,c,d_cycle=E-P-L+c`, censoring/endpoint flags and
dependency component. Also report:

- total raw and covariance-conservative quotient rank;
- number of cycle-positive panels and dependency components;
- largest component share;
- target-family and scaffold diversity;
- publication-time eligible confirmation components;
- exact construct mapping coverage;
- ChEMBL overlap exclusions.

## Gate

```text
components >= 60
largest_component_share <= 0.25
conservative_effective_quotient_rank >= 245
construct_mapping_coverage >= 0.95
ChEMBL_overlap_in_primary == 0
```

Exactly one verdict:

```text
CQ_R0_SOURCE_OR_LICENCE_FAIL_CLOSED
CQ_R0_PANEL_DATA_NOT_IDENTIFIABLE
CQ_R0_PANEL_SOURCE_IDENTIFIABLE
```

Only the last verdict authorizes a separate CQ-R1 interaction-existence
preregistration. It does not authorize affinity model training, GPU use,
mutation data, DAVIS, few-shot adaptation or biological `z`.
