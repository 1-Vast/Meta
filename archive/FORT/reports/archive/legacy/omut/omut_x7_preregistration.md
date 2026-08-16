# OpenMut `OMUT-X7` preregistration

**Frozen before execution:** 2026-07-28.
**Predecessors:** `OMUT-X6`,
`OMUT_X6_REAGENT_VERIFICATION_INSUFFICIENT_STOP`; `OMUT-X0`, BindingDB
sequence-exact upper bound with no source-native assay ID.
**Stage boundary:** no-outcome BindingDB composite-context audit.

## 1. Question

X0 recovered 37 sequence-exact BindingDB `Ki/Kd, k>=4` mutation components
but rejected them because the frozen schema has no assay identifier. X7 asks:

> Does exact equality over independent source and condition fields provide a
> sufficiently complete and diverse conservative composite assay context?

This is not permission to invent an assay ID. It is a falsifiable audit of
whether the archive already carries a lossless-enough row context.

## 2. Frozen row relation

The X0 exact construct relation is unchanged:

- one verified single substitution;
- exact mutant target sequence;
- WT target sequence obtained by exact reversion of that substitution;
- same ligand parent InChIKey and endpoint;
- at least four distinct shared ligands.

For each WT/mutant row pair, X7 additionally requires the same nonempty:

1. DOI/PMID provenance unit;
2. `Curation/DataSource`;
3. `Authors`;
4. `Institution`;
5. pH;
6. temperature.

pH and temperature strings are normalized only for whitespace and case; no
rounding or tolerance is allowed. Their values are used only to hash and
compare signatures and are never retained.

A ligand is strict-matched when at least one WT/mutant row pair has an
identical complete signature. A component is strict only with four such
ligands.

## 3. Sensitivity upper bound

One non-primary sensitivity signature is reported:

- exact document and curation source;
- exact nonempty authors or exact nonempty institution;
- pH/temperature must agree whenever either side is populated, but both may
  be missing.

This missing-condition upper bound cannot unlock training.

## 4. Combined registry and independence

The primary registry is the union of:

- X5's five strict ChEMBL components;
- X7 strict BindingDB components.

Duplicate `(accession, mutation, endpoint)` components collapse to one.
Independent provenance units are unique `(accession, document, endpoint)`
clusters; multiple mutations or ligand pairs in one cluster do not increase
independent `n`.

## 5. Gates

| gate | condition |
| --- | --- |
| `X7_X0_X6_BOUND` | X0/X6 result hashes, X0 exact-sequence firewall, and X6 stop bind |
| `X7_ARCHIVE_BOUND` | BindingDB archive checksum and 640-column header match X0 |
| `X7_NO_OUTCOME_MATERIALIZED` | endpoint values are presence-tested only; no numeric outcome is retained |
| `X7_EXACT_CONSTRUCT_PRESERVED` | every candidate uses the exact X0 sequence-reversion construct relation |
| `X7_COMPLETE_SIGNATURE_ENFORCED` | every strict pair has all six nonempty equal context fields |
| `X7_DISCOVERY_ACCESSION_EXCLUDED` | `P15056` contributes zero primary components |
| `X7_FAMILY_MAPPED` | every primary component has a nonempty family disposition |
| `X7_COMPONENT_TOPOLOGY_ADEQUATE` | at least 25 unique components, six broad families, and largest accession share at most 50% |
| `X7_PROVENANCE_TOPOLOGY_ADEQUATE` | at least 12 independent accession-document-endpoint clusters and no cluster supplies more than 25% of components |

Verdicts:

- all gates pass: `OMUT_X7_COMPOSITE_CONTEXT_REGISTRY_ADEQUATE`;
- execution/firewall gates pass but either topology fails:
  `OMUT_X7_COMPOSITE_CONTEXT_REGISTRY_INADEQUATE_STOP`;
- any binding, checksum, construct, signature, or firewall gate fails:
  `OMUT_X7_INCOMPLETE_STOP`.

## 6. Deliverables

- `research/omut_x7.py`
- `tests/test_omut_x7.py`
- `reports/active/omut_x7.json`
- `reports/active/omut_x7_decision.md`
