# RBSDD public-source recovery round decision

Date: 2026-07-29  
Decision: `STOP_RBSDD_PUBLIC_SOURCE_ROUND_NO_IDENTIFIABLE_INDEPENDENT_TEACHER`

## Binding requirement

RBSDD requires one admissible Stage-1 population that jointly provides:

1. exact target and ligand identity;
2. real ligand-conditioned complex/contact semantics;
3. repeated same-target directed ligand edits or paired complexes;
4. exact affinity with row-level source provenance and independent support.

These properties cannot be supplied by summing marginally compatible sources
unless their rows can be joined exactly and their evidence dependence can be
measured.

## Audited sources

| Source | What was recovered | Decisive deficit |
| --- | --- | --- |
| PBCNet2.0-D0 | Open CC BY 4.0 release, pair identities, pose directories, separable labels | No original BindingDB row/assay/document/source lineage; repository has no license file despite README claim |
| BioLiP2-D1C | 66,660 affinity-blind exact-complex rows, 391,262 same-target ligand pairs, ChEMBL target/ligand support `70.30%/49.79%`, PubMed concentration `0.8161%` | After generic-ligand removal and provenance collapse, the largest bipartite 2-core component still contains `79.9767%` of edges |
| PSICHIC-G0 | Apache-2.0 code, pinned weights, documented identity and task-label schema | Public XL inventory has no visible `train.csv`; documented rows expose no source dataset or stable source-record ID |

## Decision

No audited source satisfies all four requirements. RBSDD model construction,
affinity loading, and weight inference remain unauthorized. BioLiP2 is
retained only as an exact-complex/contact index for future cross-validation;
it is not counted as independent Stage-1 evidence. PBCNet2.0 and PSICHIC
remain metadata-only candidates with explicit reopening conditions in their
decisions.

This closes the current public-source recovery round. It does not prove that
the required information object cannot exist, and it does not justify another
attention, router, alignment, or representation module. A new route must
first establish a nonzero identifiable target-domain information object under
destruction and source-independence controls.

