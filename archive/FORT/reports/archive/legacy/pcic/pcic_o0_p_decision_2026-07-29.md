# PCIC-O0-P decision

Date: 2026-07-29
Decision: `STOP_PCIC_O0_PROVENANCE_OR_TOPOLOGY_INADEQUATE`

## What was executed

`research/pcic_o0.py` scanned the two frozen ChEMBL-37 raw endpoint files with
a strict flat-JSON whitelist. It decoded safe identity/covariate fields only
and skipped `value_nM` and `pK` as byte spans.

The source and run contracts bind:

- raw pKi/pKd gzip files;
- TRAIN split registry;
- source manifest and official ChEMBL-37 status;
- feature archives;
- preregistration, execution clarification, and implementation;
- generated safe-cell Parquet and document metadata cache.

The official status response reported ChEMBL 37, release date 2026-05-01.
The document endpoint returned all 9,643 requested nonempty ChEMBL document
IDs.

## Firewall result

| Item | Result |
| --- | ---: |
| Physical raw rows | 484,616 |
| Safe rows decoded | 484,616 |
| Protected value spans skipped | 969,232 |
| Protected values decoded | 0 |
| Non-TRAIN safe identities discarded | 205,581 |
| TRAIN raw rows retained | 279,035 |
| Malformed rows | 0 |
| Safe exact cells | 259,839 |
| Development outcomes read | 0 |
| Confirmation outcomes read | 0 |
| Sealed outcomes read | 0 |

The exact-cell counts reproduce the historical P0 topology without using its
outcome loader:

- pKi: 231,090;
- pKd: 28,749.

## Lineage recovery

Of 9,643 metadata documents:

- 8,901 contain a normalized DOI;
- 8,520 contain an exact PubMed ID;
- 740 contain an exact patent ID;
- 9,641 have at least one exact external lineage identifier;
- union closure gives 9,639 lineages.

Institution, experimental site, and patent-family fields are unavailable.
Empty document IDs remain unknown and are not promoted to independent
lineages. Because a blank pKi document value occurs on many cells, exact
lineage coverage by pKi cell is 91.1692%, below the frozen 95% gate. pKd cell
coverage is 100%.

## Frozen topology gates

| Endpoint | Exact cells | Joint components | Required | Largest component | Maximum | Lineage coverage | Decision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| pKi | 231,090 | 37 | at least 88 | 99.4690% | 20% | 91.1692% | fail |
| pKd | 28,749 | 65 | at least 88 | 95.5129% | 20% | 100% | fail |

The failure is unchanged after the homology-scaffold-lineage hypergraph
2-core:

- pKi: 211,599 cells, 38 components, largest 99.2410%;
- pKd: 25,412 cells, 52 components, largest 94.4436%.

Five nonempty lineage folds exist for both endpoints, but folds drawn inside
one 94-99% dependency component are not independent replication.

## Decision

Both endpoints fail two hard topology gates; pKi also fails lineage coverage.
`PCIC-O0-I` is therefore not executed. In particular:

- no CUDA operator residualization or SVD is run;
- no affinity is loaded;
- no direct operator, HCRR, or PB-CEC model is trained;
- numerical rank inside the giant component may not be used to override the
  component gate.

This is a failure of the current public measurement graph, not proof that a
target-ligand operator cannot exist.

## Required successor

The current public-data nuisance-reparameterization route is closed. Its
admissible successor is a label-blind, provenance-separated cycle-closing
acquisition design:

- complete target-by-ligand mini-blocks;
- new independent laboratory/campaign lineage;
- no reuse of database document IDs as independent replication;
- target homology and ligand scaffold components assigned as whole units;
- comparison of information-directed blocks with randomized complete blocks;
- a separate sealed provenance lineage for future confirmation.

