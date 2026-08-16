# RDIB / PD-MVR G0-A decision

Date: 2026-07-29  
Decision: `STOP_RDIB_PDMVR_G0A_REPLICATION_TOPOLOGY_INADEQUATE`

## Result

The affinity-blind screen enumerated all eligible same-target ligand pairs
from the 66,660-row BioLiP2 closed registry. A biological candidate was
counted only when two different PubMed records each contained the same
sequence-exact target and both exact ligand connectivities.

| Quantity | Observed | Frozen requirement |
| --- | ---: | ---: |
| target-PubMed panels with at least two ligands | 8,512 | descriptive |
| source-specific ligand pairs | 133,352 | descriptive |
| target-ligand observations repeated across PubMeds | 1,028 | descriptive |
| replicated exact ligand-pair difference blocks | 145 | at least 200 |
| sequence-exact targets | 59 | at least 80 |
| optimistic conflict-free packing ceiling | 56 | at least 200 |
| largest PubMed candidate-block share | 17.2414% | at most 5% |

The optimistic ceiling was already limited to 56 by the 113 distinct
ligands (`floor(113/2)`). It ignored homology, family, chemical-neighbour,
author, institution, construct, cofactor, assay, and related-publication
closure, so every stricter packing can only be smaller.

## Target-domain support

Only 8/145 candidates had target-accession support in ChEMBL TRAIN, 6/145
had both ligands present anywhere in TRAIN, and only 5/145 (`3.4483%`) had
both ligands recorded against the same ChEMBL TRAIN target. The replicated
structure differences therefore do not span the intended downstream query
domain even before contact extraction.

## Consequence

- RDIB's exact ligand-pair primary route is stopped before coordinates,
  PLIP, contact ICC, operator-rank computation, or model training.
- PD-MVR fails its independent-structure-block requirement on the same
  evidence.
- The 391,262 historical same-target pairs must not be cited as independent
  replicated differences.
- A recurring directed-edit secondary route is not inferred from this
  result. It requires a separately frozen chemical-standardization,
  atom-mapping, direction, and edit-equivalence contract.

No affinity field, development/confirmation feature or outcome, sealed-test
outcome, or structure coordinate was loaded. The audit used CPU because it
was exact set enumeration rather than a tensor operation.

Authoritative machine result:
`reports/active/rdib_pdmvr_g0a.json`.

