# RDIB-Edit-G0 decision

Date: 2026-07-29  
Decision: `STOP_RDIB_EDIT_G0_REPLICATION_TOPOLOGY_INADEQUATE`

## Corrected result

The formal reproduction used the existing frozen MMP-X single-cut rules on
the BioLiP2 closed registry. It enumerated 133,352 source-specific ligand
pairs:

- 2,530 pairs had at least one eligible MMP token;
- 537 multi-token ambiguous pairs were excluded;
- 1,993 unique-token source records represented 730 global edit classes;
- 200 edit classes appeared somewhere in at least two PubMeds and at least
  two targets.

The global count is not replicated target-conditioned supervision. Under the
correct biological unit—one `(exact target_key, edit)` supported by at least
two distinct PubMeds—only 29 units remained:

- 17 sequence-exact targets and 17 edit classes;
- 40 PubMed IDs, 95 ligands, and 61 PubMed-unit incidences;
- 6 units with an exact ligand pair repeated across PubMeds;
- 23 edit-only salvage units using different molecule pairs;
- largest PubMed unit share `10.3448%`.

The optimistic resource ceiling was only 20:

`min(29 units, 3*17 targets, floor(40/2) PubMeds, floor(95/2) ligands)`.

This is below the inherited 88-unit optimistic mechanism floor, as well as
the 200-block structural floor and approximately 423-unit predictive scale.

## Why no stricter identity run is needed

The source connectivity was generated with `isomericSmiles=False`. The
formal result is therefore an optimistic non-isomeric upper bound, not a
final chemical-identity pass. Reconstructing stereo, charge, tautomer,
scaffold, chemical-neighbour, exact construct, homology, author,
institution, and deposition-series conflicts can only remove or merge
eligible units. None can raise the ceiling from 20 to 88.

## Consequence

The recurring-directed-edit secondary route is stopped before original
ligand identity reconstruction, coordinate download, PLIP, contact ICC,
operator-rank calculation, affinity access, or model training. The existence
of a 730-token edit vocabulary must not be confused with replicated
target-conditioned interaction evidence.

No affinity field, development/confirmation feature or outcome, or sealed
outcome was loaded. The audit was CPU-only RDKit fragmentation and set
enumeration.

Authoritative machine result:
`reports/active/rdib_edit_g0.json`.

