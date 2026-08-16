# BioLiP2-D1C decision

Date: 2026-07-29  
Decision: `STOP_BIOLIP2_D1C_INDEPENDENT_TOPOLOGY_INADEQUATE`

## Closure

D1C removed 182 ligands observed against more than 50 sequence-exact targets
and 443,978 associated rows. It then collapsed 52,156 repeated
`(PubMed, target, ligand)` structures, leaving:

- 66,660 rows;
- 40,969 exact targets, 25,648 ligands, 42,887 PDB entries;
- 19,764 PubMed IDs;
- 6,284 targets and 391,262 same-target distinct-PDB ligand pairs.

Scale, same-target pairs, ChEMBL support, provenance, and firewall remained
above the unchanged D1 thresholds. ChEMBL target/ligand support was
`70.3041%/49.7900%`, and the largest PubMed block fell to `0.8161%`.

## Failure

The closed 2-core still has a largest component of 14,387/17,989 edges
(`79.9767%`), above the frozen 50% limit. The giant component therefore is
not explained solely by the 182 most promiscuous ligands or by repeated
structures from one publication.

## Consequence

BioLiP2 is retained as a valuable exact-complex/contact cross-validation
layer, but it is not authorized as an independent RBSDD Stage-1 source. D2
directed-edit construction and all structure downloads remain locked.

Authoritative machine result:
`reports/active/biolip2_d1c_seed1729.json`.

