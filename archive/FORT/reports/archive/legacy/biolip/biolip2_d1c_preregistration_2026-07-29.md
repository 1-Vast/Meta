# BioLiP2-D1C provenance/promiscuity closure preregistration

Date: 2026-07-29  
Status: frozen before safe-registry closure

## Motivation

D1 passed source scale, same-target pair, ChEMBL support, and firewall gates,
but its ordinary target-ligand core was connected by generic/promiscuous
ligands and its largest PubMed block narrowly exceeded the cap. D1C does not
change either failed D1 threshold. It asks whether a source-native
independence closure produces an adequate substrate.

## Frozen closure

Starting only from the D1 safe registry:

1. compute ligand degree as the number of sequence-exact targets;
2. remove every ligand with degree greater than 50;
3. collapse repeated structures to one deterministic lowest-lexicographic PDB
   representative per `(PubMed, target_key, ligand connectivity)`;
4. recompute same-target distinct-PDB ligand pairs, target-ligand 2-core,
   connected components, ChEMBL support, and PubMed concentration.

No affinity value or structure archive may be loaded. The removed generic
ligands and collapsed rows must be counted and hashed.

## Frozen gates

Retain the D1 thresholds without relaxation:

1. at least 20,000 closed rows, 1,000 targets, 10,000 ligands, 10,000 PDB
   entries, and 5,000 PubMed IDs;
2. at least 500 targets and 20,000 same-target distinct-PDB ligand pairs;
3. 2-core at least 5,000 edges/300 targets/2,000 ligands, more than one
   component, and largest component at most 50% of edges;
4. ChEMBL-TRAIN exact target support and ligand support each at least 20%;
5. all rows have PubMed and the largest PubMed contributes at most 2%;
6. zero affinity, protected feature/label, structure archive, or sealed value.

Failure returns `STOP_BIOLIP2_D1C_INDEPENDENT_TOPOLOGY_INADEQUATE`.

Pass returns `REQUEST_BIOLIP2_D2_DIRECTED_EDIT_AUDIT`; it still does not
authorize structure downloads or RBSDD training.

