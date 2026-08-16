# RBSDD source-recovery amendment

Date: 2026-07-29  
Status: binding amendment

## Evidence alignment

The user-supplied audit is consistent with the completed gates:

- R16 learned a protein-local response but failed ligand-pair specificity;
- R17 had stable KLIFS measurements but inadequate independent rectangle
  topology and ligand support;
- R18 retained 6,561 PDB identifiers but only 340 independent ligands, and
  official BindingDB documentation confirms that PDB links are based on
  ligand chemical identity and protein sequence identity rather than proof
  that the affinity was measured on that exact crystal complex;
- R19/R20 showed that raw combinatorial affinity rectangles collapse under
  concentration, endpoint, and provenance controls.

Therefore no current local source authorizes another PLINDER/KLIFS adapter,
structure teacher, affinity curriculum, or representation alignment.

## New admissible model

The only structural model candidate is Real Bilateral Structural-Delta
Distillation (RBSDD). Its source unit must contain one target and two real
ligand complexes with verified chain/ligand correspondence. The source must
support a target-by-directed-edit graph with repeated edits across independent
targets or families.

RBSDD model construction remains locked until a label-blind source gate
establishes:

- exact complex identity, not similarity-based PDB links;
- repeated same-target ligand pairs and directed chemical edits;
- target/edit 2-core, multiple independent components, and adequate effective
  rank;
- source lineage and document/institution closure;
- target and ligand support relative to ChEMBL TRAIN;
- wrong-target, wrong-edit, and ligand-edit-only destructive controls.

## Source priority

1. PBCNet2.0 Zenodo/GitHub source audit;
2. BioLiP2 as an exact-complex/contact identity layer, not an independent
   affinity source;
3. PSICHIC overlap and pair-specific attribution audit;
4. if these fail, close public structure-source replacement and return to the
   already-defined non-structural recovery program.

