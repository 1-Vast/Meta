# FACTOR-U U0 preregistration: strict unlabeled-corpus eligibility audit

Date: 2026-07-26  
Route origin: user-supplied post-F0-C option (“expand unlabeled chemistry after strict
firewall”); consumes no agent-proposed candidate slot.  
Role: data-only audit. No encoder or affinity model is trained in U0.

## Motivation and non-tuning boundary

F0-C1 learned real local chemistry but had participation rank 8--10, atom-decoy false coverage
0.066--0.075 and inner scaffold-OOD median 0.774--0.820. U0 does not modify that model. It tests the
substantively different hypothesis that the 2,489-molecule development corpus is too small to
identify a broad continuous carrier geometry.

The candidate expansion source is the local PLINDER 2024-06/v2 annotation snapshot, whose project
manifest records CC BY 4.0. Only ligand structure and ligand quality flags may be projected from the
Parquet file. Affinity, protein, pocket, split, source and interaction columns are forbidden.

## Frozen ligand filter

Read only:

- `ligand_rdkit_canonical_smiles`;
- `ligand_is_rdkit_loadable`;
- `ligand_is_cofactor`, `ligand_is_artifact`, `ligand_is_ion`, `ligand_is_oligo`,
  `ligand_is_invalid`.

Keep RDKit-loadable, non-cofactor, non-artifact, non-ion, non-oligo and non-invalid records. After
RDKit parsing, keep 6--80 heavy atoms and elements in
`{B,C,N,O,F,Si,P,S,Cl,Se,Br,I}`. Canonicalize without stereochemistry and deduplicate parent
connectivity.

## Global strict firewall

Build the union of all KIRHub2026, Reinecke2024 and Papyrus-Christmann2016 evaluation
connectivities and Bemis--Murcko scaffolds using the already frozen label-blind loaders. Delete a
PLINDER molecule if either its connectivity or scaffold occurs anywhere in that union. This is
stricter than per-fold removal and gives every later fold the same untouched pretraining corpus.

No existing ChEMBL confirmation record or label may be opened. PLINDER test/confirmation metadata
are irrelevant because U0 projects no PLINDER split or affinity field. The project-wide historical
ChEMBL confirmation contamination remains recorded and quarantined.

## Eligibility gates

All must pass:

1. zero parent-connectivity overlap with every evaluation source;
2. zero Bemis--Murcko scaffold overlap with every evaluation source;
3. at least 20,000 unique retained molecules (more than ten times the largest F0-C1 inner-train
   fold);
4. at least 8,000 distinct nonempty Murcko scaffolds;
5. no single scaffold contributes more than 0.5% of retained unique molecules;
6. every element class and every RDKit pharmacophore role observed in the evaluation corpus occurs
   in the retained corpus;
7. retained heavy-atom q01 is no larger than evaluation q05 and retained q99 is no smaller than
   evaluation q95;
8. source license and raw-file SHA-256 are recorded;
9. current-run activity/affinity/protein columns read = false;
10. current-run confirmation labels read = false and sealed test consumed = false.

Pass: `FACTOR_U0_PASS_AUTHORIZE_U1_PREREGISTRATION`.  
Fail: `FACTOR_U0_CORPUS_INELIGIBLE_STOP`.

If U0 passes, U1 must be separately preregistered. U1 may test public-unlabeled pretraining but may
not reuse F0-C1 checkpoints, change external coverage gates, use activity labels or claim strict
corpus-disjoint success without the global firewall above. The local pair cap remains an audit
weighting bound; a passed U1 requires uncapped replication on the larger machine.
