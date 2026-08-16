# DCST-R2 entity-aligned privileged transfer preregistration

Date: 2026-07-28  
Status: frozen before RCSB entity retrieval and R2 training

## Hypothesis

The principal Stage-1 information loss is not model width. It is the
conflation of a PLINDER split cluster with the protein actually present in a
structure–ligand row, followed by an invalid projection of PDB residue
coordinates onto that cluster representative.

R2 tests **entity-aligned privileged transfer**:

1. a PLINDER cluster remains only a firewall, split, and cross-fitting unit;
2. the model target key is the row's unique `system_pocket_UniProt`;
3. the protein input is the frozen ESM-2 representation of that exact
   accession;
4. each structural residue is mapped from PDB polymer-entity sequence
   coordinates to the same UniProt accession through RCSB's SIFTS aligned
   regions;
5. only the resulting correctly aligned 32-segment × 8-interaction-type map
   supervises Stage 1;
6. source spectral directions still require target- and ligand-destruction
   certification before Stage-2 exposure.

This separation of split identity, model identity, privileged structural
coordinates, and transfer certificate is the candidate innovation. Generic
pretraining or fine-tuning is not the novelty claim.

## Frozen source registry

- parent rows: the existing PLINDER train/development rows after the frozen
  cross-source firewall;
- admissible target: exactly one `system_pocket_UniProt`, consistent with the
  row's processed accession set and present in the local UniProt sequence
  registry;
- observed pre-retrieval coverage: 2,124 train rows and 220 development rows,
  covering 863 exact accessions;
- multi-accession or unmapped rows are not silently assigned a representative;
- RCSB reads are limited to `entry_pdb_id/polymer_entity_id`, canonical entity
  sequence, UniProt identifiers, and SIFTS aligned regions; no affinity source
  is queried;
- ESM model, maximum length 1,022, ligand representation and 32-segment rule
  are unchanged from R1.

Rows whose contact residue cannot be mapped into the exact accession and the
embedded 1,022-residue prefix are excluded from the privileged target and
reported. Source affinity training and certificates group by exact accession,
while cross-fitting remains blocked by the original transitive target
component.

## Frozen source-only gate

Seed 1729 and 4,000 source-base/interaction steps are fixed. All must hold:

1. true centered joint alignment is positive;
2. true alignment exceeds wrong-target alignment by more than `0.05`;
3. true alignment exceeds within-target wrong-ligand alignment by more than
   `0.05`;
4. true absolute joint cross-entropy is below uniform `log(256)`;
5. at least one two-direction source affinity band is certified;
6. the privileged model certifies strictly more bands than the matched
   exact-target `DCST-NoPriv` model.

R2 stops before any new downstream-label load if this gate fails. A pass
authorizes generation of the matching ChEMBL target cache and the unchanged
registered Stage-2 arm comparison; confirmation and sealed evaluation remain
closed.

