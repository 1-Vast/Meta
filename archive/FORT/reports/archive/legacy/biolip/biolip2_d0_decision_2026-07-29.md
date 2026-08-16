# BioLiP2-D0 decision

Date: 2026-07-29  
Decision: `REQUEST_BIOLIP2_D1_SAFE_ANNOTATION_DOWNLOAD`

All six remote readme/schema gates passed.

- Primary annotation: 78,570,940 bytes, ETag
  `"4aee5bc-6435c2eb61100"`.
- Ligand summary: 8,065,879 bytes, ETag
  `"7b1357-6435bad8c4540"`.
- Exact identity fields: PDB ID, receptor chain, ligand CCD ID/chain/serial,
  receptor sequence, UniProt, and PubMed.
- Contact fields: PDB-numbered and reindexed binding residues.
- Connectivity fields: InChI, InChIKey, and one or more SMILES strings.
- Mixed affinity sources are separate columns (manual literature, MOAD,
  PDBbind-CN, and BindingDB), so all can be excluded by a byte-level safe
  projection.
- Curation repository commit:
  `f56c011b39273779585e4862db1e1033be0b4234`, BSD-2-Clause with a root
  LICENSE file.

No annotation or ligand data row and no affinity value was read. No structure
archive was downloaded. D0 authorizes only a separately preregistered safe
download/projection of the annotation and ligand summary.

Authoritative machine result: `reports/active/biolip2_d0.json`.

