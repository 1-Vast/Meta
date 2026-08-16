# BioLiP2-D0 exact-complex layer preregistration

Date: 2026-07-29  
Status: frozen before remote readme/schema retrieval

## Role

BioLiP2 is evaluated only as an exact-complex/contact identity layer for
RBSDD. Its affinity annotations combine MOAD, PDBbind-CN, BindingDB, and
literature sources and are not presumed independent outcomes.

Allowed remote reads:

- official download-page HTML and HTTP headers;
- `readme.txt` and `readme_ligand.txt`;
- file names, sizes, modification metadata, and declared license;
- no `BioLiP.txt.gz`, ligand table row, structure archive, or affinity value.

## Frozen gates

All must pass before a D1 annotation download:

1. official primary annotation, ligand summary, receptor structures, and
   ligand structures are reachable with stable URLs;
2. annotation schema documents PDB ID, receptor chain, ligand CCD/chain or
   serial identity, binding residues, and UniProt or sequence identity;
3. ligand schema documents a connectivity-resolving field (canonical SMILES,
   InChI, or InChIKey) in addition to CCD ID;
4. annotation schema documents PubMed/reference and affinity-source fields so
   mixed affinity can be excluded and provenance closure can be audited;
5. the data/code license is explicit and compatible with a local derived
   metadata audit;
6. no annotation data row, affinity, development/confirmation data, or sealed
   test is read.

Failure returns `STOP_BIOLIP2_D0_SCHEMA_OR_PROVENANCE_INADEQUATE`.

Pass returns `REQUEST_BIOLIP2_D1_SAFE_ANNOTATION_DOWNLOAD`. It authorizes only
the primary annotation and ligand summary under a safe-column projection; all
structure archives and affinity values remain locked.

