# BioLiP2-D1 safe exact-complex topology preregistration

Date: 2026-07-29  
Status: frozen before annotation/ligand download

## Scope and safe projection

Download only:

- `download/BioLiP.txt.gz` (78,570,940 bytes at D0);
- `data/ligand.tsv.gz` (8,065,879 bytes at D0).

Do not download receptor or ligand structure archives.

For the 21-column annotation, decode only columns 1-13 and 18-21. Columns
14-17 (manual, MOAD, PDBbind-CN, BindingDB affinity) must not be decoded,
copied, hashed, logged, or written; only a four-bit nonempty presence mask may
be retained for audit counts.

The ligand summary may decode CCD ID, InChI, InChIKey, and SMILES. Name and
external registry IDs are not needed.

## Registry and firewall

Retain rows with:

- receptor sequence length at least 40;
- a canonicalizable small-molecule ligand with at least six heavy atoms;
- explicit PDB, receptor chain, ligand CCD/chain/serial, binding residues,
  and PubMed;
- exact sequence hash as the source target identity.

Before topology/support measurement, exclude any row matching ChEMBL
development/confirmation by:

- exact accession;
- sequence 4-mer containment at least `0.40`;
- exact ligand connectivity or scaffold;
- Morgan Tanimoto at least `0.95`.

Confirmation is used for firewall metadata only; no protected feature or
label may be loaded.

## Frozen D1 gates

All must pass:

1. at least 20,000 firewalled exact-complex rows, 1,000 sequence-exact
   targets, 10,000 ligands, 10,000 PDB entries, and 5,000 PubMed IDs;
2. at least 500 targets have two or more unique ligands in distinct PDB
   entries, yielding at least 20,000 same-target ligand pairs;
3. target-ligand bipartite 2-core has at least 5,000 edges, 300 targets,
   2,000 ligands, and more than one connected component; the largest
   component contributes at most 50% of edges;
4. exact accession support covers at least 20% of ChEMBL-TRAIN targets and
   maximum source Morgan Tanimoto at least `0.40` covers at least 20% of the
   frozen 20,000-ligand ChEMBL-TRAIN sample;
5. at least 60% of retained rows have PubMed provenance and no single PubMed
   ID contributes more than 2% of rows;
6. zero protected rows, decoded affinity values, structure archives,
   development/confirmation features or labels, and sealed-test values.

Failure returns `STOP_BIOLIP2_D1_EXACT_COMPLEX_TOPOLOGY_INADEQUATE`.

Pass returns `REQUEST_BIOLIP2_D2_DIRECTED_EDIT_AUDIT`. A pass authorizes only
construction of a label-free atom-correspondence/edit registry. It does not
authorize structure downloads or RBSDD training.

