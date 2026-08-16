# DCST-R18 BindingDB PDB-linked source preregistration

Date: 2026-07-29  
Status: frozen before safe-column scan

## Candidate

The local BindingDB 2026-07 native-article archive is an independent,
staff-curated source rather than a ChEMBL export. Its row schema contains:

- ligand SMILES and identifiers;
- target-chain sequences and UniProt identifiers;
- PDB IDs for ligand-target complexes;
- article DOI/institution provenance;
- Ki and Kd fields.

R18 asks only whether enough PDB-linked, single-chain, exact-endpoint records
exist to justify a new Stage-1 source build. It does not read numeric affinity.

## Byte-level outcome firewall

The scanner may decode only the header and these safe fields:

- row ID, ligand SMILES/InChIKey, target name and organism;
- curation/data source, DOI, institution;
- ligand HET ID, complex PDB ID;
- target-chain count, chain-1 sequence, chain-1 PDB IDs, and chain-1 UniProt.

For Ki and Kd, it may record only a boolean `field_nonempty`; the bytes must
not be decoded, parsed, copied, hashed, logged, or written. IC50, EC50,
kinetics, pH, and temperature are ignored.

Rows must be native BindingDB provenance, single-chain, have exactly one
usable chain sequence and ligand, an explicit complex PDB ID, and a nonempty
Ki or Kd field.

## Cross-source firewall

Before support counts, exclude rows with:

- ChEMBL development/confirmation exact accession;
- sequence 4-mer containment at least `0.40` to a protected target;
- exact protected ligand connectivity/scaffold;
- Morgan Tanimoto at least `0.95` to a protected ligand;
- missing DOI or a curation source containing `ChEMBL`.

PDB structures and coordinates are not downloaded during R18.

## Frozen gate

The source is eligible for an authorized construction request only if the
firewalled safe registry has:

1. at least 2,000 rows, 150 accessions, 1,500 ligands, 1,500 PDB IDs, and 300
   article DOIs;
2. at least 500 exact `(accession, ligand)` pairs with a matching ChEMBL TRAIN
   pair, spanning 100 targets and 80 ChEMBL homology components;
3. exact accession support for at least 20% of ChEMBL-train targets;
4. maximum source Morgan Tanimoto at least `0.40` for at least 20% of a
   deterministic 20,000-ligand ChEMBL-train sample;
5. no protected-entity or ChEMBL-provenance violation.

Failure returns `STOP_BINDINGDB_PDB_SOURCE_INADEQUATE`.

Pass returns `REQUEST_BINDINGDB_STAGE1_AUTHORIZATION`; it does not itself
authorize reading numeric affinities or downloading PDB coordinates. Those
are meaningful new data access and require a separately frozen construction
manifest and explicit user authorization.

