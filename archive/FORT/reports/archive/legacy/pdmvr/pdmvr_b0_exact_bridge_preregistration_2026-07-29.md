# PD-MVR-B0 exact bridge preregistration

Date frozen: 2026-07-29  
Status: frozen before exact triple-join enumeration  
Purpose: determine whether the 348 R18 ChEMBL identity links hide at
least 40 exact, provenance-separable structure-affinity bridge candidates.

## Frozen sources

- BindingDB native-article archive:
  `dataset/public/open_s/BindingDB_BindingDB_Articles_202607_tsv.zip`
  (`d2584d1519318d00ab5f46289da5ab3549affe732d598a5072f8777b6b3b5262`)
- BioLiP2 closed registry:
  `dataset/public/biolip2/processed/closed_registry.parquet`
  (`7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`)
- Safe ChEMBL document metadata:
  `dataset/public/chembl_37/processed/dualcold/pcic_o0_document_metadata.json`
  (`5c920d4b33b88389c5331879ce1b620fe9b70e5006c3b3b845b1ad7102734109`)
- ChEMBL dual-cold registry, used only for the existing
  development/confirmation entity firewall and TRAIN homology metadata:
  `dataset/public/chembl_37/processed/dualcold/registry.parquet`
  (`0e754f73f5d75913d61791d6ccd08e05662cd8015fc608ba370d4ee303e6b784`)

## Exact bridge candidate

A BindingDB row is only a candidate when all are true:

1. it passes the already implemented R18 protected-entity and provenance
   firewall;
2. the BindingDB row records Ki or Kd presence without decoding the value;
3. one listed PDB ID exactly matches a BioLiP2 PDB ID;
4. canonical ligand connectivity is exact;
5. receptor sequence is exact;
6. accession sets overlap.

The audit will additionally compare the BindingDB DOI with the DOI mapped
from BioLiP2 PubMed through the frozen local ChEMBL document metadata.
Matching DOI is a verified same-article subset. Missing local PubMed-to-DOI
coverage remains unverified rather than being treated as a mismatch.

## Frozen gates

Before any coordinates, contacts, or numeric affinity may be read, the
following optimistic upper-bound gates must all pass:

- at least 40 exact PDB-sequence-accession-ligand bridge units;
- at least 40 sequence-exact targets;
- at least 40 provenance-conflict-free units under one unit per target,
  ligand, BindingDB DOI, institution, and BioLiP PubMed;
- at least 4 broad families, if a safe family mapping is available;
- zero protected outcome values loaded.

The conflict-free ceiling is the minimum of the corresponding unique
resource counts. This is optimistic because it does not yet close homology,
authors, institutions, assays, chemical neighbours, constructs, or related
publications.

If any available upper bound is below 40, stop PD-MVR bridge anchoring. If
all available upper bounds pass but DOI/family coverage is incomplete,
request only the missing metadata audit. Coordinates and bridge-rank
calculation remain locked until exact same-article provenance and family
gates pass.

