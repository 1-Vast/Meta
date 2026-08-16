# A2S-DTA H0-S metadata inventory decision

**Date:** 2026-07-31  
**Status:** `DATA_NOT_READY`  
**Machine artifact:** `dataset/processed/a2s_h0_metadata_inventory.v1.json`  
**Artifact SHA-256:**
`b8812489f4d8535f28fbd8c3bce5c3d6c39539da03031e82634e958594ad7942`

## Firewall

**FACT:** H0-S read Parquet footer schema, JSON document metadata, ZIP central
directory metadata, and the first BindingDB TSV line only. It downloaded only
licenses, checksum lists, RDF dataset metadata, file descriptions, and ChEMBL
document fields. It read and materialized zero numeric affinity values and fit
no model. CUDA was not invoked.

## Source results

### ChEMBL HIST-S

**FACT:** The ChEMBL-37 registry SHA-256 is
`0e754f73f5d75913d61791d6ccd08e05662cd8015fc608ba370d4ee303e6b784`.
Its footer schema has no activity-row first-seen release field.

**FACT:** The original 9,643-document cache SHA-256 is
`5c920d4b33b88389c5331879ce1b620fe9b70e5006c3b3b845b1ad7102734109`.
The metadata-only API projection added publication year for 9,640 documents
and first-document-release date for all 9,643. Its SHA-256 is
`ecd0b47b58bf12171f188820c25ce77ecce82957dfe89bca460da1522a17c177`.
The three undated records are datasets `CHEMBL1201862`, `CHEMBL1909046`, and
`CHEMBL3885741`; they are ineligible for temporal episodes.

**FACT:** The latest official database files publicly available by the frozen
index dates are ChEMBL 24.1 for 2018-12-31, ChEMBL 27 for 2020-12-31, and
ChEMBL 31 for 2022-12-31. Their official licenses, attribution records, and
checksum manifests are frozen locally. The database snapshots themselves are
not present.

**INFERENCE:** A document's first ChEMBL release cannot prove that every
current ChEMBL-37 activity under that document existed in that release. Using
document release as activity time would admit database backfill and invalidate
the historical deployment estimand. HIST-S therefore cannot enter H1-S yet.

### BindingDB HIST-L candidate

**FACT:** The local BindingDB native-article archive is version 202607, has
SHA-256
`d2584d1519318d00ab5f46289da5ab3549affe732d598a5072f8777b6b3b5262`,
and exposes 640 header fields. The header includes record ID, InChIKey,
curation source, DOI, PMID, patent, publication date, BindingDB date,
institution, PubChem assay ID, UniProt fields, and target sequences. No data
row was read.

**FACT:** The version-matched reaction-set-to-assay and assay-description
companion archives are remotely available. They are not yet frozen locally.
The main TSV alone has no generic assay-lineage key.

**INFERENCE:** BindingDB remains a HIST-L candidate. It cannot repair HIST-S,
authorize H1-S, or be counted as independent evidence until companion assay
metadata and upstream-origin coverage are audited.

### GtoPdb HIST-L candidate

**FACT:** Official GtoPdb RDF metadata identifies version 2026.2, issued
2026-06-15, with ODbL database terms and CC BY-SA 4.0 content terms. The
frozen metadata hashes are:

- `gtp2026.2.ttl`:
  `27e5435c7ea6632d5d72d9c21ca9cb135703bbcf01fc7ca98a49569b909358d4`;
- `file_descriptions.txt`:
  `7bf4e6876ca11b693db9f6f6a5ca21f8273ac0b7b071acd149b447c885442daa`.

**FACT:** The schema documentation exposes PubMed, patent, assay-description,
and InChIKey fields, but the current release metadata does not provide a
per-interaction first-seen release history.

**INFERENCE:** GtoPdb is localized for later L0-L audit only. It does not
authorize current label access or model training.

## Decision

**DATA STOP:** Document dates are now largely complete, but true activity-row
availability at the historical index dates is still unidentifiable from the
ChEMBL-37 final snapshot. H0-S does not authorize H1-S, S0-S, A2S-MAP, MAML,
AdaMBind, or any other affinity fit.

The unique next action is to acquire and checksum the official ChEMBL 24.1,
27, and 31 historical databases, then project only stable activity identities
and provenance metadata needed to test first appearance. The projection must
not select, print, aggregate, or materialize affinity values. Large-model
training remains prohibited on this workstation.
