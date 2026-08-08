# Data Availability and Repository Scope

This Git repository contains the MetaSieve source code, frozen mathematical
theory, contracts, tests, registered research protocols, terminal research
artifacts, compact PASS evidence, and provenance manifests.

Large upstream and generated datasets are intentionally not stored in Git.
The local workspace contains more than 50 GiB of third-party releases,
embeddings, caches, and compiled data. Those files include individual objects
larger than GitHub's 100 MiB limit and must be obtained from their original
sources or regenerated from the tracked workflows.

## Tracked data

- release, acquisition, corpus, split, governance, and cache manifests;
- upstream checksum and license records;
- compact P1B PASS evidence and checkpoint;
- registered research artifacts under `research/e0_identifiability/artifacts/`;
- reports required to reproduce the current evidence ledger.

Research artifacts are evidence for the registered experiments only. They do
not constitute a general-purpose affinity dataset or a validated end-to-end DTA
model.

## Data not redistributed

- ChEMBL 37 SQLite archive and extracted database;
- raw DAVIS and KIBA benchmark labels;
- BioLiP2/RCSB structure downloads and compiled structure corpora;
- ESM and ligand feature banks;
- OntoProtein, KeAP, ProteinKG25, PLINDER, and other third-party archives;
- downloaded MMseqs2 binaries and model caches.

## Primary upstream sources

- ChEMBL 37 static release: tracked in
  `dataset/raw/source_affinity/chembl37_sqlite_v1/release_manifest.json`.
  The manifest records the official archive SHA-256, schema hash, source URL,
  and CC BY-SA 3.0 license URL.
- Protein structures and chemical components: RCSB PDB/mmCIF and CCD, acquired
  by `scripts/acquire_open_structures.py`.
- BioLiP2 annotations: acquired and parsed through
  `scripts/structure_sources/biolip.py`.

## Reproduction boundary

The tracked manifests bind the expected release and generated artifacts, but
they do not grant authorization to run frozen downstream stages. Scientific
authorization is defined by `project_state.json`, `AGENT_HANDOFF.md`, and the
registered Gate documents. DAVIS and recipient-label access remains governed
by those contracts.

No private key, access token, local credential, or machine-specific SSH target
is included in this repository.
