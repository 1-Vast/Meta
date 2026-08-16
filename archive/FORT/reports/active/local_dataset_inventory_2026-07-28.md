# Local Public Dataset Inventory

Date: 2026-07-28

Status: metadata-only inventory. This report records what is physically present
under `D:\FORT\dataset\public`, what each source may currently be used for, and
what it must not be used for. It is not an authorization to train on real
affinity labels.

## 1. Audit boundary

This inventory inspected:

- directory and file names, byte sizes, and file formats;
- JSON manifests and schema records;
- README, license, checksum, and acquisition records;
- Parquet footer schemas and footer row counts;
- XLSX workbook and sheet metadata;
- text-file header rows; and
- ZIP central-directory metadata.

No affinity row or affinity value was opened, printed, sampled, aggregated, or
materialized during this inventory. Counts and prior consumption flags below
come from existing manifests, Parquet footers, and frozen decision reports.
They describe earlier program state, not label access by this audit.

The inventory root contains 115 files totaling 1,470,606,023 bytes
(1,402.48 MiB; 1.37 GiB).

## 2. Physical inventory

| Top-level directory | Exact local root | Files | Bytes | MiB | Current program role |
| --- | --- | ---: | ---: | ---: | --- |
| `chembl_37` | `D:\FORT\dataset\public\chembl_37` | 20 | 123,867,585 | 118.13 | Historical exact-pK registry and frozen engineering substrates |
| `kirhub_2026` | `D:\FORT\dataset\public\kirhub_2026` | 6 | 17,824,560 | 17.00 | Historical ordinal mutation and ligand-reordering evidence |
| `klifs_2026_07_22` | `D:\FORT\dataset\public\klifs_2026_07_22` | 11 | 3,748,559 | 3.57 | Pocket alignment and kinase mechanism infrastructure |
| `open_s` | `D:\FORT\dataset\public\open_s` | 1 | 18,114,757 | 17.28 | BindingDB source-recovery archive; affinity rows remain unopened |
| `openbind_ev_a71_2026` | `D:\FORT\dataset\public\openbind_ev_a71_2026` | 44 | 12,867,576 | 12.27 | Single-protein local-SAR and structure benchmark |
| `papyrus_05_7` | `D:\FORT\dataset\public\papyrus_05_7` | 3 | 58,552,385 | 55.84 | Aggregated-source provenance audit; stopped before training |
| `plinder_2024_06_v2` | `D:\FORT\dataset\public\plinder_2024_06_v2` | 8 | 1,175,433,550 | 1,120.98 | Historical structure-affinity engineering diagnostics |
| `reinecke_2024` | `D:\FORT\dataset\public\reinecke_2024` | 11 | 29,733,405 | 28.36 | Development-only kinase pKd_app panel |
| `spd_2023` | `D:\FORT\dataset\public\spd_2023` | 7 | 28,871,214 | 27.53 | Systematic, censored multi-family source and power audit |
| `toxcast_invitrodb_v4_3` | `D:\FORT\dataset\public\toxcast_invitrodb_v4_3` | 4 | 1,592,432 | 1.52 | Target-contrast pretext metadata and lightweight projection |
| **Total** | `D:\FORT\dataset\public` | **115** | **1,470,606,023** | **1,402.48** | Mixed evidence store; not one exchangeable DTA table |

File counts are recursive and count files only. In particular,
`D:\FORT\dataset\public\openbind_ev_a71_2026\raw` exists but is empty.

## 3. Dataset records

### 3.1 ChEMBL 37

**Physical content**

- Exact root: `D:\FORT\dataset\public\chembl_37`.
- Formats: 2 GZIP files, 8 JSON files, 6 NPZ files, and 4 Parquet
  files.
- There is no protected `raw` source directory in this snapshot. The local
  object is a processed API extract plus derived registries and features.
- Source record:
  `D:\FORT\dataset\public\chembl_37\processed\api_source_manifest.json`.
- Current main registry:
  `D:\FORT\dataset\public\chembl_37\processed\dualcold\registry.parquet`.
- Nested panel registries:
  `D:\FORT\dataset\public\chembl_37\processed\panel_metz\registry.parquet`
  and
  `D:\FORT\dataset\public\chembl_37\processed\panel_davis\registry.parquet`.

**Representative schema metadata**

- The main Parquet footer reports 343,211 rows with target, ligand parent
  connectivity, endpoint, affinity, replicate reliability, scaffold,
  organism, assay, document, UniProt accession, homology cluster, chemical
  similarity, and frozen split fields.
- Its manifest reports 201,827 TRAIN pairs, 559 TRAIN targets, 517 TRAIN
  homology clusters, 121,401 TRAIN ligand parents, and 48,234 TRAIN
  scaffolds.
- The Metz registry footer reports 29,291 cells from one pKi document. Its
  historical TRAIN block contains 12,574 cells, 112 targets, 101 homology
  components, and 619 ligands.
- The Davis registry footer reports 3,429 cells from one pKd document. It is
  marked `single_use=true` and `role=confirmation` in its manifest.
- NPZ objects hold frozen ligand features and pooled target ESM-2 features;
  they are derived features, not independent measurements.

**Provenance and license**

- The recorded upstream endpoint filter is human single protein, exact Ki or
  Kd, nM units, equality relation, and positive value. pIC50 is explicitly
  excluded.
- The manifest records ChEMBL CC BY-SA 3.0 and requires attribution,
  share-alike handling, and a terms check before redistribution.
- Absence of a local raw database dump means the API source manifest and
  processed-object hashes are essential reconstruction evidence.

**Authorized role**

- Main dual-cold TRAIN: historical estimator engineering and label-blind
  topology audits only under the current program stop.
- Metz: historical, one-document, kinase-only train/development mechanism
  substrate. It is not fresh confirmation evidence.
- Davis: protected single-use target-conditioned confirmation candidate,
  pending the irreversible Davis overlap decision in the active task ledger.

**Prohibited use**

- Do not restart real affinity training on this observational graph to rescue
  a failed coordinate.
- Do not treat database rows, documents, assays, or replicated database
  exports as independent biological measurement sites.
- Do not merge pKi, pKd, pIC50, or other endpoint types into one unqualified
  label.
- Do not use Metz or previously inspected main-registry development material
  as independent confirmation.
- Do not read Davis target-conditioned labels before its single-use protocol
  and the DAVIS-Complete overlap policy are frozen.

**Firewall and consumption state**

- The main manifest records target, UniProt, homology-cluster, scaffold,
  ligand-connectivity, high-similarity, assay, and document separation.
- The current research ledger also records a giant ligand/provenance
  component and residual observational-source confounding. Manifest-clean
  axes therefore do not establish a compliant new factorial substrate.
- The local main manifest says `confirmation_labels_read=false`, but the
  program-level ledger is stricter and records historical ChEMBL confirmation
  access as true because five confirmation rows were displayed during an
  earlier schema inspection. The program-level record takes precedence.
- Davis target-conditioned confirmation remains unconsumed. Historical
  arm-blind Davis power labels were read; this does not authorize a
  target-conditioned run.
- `sealed_test_consumed=false`.

### 3.2 KirHub 2026

**Physical content**

- Exact root: `D:\FORT\dataset\public\kirhub_2026`.
- Formats: 1 XLSX workbook, 4 JSON files, and 1 NPZ feature file.
- Raw workbook:
  `D:\FORT\dataset\public\kirhub_2026\raw\41587_2026_3090_MOESM4_ESM.xlsx`.
- Derived registry:
  `D:\FORT\dataset\public\kirhub_2026\processed\h0_component_registry.json`.
- Strict connected components:
  `D:\FORT\dataset\public\kirhub_2026\processed\strict_components.json`.

**Representative schema metadata**

- The workbook has 15 sheets, including `Table S3`, `Table S4`, and
  `Table S13`.
- Frozen reports identify Table S3 as a dose-response block, Table S4 as a
  409-preparation by 92-inhibitor wild-type single-dose block, and Table S13
  as a mutant/fusion single-dose block.
- The H0 registry declares 409 targets and 92 ligands. Derived sequence maps
  cover 377 resolved gene entries, and pooled ESM-2 coordinates are stored in
  `target_esm2.npz`.
- The released activity cells are aggregate residual-activity measurements;
  raw duplicate observations are not present locally.

**Provenance and license**

- The source is the supplementary workbook for Saifudeen et al., Nature
  Biotechnology (2026).
- No explicit reusable data-license identifier is frozen in the local dataset
  tree. The active preregistration permits internal academic analysis under
  source terms and explicitly forbids redistribution. License status is
  therefore **not verified for redistribution**.

**Authorized role**

- Historical ordinal mutation double-difference and target-specific
  ligand-reordering audits.
- Historical frozen-coordinate necessity and pocket-alignment oracle audits.
- Evidence that biological reordering exists locally, not an affinity
  training source and not cross-family confirmation.

**Prohibited use**

- Do not convert percent residual activity at one concentration into Ki, Kd,
  or a continuous exact-pK endpoint.
- Do not count construct-ligand pairs as independent samples; the defensible
  units are genes, families, or homology components.
- Do not claim cross-assay, cross-document, or cross-family transfer from this
  single kinase source.
- Do not add model capacity or tune a protein coordinate after its registered
  fail-stop results.

**Firewall and consumption state**

- The authoritative strict split uses full-sequence homology components and
  chemical connected components joining parent identity, equal Murcko
  scaffold, or Morgan Tanimoto at least 0.50.
- Assay, document, and measurement-site independence cannot be tested within
  this one-source release.
- Activity labels were consumed in historical mechanism audits. This source
  has no protected independent confirmation role.
- Program sealed tests remain unconsumed.

### 3.3 KLIFS 2026-07-22

**Physical content**

- Exact root: `D:\FORT\dataset\public\klifs_2026_07_22`.
- Formats: 5 GZIP JSON snapshots, 3 JSON files, 2 Parquet registries, and 1
  NPZ ligand-feature file.
- Raw API snapshot:
  `D:\FORT\dataset\public\klifs_2026_07_22\raw`.
- Authoritative manifest:
  `D:\FORT\dataset\public\klifs_2026_07_22\processed\manifest.json`.
- Full processed registry:
  `D:\FORT\dataset\public\klifs_2026_07_22\processed\registry.parquet`.

**Representative schema metadata**

- The full registry footer reports 13,325 structure rows.
- It includes structure, kinase, species, PDB, chain, pocket quality,
  structural state, subpocket, kinase-family, canonical 85-residue pocket,
  ligand identity, scaffold, and interaction-fingerprint fields.
- The manifest defines an 85 aligned-residue by 7 interaction-type contract,
  or 595 interaction bits.
- The filtered MNI-0 registry footer reports 2,091 eligible complexes.

**Provenance and license**

- Source: KLIFS v3.2 API snapshot, site update 2026-07-22.
- The local manifest records SHA-256 values for every raw API object.
- The KLIFS FAQ states that KLIFS data are freely available for academia and
  industry, but no SPDX or Creative Commons identifier is recorded.
  Attribution is required, and redistribution of derived PDB or MOE content
  needs a separate terms review.

**Authorized role**

- Kinase pocket alignment, family/group taxonomy, structural-state metadata,
  and mechanism-supervision infrastructure.
- It is explicitly not an affinity-label source.

**Prohibited use**

- Do not interpret structures or interaction fingerprints as independent
  affinity measurements.
- Do not use KLIFS family labels as proof that a residue coordinate exceeds
  coarse taxonomy.
- Do not reuse structures, ligands, families, or templates that overlap a
  held downstream panel without fold-specific exclusion.
- Do not revive MNI-0 by increasing capacity; its registered target-
  conditioning audit is a fail-stop.

**Firewall and consumption state**

- The existing Reinecke firewall excludes all Reinecke development kinase
  families and exact connectivity/scaffold/near-neighbor overlaps before
  mechanism scoring.
- `reinecke_affinity_labels_read=false` in the KLIFS audit.
- `confirmation_labels_read=false` and `sealed_test_consumed=false`.

### 3.4 BindingDB native-article archive (`open_s`)

**Physical content**

- Exact root: `D:\FORT\dataset\public\open_s`.
- One file is present:
  `D:\FORT\dataset\public\open_s\BindingDB_BindingDB_Articles_202607_tsv.zip`.
- ZIP size: 18,114,757 bytes.
- ZIP central-directory metadata reports one member,
  `BindingDB_BindingDB_Articles.tsv`, with an uncompressed size of
  328,109,536 bytes.
- Source record:
  `D:\FORT\manifests\open_sources.json`.

**Representative schema metadata**

- Release: BindingDB 2026-07 native-article subset.
- The archive is intended to contain BindingDB staff-curated primary-article
  records rather than records imported from ChEMBL.
- No extracted TSV or derived affinity registry exists in this directory.

**Provenance and license**

- The frozen download URL, release, byte size, SHA-256, and timestamp are
  recorded in `D:\FORT\manifests\open_sources.json`.
- BindingDB data are recorded as CC BY 3.0 US. Terms must be verified at the
  BindingDB source before redistribution.

**Authorized role**

- Metadata, DOI, source-lineage, and candidate evidence recovery for the
  conditional OpenMut route.
- A future blind projection may be built only after the current acquisition
  and provenance gate passes.

**Prohibited use**

- Do not train on or inspect affinity rows under the current gate.
- Do not pool measurements across articles to manufacture query depth.
- Do not treat a BindingDB copy of a ChEMBL or publication measurement as an
  independent source.
- Do not infer exact endpoint, protein construct, mutation, or assay
  comparability from archive membership alone.

**Firewall and consumption state**

- Required closure axes are DOI/experiment lineage, database import lineage,
  exact endpoint, target construct and mutation, ligand parent/scaffold/
  neighbor, assay, document, and provenance.
- The active ledger records the affinity rows as unopened.
- This metadata inventory did not extract the archive.

### 3.5 OpenBind EV-A71 2A

**Physical content**

- Exact root: `D:\FORT\dataset\public\openbind_ev_a71_2026`.
- Formats: 18 CSV, 2 Parquet, 2 TSV, 5 Markdown, 1 Apache license, 7 Python,
  1 shell, 1 YAML, 5 PNG, and 2 small text files.
- `D:\FORT\dataset\public\openbind_ev_a71_2026\raw` is empty.
- The material present is a benchmark repository snapshot under
  `D:\FORT\dataset\public\openbind_ev_a71_2026\benchmark`.

**Representative schema metadata**

- The source README describes 925 crystallographic binding events, 699
  compounds, and 601 compounds with affinity measurements for one EV-A71 2A
  protease campaign.
- `all_affinity_data_release_v1.csv` is documented as row-level measurements;
  the reference table is one row per Fragalysis structure.
- The docking and cofolding Parquet footers report 277,500 and 158,675 pose
  rows. Their fields describe methods, ranks, pose scores, RMSD, LDDT-PLI,
  validity, ligand structure, and similarity. These pose rows are benchmark
  outputs, not independent binding measurements.

**Provenance and license**

- Data DOI: `10.5281/zenodo.20026661`.
- Data license: CC0 1.0 Universal.
- Benchmark repository code and documentation: Apache 2.0.

**Authorized role**

- Single-target local-SAR, support-label, pose, and structure engineering
  audits.
- The historical O0 result is evidence that strict local chemical support can
  help on this one protein.

**Prohibited use**

- Do not use this source to claim unseen-protein, multi-family, or strict
  target-cold transfer.
- Do not count many structures, poses, predictions, or repeated structures as
  independent proteins or measurement sites.
- Do not use benchmark prediction files as experimental supervision.
- Do not promote the source to an independent confirmation panel for a model
  designed after inspecting it.

**Firewall and consumption state**

- The O0 chemical firewall used parent identity, Murcko scaffold, and maximum
  Morgan Tanimoto below 0.50 between external training and held chemical
  series.
- Affinity labels were read in the historical O0 development audit. The source
  has no sealed target-conditioned confirmation role.
- Reopening requires independent targets or families, not another model on
  this campaign.

### 3.6 Papyrus 05.7

**Physical content**

- Exact root: `D:\FORT\dataset\public\papyrus_05_7`.
- The local snapshot contains only 3 of the 14 files described by the Papyrus
  README:
  - `D:\FORT\dataset\public\papyrus_05_7\raw\05.7++_combined_set_without_stereochemistry.tsv.xz`
    (56,759,540 bytes);
  - `D:\FORT\dataset\public\papyrus_05_7\raw\05.7_combined_set_protein_targets.tsv.xz`
    (1,780,032 bytes); and
  - `D:\FORT\dataset\public\papyrus_05_7\raw\README.txt`
    (12,813 bytes).
- No local processed registry is present.

**Representative schema metadata**

- The README documents ligand identity and standardized structure, protein
  target/accession/sequence, assay and document identifiers, endpoint flags,
  activity class, censoring relation, aggregate pChEMBL values, dispersion,
  and replicate-count fields.
- Target metadata include Papyrus target ID, source target ID, UniProt,
  reviewed status, organism, classification, sequence length, and sequence.
- The activity table is already aggregated at compound-protein level.

**Provenance and license**

- Version 5.7, README modified 2024-10-24.
- Papyrus aggregates ChEMBL, Guide to Pharmacology, and other sources.
- The README says the license is in `05.7_additional_files.zip`, but that
  archive and `LICENSE.txt` are absent locally. The local redistribution
  license is therefore **not verified**.

**Authorized role**

- Historical source/provenance and topology audit only.
- Target and sequence metadata may inform a future evidence registry after
  source lineage is closed.

**Prohibited use**

- No Papyrus pretraining, DTA training, or architecture comparison is
  authorized from this local corpus.
- Do not split semicolon-combined values or document lists into
  pseudo-observations.
- Do not treat one aggregated compound-target row as independent replication.
- Do not use it as independent confirmation of ChEMBL-derived training.

**Firewall and consumption state**

- The F0 audit required exact endpoint separation, unambiguous single-source
  and single-document provenance, human wild-type targets, and atomic row
  metadata.
- It found no document-replicated parent-target cells after strict resolution;
  the route stopped for provenance, not model capacity.
- Activity data were read in the historical provenance audit. No FORT
  development, confirmation, or sealed labels were consumed by that audit.

### 3.7 PLINDER 2024-06/v2

**Physical content**

- Exact root: `D:\FORT\dataset\public\plinder_2024_06_v2`.
- Formats: 3 Parquet files, 3 JSON files, and 2 NPZ files.
- The raw snapshot contains only:
  - `D:\FORT\dataset\public\plinder_2024_06_v2\raw\annotation_table.parquet`
    (1,018,256,581 bytes); and
  - `D:\FORT\dataset\public\plinder_2024_06_v2\raw\split.parquet`
    (10,651,188 bytes).
- No local PLINDER coordinate-system tree is present.
- Derived dual-cold objects are under
  `D:\FORT\dataset\public\plinder_2024_06_v2\processed\dualcold`.

**Representative schema metadata**

- The annotation-table footer reports 1,357,906 rows and extensive entry,
  system, pocket, ligand, structure-validation, affinity-availability,
  sequence-similarity, and chemical-similarity metadata.
- The split-table footer reports 409,726 rows with system, uniqueness, split,
  cluster, validation, ligand/pocket, affinity-availability, and apo/predicted
  structure fields.
- The processed registry footer reports 9,128 target-unit by ligand-parent
  rows. It carries target unit, cluster, accession, ligand parent, scaffold,
  affinity, reliability, source date/system, chemical similarity, and frozen
  split fields.
- Derived NPZ files contain ligand and pooled target features; JSON files
  contain sequence maps and the manifest.

**Provenance and license**

- Raw object URLs, generations, ETags, SHA-256 values, and byte sizes are
  frozen in `D:\FORT\manifests\raw_datasets.json`.
- License is recorded as CC BY 4.0 with attribution verification required
  before redistribution.
- The processed manifest describes affinity as PLINDER
  `ligand_binding_affinity`, BindingMOAD-derived, and treated as exact by the
  historical registry.

**Authorized role**

- Historical structure-affinity engineering and native-pose diagnostic
  substrate.
- Metadata-only structure, topology, and firewall planning under the current
  stop.

**Prohibited use**

- Do not interpret the local metadata snapshot as a full PLINDER structure
  installation.
- Do not perform new label training or confirmation scoring on this historical
  registry under the current information gate.
- Do not assume that target/ligand novelty establishes independent assay,
  document, or provenance.
- Do not use sparse per-target coverage as a powered factorial
  target-ligand-reordering test.

**Firewall and consumption state**

- The processed manifest records separation by target unit, cluster,
  accession, system, ligand parent, scaffold, and chemical similarity.
- `affinity_labels_read=true`: historical engineering audits consumed
  non-test affinity labels.
- PLINDER test components were excluded before affinity scanning.
  `confirmation_labels_read=false` and `sealed_test_consumed=false`.

### 3.8 Reinecke 2024

**Physical content**

- Exact root: `D:\FORT\dataset\public\reinecke_2024`.
- Formats: 5 XLSX workbooks, 1 checksum file, 2 JSON files, 1 Parquet registry,
  and 2 NPZ feature files.
- Raw supplementary workbooks and hashes:
  `D:\FORT\dataset\public\reinecke_2024\raw`.
- Processed panel:
  `D:\FORT\dataset\public\reinecke_2024\processed\panel_reinecke`.

**Representative schema metadata**

- Workbook metadata identifies compound annotation, target annotation,
  Kinobeads drug matrices, selectivity tables, orthogonal validation, and
  unrelated phosphosite sheets.
- The processed Parquet footer reports 9,333 rows with gene/accession, ligand
  parent, apparent affinity, reliability, scaffold, compound identifiers,
  homology cluster, external-novelty checks, chemical similarity, endpoint,
  document, target, and frozen split.
- The frozen development block contains 826 cells, 109 targets, 104 homology
  components, 171 ligands, and 121 scaffolds.
- Endpoint: `pKd_app = 9 - log10(Kd_app,nM)`.

**Provenance and license**

- Source: Reinecke et al., Nature Chemical Biology 2024,
  DOI `10.1038/s41589-023-01459-3`.
- Article and supplementary data are recorded as CC BY 4.0.
- SHA-256 values for all five source workbooks are frozen locally.

**Authorized role**

- Development-only kinase panel for powered, preregistered comparisons.
- It is useful for apparent-Kd ranking and reliability diagnostics within its
  stated measurement system.

**Prohibited use**

- Do not claim independent confirmation; aggregate label shape and resolution
  were inspected before registry construction.
- Do not silently equate pKd_app with biochemical Ki/Kd from another source.
- Do not tune a candidate to this panel after observing historical failures.
- Do not treat its sparse measurement graph as a complete factorial design.

**Firewall and consumption state**

- Evaluation target homology components are novel relative to historical
  ChEMBL/Metz; query ligand parents and scaffolds are novel relative to
  historical substrates; exact Morgan similarity is below 0.95; held target
  components are removed at fitting time.
- Historical development labels were consumed in registered audits.
- `confirmation_labels_read=false` and `sealed_test_consumed=false`.

### 3.9 Novartis SPD 2023

**Physical content**

- Exact root: `D:\FORT\dataset\public\spd_2023`.
- Formats: 5 tab-delimited text files, 1 XLSX workbook, and 1 JSON manifest.
- Raw source files:
  `D:\FORT\dataset\public\spd_2023\raw`.
- Authoritative source manifest:
  `D:\FORT\dataset\public\spd_2023\manifest.json`.

**Representative schema metadata**

- The activity header contains DrugCentral structure ID, InChIKey, assay and
  assay-group identifiers, censoring prefix, summarized IC50, concentration-
  response counts, source-specific AC50 summaries, single-concentration
  summaries, mechanism-of-action flag, gene mapping, exposure margins, and row
  ID.
- Separate files map assay groups to genes, assays to thresholds, codes to
  names, and preferred assay annotations.
- The frozen manifest reports 1,948 compounds, 144 gene-mapped assay groups,
  101 genes across multiple protein families, and 87,757 gene-mapped tested
  cells. These are prior manifest statistics, not values read by this
  inventory.
- The endpoint includes summarized IC50/AC50 with explicit censoring and many
  retained tested-negative cells.

**Provenance and license**

- Citation: Brennan et al., Nature Communications 2023,
  DOI `10.1038/s41467-023-40064-9`.
- Data DOI: `10.5281/zenodo.8103950`.
- Data license: CC BY 4.0. Source code license: MIT.
- Source URLs, download date, MD5, and SHA-256 are frozen in the local
  manifest.

**Authorized role**

- Source-shape, censoring, multi-family, and pseudo-prospective design audit.
- Template for the required systematic inactive-retaining prospective
  substrate.

**Prohibited use**

- Do not convert censored assay-conditioned IC50/AC50 into exact universal pK
  labels.
- Do not run the planned full dual-cold model: the prior power audit found
  inadequate query depth and severe inactive-floor dominance.
- Do not use SPD as independent confirmation of ChEMBL-trained models without
  closing heavy marketed-drug and source overlap.
- Do not ignore assay group, source contribution, or censoring in any future
  likelihood.

**Firewall and consumption state**

- Any self-contained future use requires target-homology, ligand-parent,
  scaffold, chemical-neighbor, assay, gene, source, endpoint, and provenance
  closure.
- Activity records were read in the historical source and power audit, not in
  this inventory. That audit did not consume FORT development, confirmation,
  or sealed labels.
- Terminal source verdict: systematic and valuable, but not dual-cold-capable
  at the registered power target.

### 3.10 ToxCast invitrodb v4.3 projection

**Physical content**

- Exact root: `D:\FORT\dataset\public\toxcast_invitrodb_v4_3`.
- Formats: 2 XLSX workbooks, 1 GZIP CSV projection, and 1 JSON manifest.
- Official annotation workbooks:
  `D:\FORT\dataset\public\toxcast_invitrodb_v4_3\raw`.
- Lightweight projection:
  `D:\FORT\dataset\public\toxcast_invitrodb_v4_3\raw\toxcast_data_deepchem_projection.csv.gz`.

**Representative schema metadata**

- The assay-annotation workbook contains `annotations_combined`, `assay`,
  `assay_component`, and `assay_component_endpoint` sheets.
- The target-mapping workbook contains one mapping sheet.
- The projection header contains `smiles` followed by many assay-response
  columns. It is a wide multi-task projection, not a continuous affinity
  registry.
- The full official invitrodb activity release is not present locally.

**Provenance and license**

- Official source: US EPA ToxCast invitrodb version 4.3,
  DOI `10.23645/epacomptox.6062623.v14`, recorded release date 2025-09-03.
- Official release license: CC0.
- The activity projection is the DeepChem MoleculeNet ToxCast CSV. Its labels
  originate from public EPA ToxCast and are not a new independent measurement
  source.
- File hashes and official source identifiers are frozen in
  `D:\FORT\dataset\public\toxcast_invitrodb_v4_3\manifest.json`.

**Authorized role**

- Lightweight target-contrast or observation-head pretext engineering.
- Assay/target mapping and source-closure planning.
- Full-scale work, if ever authorized, must return to the official invitrodb
  release rather than treating this projection as complete.

**Prohibited use**

- Do not use the projection as a continuous DTA affinity substrate.
- Do not interpret assay-response tasks as direct protein-ligand binding
  measurements without assay-mechanism review.
- Do not count the DeepChem projection and EPA release as independent
  sources.
- Do not use it for confirmation of a model designed on overlapping public
  chemistry.

**Firewall and consumption state**

- A future fold must close ligand identity/scaffold/chemical neighbor,
  assay-component and target mapping, source lineage, endpoint, and downstream
  target/homology/provenance overlap.
- No current ToxCast label use is recorded in the active research ledger.
- The source manifest records `confirmation_labels_read=false` and
  `sealed_test_consumed=false`.

## 4. Cross-source interpretation

The local tree is a layered evidence system, not a single training corpus:

| Evidence layer | Local examples | What it can establish | What it cannot establish |
| --- | --- | --- | --- |
| Historical exact or apparent affinity | ChEMBL, Metz, Davis, Reinecke, PLINDER | Engineering behavior and within-source rankings under frozen protocols | New independent multi-family confirmation |
| Mutation or ordinal response | KirHub | Local mutation-dependent ligand reordering | Exact affinity or cross-source generalization |
| Pocket and structure metadata | KLIFS, PLINDER, OpenBind | Alignment, structure quality, pose and mechanism audits | Independent target-ligand measurements |
| Aggregated open activity | Papyrus, BindingDB archive | Candidate discovery and provenance recovery | Replication until DOI/experiment lineage is closed |
| Systematic/censored screening | SPD, ToxCast projection | Inactive-retaining and multi-task observation design | Powered continuous dual-cold affinity prediction |

No source currently satisfies all of the required substrate conditions at once:
multiple independent protein families, adequate scaffold-diverse query depth,
complete or known-probability target-ligand measurement, endpoint consistency,
independent provenance, and adequate power for target-specific reordering.

## 5. Mandatory use rules

1. File presence does not authorize label access or model training.
2. Every future data projection must be built independently inside each
   training fold after downstream target, homology, binding-profile, ligand,
   scaffold, chemical-neighbor, assay, document, and provenance closure.
3. Endpoint and censoring semantics must remain explicit. Exact Ki, exact Kd,
   apparent Kd, IC50/AC50, single-dose residual activity, and binary assay
   responses are not interchangeable labels.
4. Database mirrors and records derived from the same DOI, experiment, or
   deposited measurement count as one evidence source.
5. Development access permanently prevents an independent-confirmation claim
   for that material.
6. Davis remains a single-use protected asset for target-conditioned
   confirmation. The program-level historical ChEMBL confirmation flag and
   arm-blind Davis power access must not be reset by a local manifest.
7. Raw-source redistribution is prohibited until each source-specific license
   and attribution requirement is verified. KirHub and Papyrus are the clearest
   unresolved local cases.
8. The present research stop remains in force: real-label pretraining or
   complex-model rescue requires the active information, topology, provenance,
   and power gates to pass first.

## 6. Local authority records

The following files are the principal local authorities for this inventory:

- `D:\FORT\manifests\raw_datasets.json`
- `D:\FORT\manifests\open_sources.json`
- `D:\FORT\dataset\public\chembl_37\processed\dualcold\manifest.json`
- `D:\FORT\dataset\public\chembl_37\processed\panel_metz\manifest.json`
- `D:\FORT\dataset\public\chembl_37\processed\panel_davis\manifest.json`
- `D:\FORT\dataset\public\klifs_2026_07_22\processed\manifest.json`
- `D:\FORT\dataset\public\plinder_2024_06_v2\processed\dualcold\manifest.json`
- `D:\FORT\dataset\public\reinecke_2024\processed\panel_reinecke\manifest.json`
- `D:\FORT\dataset\public\spd_2023\manifest.json`
- `D:\FORT\dataset\public\toxcast_invitrodb_v4_3\manifest.json`
- `D:\FORT\reports\active\kirhub_dd_decision.md`
- `D:\FORT\reports\active\kirhub_spkop_a1_strict_decision.md`
- `D:\FORT\reports\active\kirhub_pocket_oracle_c2_decision.md`
- `D:\FORT\reports\active\openbind_o0_decision.md`
- `D:\FORT\reports\active\papyrus_f0_decision.md`
- `D:\FORT\reports\active\reinecke2024_power_decision.md`
- `D:\FORT\reports\active\ramci_s0_audit.md`
- `D:\FORT\history.md`
- `D:\FORT\task.md`

When a dataset-local manifest conflicts with the program-level consumption
ledger, the stricter program-level state controls.
