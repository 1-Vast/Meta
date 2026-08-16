# OpenMut `OMUT-D0` decision

**Verdict:** `OMUT_D0_SOURCE_FREEZE_COMPLETE`

**Date:** 2026-07-28.
**Preregistration:** `reports/active/omut_d0_preregistration.md`,
SHA-256 `01d3f93ac65df3e3541b48c1b2db6865c9c3ccdd62edf533ada943ba56ca9b7a`
(amendments A1-A3, all frozen before the accepted run).
**Runner:** `research/omut_d0.py`. **Result:** `reports/active/omut_d0.json`.
**Registry:** `manifests/omut_d0_source_registry.v1.json`.
**Environment:** `D:\anaconda\envs\drug\python.exe`, torch 2.6.0+cu124, CUDA 12.4,
RTX 4060 Laptop GPU. No GPU computation was performed; D0 is a metadata freeze and
the environment is recorded so later stages inherit a verified binding.

## 1. Gates

| gate | result |
| --- | --- |
| `D0_PROBE_COMPLETE` | pass |
| `D0_VERSION_FROZEN` | pass |
| `D0_RIGHTS_RESOLVED` | pass |
| `D0_SCHEMA_FROZEN` | pass |
| `D0_CHECKSUM_FROZEN` | pass |
| `D0_ACQUISITION_PLAN` | pass |
| `D0_NO_VALUES_READ` | pass |
| `D0_REPRODUCIBLE` | pass |

Reproducibility anchor `substantive_registry_sha256`
`1953cff4c1d7301c51d1ef934c0c5c913f7c154022ad515c53e416bbce8f82f9`.
`research/omut_d0.py --verify` re-executed the entire freeze against the live
services and agreed with the stored artefact on every field
(`differing_field_count = 0`).

`registry_sha256` is `9267cba41fa3880dc244a7d8618509aa3ec5f176049b37a759d68301d994dc4d`.
It covers live transport sizes and is *not* expected to be stable: repository
metadata payloads change size as counters change. Cite the substantive anchor.

## 2. Firewall

25 reads, 1,172,945 bytes materialised, **0 rows materialised**, **0 violations**.

Exactly two reads touched an affinity-bearing object, both header-only and both
capped:

| source | mode | bytes | rows |
| --- | --- | ---: | ---: |
| `bindingdb_local_archive` (`BindingDB_BindingDB_Articles.tsv`) | `header_only` | 65,536 | 0 |
| `platinum_flatfile` (`platinum_flat_file.csv`) | `header_only` | 65,536 | 0 |

`davis_complete.tab` was never requested; it is name-guarded in the transport
layer, so the runner cannot fetch it under any argument. The two DAVIS-Complete
FASTA members carry sequence only and were streamed in memory, never written to
disk; both reproduced their Dataverse-published MD5 exactly
(`48dd9b67...`, `8d49c499...`), 444 records each.

The firewall aborted the first execution attempt on
`data.latestVersion.files[].dataFile.checksum.value`. That attempt wrote nothing.
Rather than extend the allowlist one abort at a time, a non-aborting sweep
enumerated every key path in every payload where a frozen token carries a value;
amendment A2 then replaced the guessed allowlist with the enumerated set. The
complete set is four paths, each an identifier and none a measurement: Dataverse
file digests, and ORCID author identifiers in the two Europe PMC records.

## 3. What is now frozen

| source | status | rights | derived use | plan |
| --- | ---: | --- | --- | --- |
| `davis_complete_dataverse` | 200 | CC0-1.0 | yes | `deferred_stage` |
| `davis_complete_github` | 200 | none declared | no | `safe_metadata_only` |
| `davis_complete_zenodo` | 200 | CC-BY-4.0 | yes | `deferred_stage` |
| `bindingdb_download_index` | 206 | CC-BY-3.0-US | yes | `safe_metadata_only` |
| `bindingdb_local_archive` | local | CC-BY-3.0-US | yes | `already_local_verified` |
| `chembl_status` | 200 | CC-BY-SA-3.0 | yes | `safe_metadata_only` |
| `chembl_variant_layer` | 200 | CC-BY-SA-3.0 | yes | `bounded_api_projection` |
| `chembl_assay_variant_link` | 200 | CC-BY-SA-3.0 | yes | `bounded_api_projection` |
| `chembl_activity_variant_projection` | 200 | CC-BY-SA-3.0 | yes | `bounded_api_projection` |
| `platinum_flatfile` | 206 | none declared | **no** | `blocked_rights` |
| `platinum_rights` | 200 | none declared | **no** | `blocked_rights` |
| `mdrdb_rights` | 200 | academic only | yes (academic) | `safe_metadata_only` |
| `binder2030` | 200 | unresolved | unknown | `blocked_unresolved` |
| `proteingym` | 200 | unresolved | unknown | `blocked_unresolved` |
| `mavedb` | 200 | unresolved | unknown | `blocked_unresolved` |
| `europepmc_supplements` | 200 | unresolved | unknown | `deferred_stage` |

Versions frozen: DAVIS-Complete Dataverse **v3.0**, released 2025-11-27, four files
with published MD5s; GitHub head `799ac5696e7afb9a23d2767cace2352e243b353e`; Zenodo
record 15391611 (CC BY 4.0, 39.3 GB + 2.8 GB archives, not fetched); BindingDB
**202607** Articles, local copy SHA-256 verified, single TSV member of 328,109,536
bytes with **640 columns**; **ChEMBL_37**, released 2026-05-01; PLATINUM flat file
687,373 bytes with **51 columns**; MdrDB download page as retrieved; Binder2030
`10.1016/j.slasd.2026.100299`; ProteinGym **PG_v1.3**; MaveDB API **2026.2.7**;
Europe PMC REST **6.9**.

Every count below is an unverified vendor or API claim. None is an independent
component count, and none has survived a provenance collapse.

## 4. The three findings that matter for the programme

**4.1 The ChEMBL variant layer is small enough to acquire, and the acquisition is
now bounded.** The layer has no dedicated resource — `/variant_sequence/schema.json`
returns 404 — and is exposed as a nested related schema on `assay` with seven
fields (`accession`, `isoform`, `mutation`, `organism`, `sequence`, `tax_id`,
`version`). Server-side filters resolve the slice: **20,150 assays** carry a
`variant_sequence`, and **119,801 activities** carry an `assay_variant_mutation`,
against 24,527,044 activities in ChEMBL_37. That is 0.49% of the table, so no full
database dump is required and the `OMUT-D0` blocker "the complete ChEMBL-37 variant
layer or a bounded API acquisition" is resolved in favour of the bounded route.
The count was obtained under a blind projector: the server returned `activity_id`
only, dropping 46 fields including every numeric affinity, relation, unit and censor
column, and the runner asserts that no unrequested field came back. ChEMBL's own
warning stands — variant sequences are not referentially linked to component
sequences, so engineered versus disease variants need document review at `OMUT-X0`.

**4.2 PLATINUM has exactly the right schema and no licence, so rights are the
binding constraint rather than a formality.** Its 51 columns include the paired
mutant contrast directly — `affin.k_wt`, `affin.k_mt`, `affin.delta_k`,
`affin.fold_change`, `affin.unit` — with construct fields `mutation`, `mut.uniprot`,
`mut.is_single_point`, `mut.wt_pdb`, `mut.mt_pdb`, and document-level provenance
`mut.doi`, `mut.pmid`, `affin.exptal_method`, `affin.temperature`, `affin.ph`. This
is the shape the programme has been unable to find in open data. The data page
carries an Open Knowledge badge and states no licence identifier, no terms of use
and no redistribution grant. Availability is not a licence, so the source stays
`blocked_rights`: rights and schema inspection only, no acquisition, no gold subset,
no labels. Resolving PLATINUM's terms is now the single highest-value non-scientific
action available to this programme.

**4.3 No newly resolved source changes the powered-substrate verdict.** Binder2030
is indexed but its dataset is unreachable through the open route: Europe PMC reports
the article as not open access with no supplementary files indexed, the version of
record is CC BY-NC-ND (an article licence, not a data or model-training licence), and
eight fields remain blocking — including whether any WT-to-single-substitution
construct pair exists in it at all. MdrDB grants academic use but aggregates GDSC,
DepMap, AIMMS, KinaseMD, PLATINUM, TKI and RET, so it is not an independent
provenance lineage and re-exports PLATINUM; its cell-response, predicted-structure,
docking and simulated outcomes remain permanently excluded. ProteinGym and MaveDB
carry software licences (MIT, AGPL-3.0) that say nothing about the aggregated assay
data, and their fitness, expression, stability and function outcomes are not affinity
labels under any circumstance.

## 5. Claim boundary

D0 froze sources. It is **not** evidence about substitution semantics, about
target-specific ligand reordering, about topology, about power, or about strict
dual-cold transfer. The programme verdict is unchanged:

> **3** - current data do not identify the substitution-geometry or tau-teacher
> mechanism; new prospective measurement conditions or a newly recovered,
> source-resolved public substrate are required.

Nothing here authorises reading an affinity value, downloading an affinity-bearing
table, materialising topology or effective rank, fitting a coordinate or estimator,
pretraining on open data, or reclassifying `panel_davis` or any sealed asset.
Declared record counts are vendor claims. A row count is not a component count; a
database is not a provenance unit; one DOI remains one provenance unit however many
databases re-export it.

## 6. What this unlocks, and what blocks next

`OMUT-D0` is complete. It unlocks **`OMUT-F0` and nothing else**.

`OMUT-F0` cannot be executed by a runner. It requires the irreversible Davis role
decision recorded in `task.md`, which is a human choice:

1. **preserve** `panel_davis` as single-use confirmation, which excludes every
   overlapping WT value and therefore generally prevents the WT-mutant contrast that
   DAVIS-Complete exists to supply; or
2. **retire** `panel_davis` from confirmation before the read, after which all
   Davis-derived records are development-only forever.

The choice is one-way and silent reassignment is prohibited. Until it is made and
recorded, `OMUT-D1` (label-free topology and projected power) and every later stage
stay blocked, and DAVIS-Complete affinity access is not permitted.

Independently of that choice, two actions are now well defined and unblocked:

- resolve PLATINUM's usage and redistribution terms (see 4.2), which is the only
  route by which an open, document-resolved, WT/mutant paired-affinity schema
  becomes admissible at all;
- the bounded ChEMBL variant acquisition (see 4.1) is specified and sized, and
  becomes executable at `OMUT-X0` under the existing DOI-level collapse and
  endpoint-stratification rules.
