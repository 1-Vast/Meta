# OpenMut `OMUT-D0` preregistration

**Frozen:** 2026-07-28, before `research/omut_d0.py` existed and before any
`OMUT-D0` probe was executed.
**Environment:** `D:\anaconda\envs\drug\python.exe` (Python 3.11.15,
torch 2.6.0+cu124, CUDA 12.4, RTX 4060 Laptop GPU, capability 8.9).
**Parent design:**
`reports/active/openmut_delta2rank_feasibility_crossreview_2026-07-28.md`,
ordered execution decision, stage 1.
**Authority:** `task.md`, "Only active stage: `OMUT-D0`".

**Amendment A1 (2026-07-28, before `research/omut_d0.py` existed and before any
accepted or unaccepted formal run).** Exploratory metadata probing showed that
the section-3 field assertion, applied literally to raw payloads, fires on two
classes of object that contain no measurement:

1. **Resource-schema documents.** A Tastypie/ChEMBL `*/schema.json` document
   carries `standard_value`, `value`, `units`, and `relation` as *field names*
   (dictionary keys) whose payloads are field *definitions*. Scanning it as a
   record would abort the run on a document that by construction holds no data.
2. **Repository metadata blocks.** Harvard Dataverse encodes dataset citation
   metadata as `{"typeName": ..., "value": ...}` pairs, so the generic token
   `value` appears throughout a payload that contains only dataset metadata.

A1 therefore fixes three definitions and adds nothing else. Sources, modes,
gates, thresholds, vocabularies, verdicts, and stop rules are unchanged.

- **`rows_materialized`** counts rows read from an *affinity-bearing tabular
  payload* only. Metadata entries, file listings, schema field definitions, and
  header lines are never rows. Every ledger entry additionally records
  `payload_kind` in `{schema, metadata, header, digest, sequence, projection}`.
- **Schema payloads** are validated by a *structural* assertion instead of the
  record field assertion: the top-level keys must be a subset of the frozen
  Tastypie schema vocabulary, and every entry under `fields` must be a mapping
  of definition attributes. A schema document that returns anything else aborts
  the run. This is strictly stronger than the record assertion for this object.
- **Generic-token allowlist.** The generic tokens `value`, `relation`, `units`,
  and `text_value` are excused only at frozen key paths, declared in the runner
  and reproduced verbatim in `reports/active/omut_d0.json` under
  `firewall.generic_token_allowlist`. A generic-token hit at any other path,
  and any unambiguous-token hit (`standard_value`, `pchembl_value`, `ki`, `kd`,
  `ic50`, `ec50`, `affinity`, `potency`, `upper_value`, `activity_comment`) at
  any path whatsoever, aborts the run. The allowlist is auditable evidence, not
  a discretionary exemption: it is frozen before the run and printed in full.

**Amendment A2 (2026-07-28, after a first execution attempt aborted on the firewall and
before any result was written or accepted).** The A1 allowlist was written from a
partial reading of the payloads and contained two guessed paths that no payload
uses, while missing four real ones. A non-aborting sweep of every JSON payload the
runner touches was then run to enumerate the complete set, so that the allowlist is
derived from the payloads rather than extended one abort at a time. The first
attempt wrote no artefact; there is no invalidated result to retain.

A2 makes two changes and nothing else.

1. **The allowlist is replaced by the enumerated set.** The two unused guessed paths
   (`...metadataBlocks.*.fields[].value[].*.value` was retained, `...license.*` and
   `binder2030:...journalInfo.*` were removed) give way to the four paths that
   actually occur, each an identifier and none a measurement:
   Dataverse file digests at `data.latestVersion.files[].dataFile.checksum.value`,
   and ORCID author identifiers at `resultList.result[].authorList.author[].authorId.value`
   and `resultList.result[].authorIdList.authorId[].value` for both Europe PMC sources.
   The complete sweep is recorded in `reports/active/omut_d0_decision.md`.
2. **OpenAPI documents are schema documents.** The MaveDB payload is an OpenAPI
   description whose token hits are property *names* under `components.schemas`.
   The A1 structural assertion is generalised to take the allowed top-level key
   vocabulary and the definition-container key as parameters, and is applied to
   OpenAPI documents with the OpenAPI vocabulary. This closes an A1 gap in which
   that payload was neither record-scanned nor structurally asserted.

Sources, modes, gates, thresholds, rights and plan vocabularies, verdicts, and stop
rules remain unchanged. No affinity value was read at any point, including during the
aborted attempt and the sweep.

**Amendment A3 (2026-07-28, after a first written A2 artefact and before acceptance).**
The A2 run passed `D0_REPRODUCIBLE` while its `registry_sha256` nevertheless changed
between two executions. Both facts are correct: the frozen gate tests that
canonicalisation is deterministic, which it is, and the headline hash covers
`probe.content_length`, which is a live-service response size. GitHub repository
metadata grows and shrinks as counters change, so that field drifts for reasons that
have nothing to do with the freeze. A stage whose pass condition is "all sources
reproducible" must not publish an anchor that moves on a stargazer count.

A3 adds one reported quantity and changes no gate, threshold, source, mode,
vocabulary, verdict, or stop rule.

- `substantive_registry_sha256` is the canonical hash of the registry with
  `content_length` additionally stripped. It covers every version, rights
  determination, schema, checksum, declared count, and acquisition plan — that is,
  everything the freeze asserts — and excludes transport sizes, which it does not.
  This is the reproducibility anchor that later stages and audits must cite.
- `research/omut_d0.py --verify` re-executes the freeze and compares its
  `substantive_registry_sha256` against the stored artefact, reporting agreement or
  the differing fields. Its outcome is recorded in the decision as evidence.
- `registry_sha256` is retained and still reported, now explicitly labelled as
  including live transport sizes and therefore not expected to be stable.

`D0_REPRODUCIBLE` keeps its frozen definition and is not restated in terms of the new
quantity. Cross-run agreement is reported as evidence alongside the gate, not
substituted for it.

## 1. Question and claim boundary

D0 asks exactly one question:

> Can every candidate public source for a WT-to-single-substitution
> continuous-affinity substrate be frozen — version, URL, HTTP outcome,
> checksum, schema, rights, and a safe acquisition plan — reproducibly and
> **without reading a single affinity value**?

D0 is a *source-freeze* stage. It is not evidence about biology, about
substitution semantics, about ligand reordering, or about dual-cold transfer.

A D0 pass authorizes **only** entry into `OMUT-F0` (the one-way Davis role
decision). It does not authorize:

- downloading any affinity-bearing table;
- reading any affinity, relation, or censor value;
- topology, effective-rank, or power statements (those are `OMUT-D1`/`I0`);
- any coordinate, estimator, pretraining, or predictive claim;
- reclassifying `panel_davis` or any sealed asset.

Declared upstream record counts recovered from publisher or database metadata
are **unverified vendor claims**, never measurements, and never independent
component counts.

## 2. Registered sources

The registry is closed at freeze time. These and only these source IDs may
appear in the result:

| id | object |
| --- | --- |
| `davis_complete_dataverse` | Harvard Dataverse `doi:10.7910/DVN/RTQGP1` dataset record |
| `davis_complete_github` | `ZhiGroup/DAVIS-complete` repository and license state |
| `davis_complete_zenodo` | Zenodo record `15391611` |
| `bindingdb_download_index` | BindingDB public download endpoint for the 202607 Articles archive |
| `bindingdb_local_archive` | local `dataset/public/open_s/BindingDB_BindingDB_Articles_202607_tsv.zip` |
| `chembl_status` | ChEMBL web-services release record |
| `chembl_variant_layer` | ChEMBL `variant_sequence` resource |
| `chembl_assay_variant_link` | ChEMBL `assay` resource restricted to assays carrying a variant |
| `chembl_activity_variant_projection` | ChEMBL `activity` resource under a blind field projector |
| `platinum_flatfile` | PLATINUM mutation-affinity flat file |
| `platinum_rights` | PLATINUM usage/redistribution terms |
| `mdrdb_rights` | MdrDB academic and commercial terms |
| `binder2030` | Binder2030 dataset release |
| `proteingym` | ProteinGym reference release |
| `mavedb` | MaveDB API release |
| `europepmc_supplements` | Europe PMC supplementary-file service, as the source-bound supplement mechanism |

Adding a source after freeze invalidates the run. Removing one is recorded as
an explicit unreachable/blocked outcome, never as a silent deletion.

## 3. The no-value firewall

Every network or local read passes through one ledger. Each entry records
`source_id`, `url_or_path`, `mode`, `bytes_materialized`, `rows_materialized`,
`affinity_bearing`, and the SHA-256 of exactly the bytes that were
materialized.

Permitted modes, and nothing else:

| mode | allowed on | cap |
| --- | --- | --- |
| `metadata_json` | records with no affinity field | full response |
| `header_only` | affinity-bearing delimited files | first line only, `<= 65536` bytes |
| `zip_directory` | archives | central directory only, no member payload beyond a header line |
| `sequence_stream` | FASTA/sequence files with no affinity field | full stream, never written to disk |
| `local_digest` | already-present local files | bytes hashed, no field parsed |
| `blind_projection` | ChEMBL `activity` | server-side `only=` field restriction; the projector must drop every numeric affinity, relation, unit, and censor field before the row leaves the server |

Two mechanical assertions run on every materialized JSON record and on every
header line:

1. **field assertion** — if any key matches the frozen affinity vocabulary
   (`standard_value`, `standard_relation`, `standard_units`, `pchembl_value`,
   `value`, `relation`, `units`, `activity_comment`, `ki`, `kd`, `ic50`,
   `ec50`, `affinity`, `potency`, `upper_value`, `text_value`) **and** carries
   a non-null payload, the run aborts and the violation is recorded;
2. **row assertion** — for any source flagged `affinity_bearing`,
   `rows_materialized` must be exactly `0`.

`davis_complete.tab` is flagged affinity-bearing. Its content must never be
requested; only its Dataverse metadata (name, bytes, MD5) may be recorded.
The two DAVIS-Complete FASTA files carry sequence only and may be streamed in
memory to reproduce their published MD5. Nothing is written under
`dataset/**/raw`.

## 4. Gates

All eight are evaluated; all eight must pass.

| gate | condition |
| --- | --- |
| `D0_PROBE_COMPLETE` | every registered source has a probe record with a UTC timestamp and either an HTTP status or a local-file outcome |
| `D0_VERSION_FROZEN` | every source with `reachable == true` carries a non-null version/release identifier |
| `D0_RIGHTS_RESOLVED` | every source carries `rights.determination` from the closed vocabulary and an explicit `rights.permits_derived_research_use` in `{true, false, null}` |
| `D0_SCHEMA_FROZEN` | every source is either `schema_frozen == true` with hashed field names, or carries an explicit `schema_deferred_reason` |
| `D0_CHECKSUM_FROZEN` | every source with an acquirable file carries a published or computed checksum, or an explicit `checksum_absent_reason` |
| `D0_ACQUISITION_PLAN` | every source carries `acquisition.plan` from the closed vocabulary and a `blocking_fields` list |
| `D0_NO_VALUES_READ` | the ledger contains zero firewall violations and zero materialized rows from any affinity-bearing source |
| `D0_REPRODUCIBLE` | a second in-process construction of the canonical registry, with volatile fields (timestamps, latencies, network transcripts) removed, hashes identically |

`D0_RIGHTS_RESOLVED` requires that rights are *determined and recorded*, not
that they are permissive. A source whose license is absent is a recorded
`no_explicit_license` outcome with `permits_derived_research_use = false`;
that is a passing D0 record and a blocking downstream fact.

Closed vocabulary for `rights.determination`:
`explicit_open_license`, `explicit_restricted_license`, `academic_only`,
`no_explicit_license`, `unresolved`.

Closed vocabulary for `acquisition.plan`:
`already_local_verified`, `safe_metadata_only`, `bounded_api_projection`,
`blocked_rights`, `blocked_unresolved`, `deferred_stage`.

## 5. Verdicts

```text
OMUT_D0_SOURCE_FREEZE_COMPLETE     all eight gates pass
OMUT_D0_INCOMPLETE_STOP            any gate fails
OMUT_D0_FIREWALL_VIOLATION_ABORT   the no-value firewall fired
```

`OMUT_D0_SOURCE_FREEZE_COMPLETE` unlocks `OMUT-F0` and nothing else. `OMUT-F0`
additionally requires the irreversible Davis role decision recorded in
`task.md`; that decision is a human choice and is not made by this runner.

## 6. Stop rules

- A source that is unreachable is recorded as unreachable. It is not replaced
  by a mirror, a scraped copy, an aggregator re-export, or a private substitute.
- A source with unresolved rights stays `blocked_rights` /
  `blocked_unresolved`. Rights are not inferred from availability, from an
  open article license, or from another database's redistribution of the same
  records.
- A database name is never a provenance unit. BindingDB and ChEMBL copies of
  one paper are one provenance unit; that collapse is executed at `OMUT-X0`,
  and D0 only records the fields needed to execute it.
- D0 failure is not rescued by widening the source list, relaxing a
  vocabulary, or reclassifying an affinity-bearing file as label-free.
- No GPU computation is performed at D0. The CUDA environment is recorded so
  that later stages inherit an identical, verified binding; recording it is
  not a claim that any modelling occurred.

## 7. Deliverables

- `research/omut_d0.py` — the runner, deterministic apart from network transcripts.
- `tests/test_omut_d0.py` — firewall, vocabulary, gate, and reproducibility tests.
- `reports/active/omut_d0.json` — machine-readable frozen registry and gates.
- `manifests/omut_d0_source_registry.v1.json` — the canonical registry alone.
- `reports/active/omut_d0_decision.md` — verdict and claim boundary.
