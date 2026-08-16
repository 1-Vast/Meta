# Current Research Index

**Purpose:** navigate the active evidence ledger without mistaking preserved
artifacts for permission to continue a research branch.

## Authority Order

1. `task.md` defines current operational authorization.
2. `history.md` is the compact chronological experiment ledger.
3. A later report that explicitly says `supersedes`, `correction`, `withdrawn`,
   or `STOP` controls over the affected earlier report.
4. JSON, Parquet, NPZ, and model artifacts are the numerical evidence for the
   associated report; they do not independently authorize a new experiment.

When sources conflict, preserve both artifacts and record the conflict. Do not
silently select the most optimistic result or resume implementation.

## Dataset Integrity

| Question | Primary record | Current handling |
|---|---|---|
| Quarantined ChEMBL 24.1 acquisition failure | `DATASET_FAILURE_CLEANUP_2026-08-03.md`, `configs/a2s_dataset_processing/retrospective_events.v1.json` | The hash-mismatched `.corrupt` archive was removed after a fresh hash check. The verified archive, extracted snapshot, and processing provenance remain. |

## A2S / CMAL Evidence

| Question | Primary records | Current handling |
|---|---|---|
| Why CMAL did not establish a transferable operator | `CMAL_FAILURE_HANDOFF.md`, `A2S_DTA_INDEPENDENT_COMMITTEE_REVIEW_2026-08-01.md` (external), `a2s_cmal_v7_disjoint_source_meta300_seed1729.json` | Historical failure evidence; no CMAL optimization is authorized by these records. |
| Source information gate | `A2S_POST_REVIEW_V2_GATE_DECISION_2026-08-01.md`, `a2s_source_information_gate_v2_2026-08-01.json`, `a2s_source_information_gate_lock_v2_2026-08-01.json` | The recorded v2 decision is `NO_GO_INFORMATION_NOT_ADMITTED`. |
| Later A2S physical/state hypotheses | `A2S_NEA_PRECONDITIONS_D0_N0_N1_DECISION_2026-08-02.md`, `A2S_TRANSFER_OBJECT_GATE_T0_DECISION_2026-08-02.md`, and their named JSON/Parquet artifacts | Read with `history.md`; a stopped predecessor cannot be reopened by a proposal alone. |
| Mechanism-design context | `A2S_IDA_TAMSK_ADRO_DEEP_RESEARCH_SYNTHESIS_2026-08-01.md`, `A2S_FINAL_PI_META_MECHANISM_REDESIGN_2026-08-01.md` | Hypothesis/analysis only unless a later registered gate admits it. |

## PSEP Evidence and Corrections

| Record | Status for navigation |
|---|---|
| `PSEP_CORE_MECHANISM_DECISION_2026-08-02.md` | Preserves the D0/M0/M2/M3/M4 measurement and its proposed Stage-1 falsifier. It is not a license to build a mechanism. |
| `PSEP_OPERATOR_RESULTS_2026-08-02.md` | Registered result: `NO_OPERATOR_PASSES_ADMISSION_GATE`. |
| `PSEP_ABUNDANT_TO_SCARCE_DECISION_2026-08-02.md` | Registered result: `HEADROOM_IS_A_DOCUMENT_LOTTERY_NOT_TARGET_STRUCTURE`; support routing/selection is closed on this substrate. |
| `PSEP_REPRESENTATION_R1_R2_DECISION_2026-08-02.md` | Superseded only where the later R3 correction says so; preserve unmodified. |
| `PSEP_R3_CORRECTION_2026-08-02.md` | Latest PSEP correction. It withdraws the representation-mechanism claim, records `OBJECTIVE_EFFECT_UNRESOLVED_DUE_TO_CROSS_TASK_PAIRING`, and sets programme status to `NO CORE MECHANISM IDENTIFIED`. |

The only directed PSEP follow-up recorded by the correction is a role-safe,
within-unit/document/assay discovery gate. It must be separately authorized;
`validate` remains closed according to that document.

## Artifact Pairing Convention

For every report family, retain the report and the machine-readable companion in
the same directory. Common suffixes are:

- decision/proposal/preregistration: `.md`
- aggregate result: `.json`
- row or component evidence: `_records_*.parquet`
- numerical arrays: `.npz`
- checkpoints: `.pt`

Do not archive one member of a pair without the other. Do not change a result
filename to make it look current; add a correction or index entry instead.

## Archive Rule

`reports/archive/` is reserved for record groups that are both complete and
explicitly superseded. A future archive move must include a manifest listing all
members, relative paths, hashes, predecessor/successor documents, and the move
date. No such bulk move was made in this organization pass.
