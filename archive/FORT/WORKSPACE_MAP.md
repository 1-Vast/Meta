# FORT Workspace Map

**Inventory date:** 2026-08-03

This is a non-destructive navigation index for `D:\FORT`. It records the
current layout and retention policy; it does not change scientific decisions,
dataset roles, or Git history.

## Start Here

| Need | Canonical location | Rule |
|---|---|---|
| Current work authorization | `task.md` | Treat a later explicit correction or STOP decision as superseding an earlier task entry. |
| Experiment ledger | `history.md` | Record outcomes, including failures, before reopening a branch. |
| Data provenance summary | `DATASET_RECORD_SUMMARY.md` | Do not infer training permission from the presence of a dataset. |
| Runnable entry points | `main.py`, `scripts/` | Read the command contract before running a research module. |
| Research-report navigation | `reports/active/CURRENT_INDEX.md` | `reports/active/` is evidence storage, not an authorization queue. |

## Directory Roles

| Path | Contents | Retention and handling |
|---|---|---|
| `configs/` | Versioned experiment configuration. | Keep beside the code that consumes it. |
| `dataset/` | Raw, public, processed, formal-training, and innovation-test data. About 107.8 GiB. | Never bulk-move, deduplicate, or delete without a provenance migration and checksum audit. |
| `manifests/` | Contracts, source registries, and run plans. | Preserve paths and hashes; these are reproducibility records. |
| `model/` | Consolidated reusable model modules. | Do not add experimental branches here before a written gate passes. |
| `research/` | Isolated research runners and staging code. | Failed branches remain represented by reports/history, then can be removed only through an explicit cleanup decision. |
| `reports/active/` | Reports and machine-readable evidence. | Use `CURRENT_INDEX.md`; a file remains here for provenance even after a route is stopped. |
| `reports/archive/` | Explicitly archived historical snapshots. | Archive only complete, self-contained record groups with their manifests and hashes. |
| `scripts/` | Reusable preprocessing, audit, and training utilities. | Preserve public command-line contracts. |
| `tests/` | Regression and gate tests. | Keep tests adjacent to the active implementation surface. |
| `tmp/` | Local mirrors, third-party source, literature copies, and bounded scratch evidence. About 1.24 GiB. | See `tmp/README.md`; do not treat this as disposable wholesale. |
| `.git/` | Repository history and index. | Never manually clean or rewrite during file organization. |

## Data Layout

The dataset tree is already divided by role and should remain so:

| Path | Role |
|---|---|
| `dataset/public/` | Downloaded and historical public sources (about 107.3 GiB). |
| `dataset/raw/` | Original local input records. |
| `dataset/processing_history/` | Transformation provenance. |
| `dataset/processed/` | Derived datasets and audit outputs. |
| `dataset/formal_training/` | Formal training packages and rosters. |
| `dataset/innovation_tests/` | Narrow test substrates for registered gates. |
| `dataset/ready/` | Prepared material with a separate readiness contract. |

The large ChEMBL snapshots and archives under `dataset/public/chembl_historical/`
are deliberate provenance assets, not duplicate caches. They must remain in
place until a checksum-preserving storage migration is explicitly approved.

## Report Reading Order

1. Read `task.md` and the latest relevant entries in `history.md`.
2. Use `reports/active/CURRENT_INDEX.md` to identify the report family and any
   explicit corrections.
3. Read the machine-readable artifact named by that report before relying on a
   numerical conclusion.
4. Treat reports without a corresponding artifact, code path, or stated gate as
   hypotheses or context, not executable authority.

The PSEP family illustrates why this order matters: the early core decision is
preserved as evidence, while `PSEP_R3_CORRECTION_2026-08-02.md` explicitly
withdraws the representation-mechanism claim and records `NO CORE MECHANISM
IDENTIFIED`.

## Cleanup Boundary

Only the following directories are treated as regenerable local cache in this
organization pass:

- `__pycache__/`
- `.pytest_cache/`
- `model/__pycache__/`
- `research/__pycache__/`
- `scripts/__pycache__/`
- `tests/__pycache__/`
- `tests/core/__pycache__/`

No dataset, report artifact, manifest, checkpoint, third-party checkout, or
Git-tracked/deleted user work is included in this boundary. In particular,
`tmp/adambind-source/` may contain local compatibility work and is not a cache.

## External Review Register

The independent committee review supplied outside the repository is retained at:

`C:\Users\59964\Desktop\编程\1\A2S_DTA_INDEPENDENT_COMMITTEE_REVIEW_2026-08-01.md`

It is external review input. It should be cited from a decision record when
used, but is not copied into the repository automatically because no canonical
import/provenance contract has been established for external attachments.
