# Dataset Failure Cleanup Record

**Date:** 2026-08-03
**Scope:** one quarantined failed archive only. Raw inputs, verified public
archives, snapshots, processed datasets, formal-training packages, and all
processing-history run records are retained.

## Verified Failed Artifact

| Field | Value |
|---|---|
| Path | `dataset/public/chembl_historical/archives/chembl_24_1_sqlite.tar.gz.sha256-46937b804020de62714fd7791de1a7b2d303652c094009e1b2bebfdf42681922.corrupt` |
| Size | 3,659,492,620 bytes |
| Observed SHA-256 | `46937b804020de62714fd7791de1a7b2d303652c094009e1b2bebfdf42681922` |
| Official SHA-256 | `6bb1030408c68b26ad8e9e9ae34ca7226e55d4214a6a057d93be987bbba5ea8c` |
| Historical status | `FAILED_HASH_MISMATCH`, then `QUARANTINED` |
| Evidence | `configs/a2s_dataset_processing/retrospective_events.v1.json`, events `chembl241_corrupt_concurrent` and `chembl241_quarantine`; `task.md` records the same STOP. |

The verified archive at
`dataset/public/chembl_historical/archives/chembl_24_1_sqlite.tar.gz` was
re-hashed before cleanup and matches the official SHA-256. The extracted
snapshot remains untouched.

## Retained Historical Records

The following run directories are retained. Their names contain `partial` or
`corrupt`, but they are processing-provenance records, not failed data folders;
two are explicitly `SUCCESS` and the others document the failed-transfer event.

- `20260731T125140855280Z_retro-chembl241-corrupt_63e63957`
- `20260731T125146802293Z_retro-chembl27-initial-partial_4c3cb40b`
- `20260731T125149797449Z_retro-chembl31-initial-partial_d1c00275`
- `20260731T125206702540Z_retro-prelaunch-delete-corrupt_e78ef5a3`
- `20260731T125242959209Z_relocate-chembl31-prelogging-partial_81d863c1`
- `20260731T125335097234Z_stage-chembl31-prelogging-partial_9574e76d`

## Result

The quarantined failed file was deleted after the checks above. Expected space
reclaimed: 3,659,492,620 bytes. No other dataset object was deleted in this
cleanup.
