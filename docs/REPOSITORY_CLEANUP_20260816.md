# Repository Cleanup Decision, 2026-08-16

The repository is an executable research program, not a storage mirror. This
cleanup keeps current code, governed datasets and compact numerical evidence;
Git history carries superseded source and raw leaf artifacts.

## Keep

- active DTA data, training and evaluation code;
- the incumbent and current Pareto-frontier model implementations;
- R0-R14 reports, compact `RESULT.json` authorities and comparison tables;
- the record audit, split seals and current checkpoints needed for strict
  loading or the next A2-moment experiment;
- `history.md`, `task.md`, `CURRENT_MODEL_EVIDENCE.md` and
  `EVIDENCE_LEDGER.md` as the four authority documents.

## Remove from the working tree

- Python/test caches and temporary test directories;
- ignored checkpoints and compressed predictions for failed or superseded
  experiments;
- pre-R0 2026-08-12 through 2026-08-14 loose experiment directories after
  consolidation in `LEGACY_PRE_R0_SUMMARY.md`;
- retired `research/` prototypes and their dedicated tests;
- crossed-interaction/source-affinity/mechanism reports whose decisions are
  already represented in history;
- physical archive snapshots. Their recovery commits and meaning are kept in
  `archive/README.md`;
- closed code families that are neither active nor required to reproduce a
  current frontier arm.

## Recovery

- `c61cb8a`: complete retained R1-R13 experiment evidence tree before pruning.
- `40d68ba`: R1-R13 active model, trainer and gate surface before pruning.
- `a249a2c`: archive relocation before physical archive removal.
- earlier commits remain the authority for pre-R0 and predecessor material.

Deletion does not promote or erase a scientific claim. A compact result and
its terminal decision are retained; raw training payloads are not.

## Completed outcome

- physical archive mirrors removed; `archive/` now contains only its recovery index;
- `report/` reduced from about 2.18 GB to about 180 MB;
- checkpoints reduced from 296 historical payloads before consolidation to 18:
  A0/B3/C2 Pareto checkpoints and the R14 three-arm diagnostic set;
- `report/meta_fewshot/` reduced to formal R0-R14, independent M0 and one legacy summary;
- retired research packages, closed model families and their dedicated tests removed;
- active suite: 227 passed, 6 opt-in tests skipped in the `drug` environment;
- record audit: zero seal/hash/loading violations and 9 formal comparator checkpoints
  strictly reloaded (R14 diagnostics are retained separately).
