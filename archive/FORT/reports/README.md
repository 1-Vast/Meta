# Report Layout

| Directory | Use |
|---|---|
| `active/` | Current A2S/CMAL/PSEP evidence, current indexes, and machine-readable companions. Start with `active/CURRENT_INDEX.md`. |
| `archive/legacy/` | Complete stopped historical branches, grouped by family. Start with `archive/legacy/README.md`. |
| `archive/` | Earlier immutable task/history snapshots and integrity hashes. |

Reports are retained as evidence even when an experiment fails. `active` means
the artifact remains in the current evidence ledger; it does not authorize a
new experiment. A later correction or STOP decision takes precedence over an
earlier positive-looking report.
