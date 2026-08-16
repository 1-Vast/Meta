# MetaSieve project context

This is the compact root-level context file. The current authority is
`README.md`, `PROJECT_SUMMARY.md`, `task.md`, and
`docs/PROJECT_FILE_ORGANIZATION.md`; `history.md` is provenance.

## Status

```text
Task: unseen-target few-shot drug-target affinity prediction
Stage: QPSMP_BPSF_DEVELOPMENT
Status: TRAINABLE_INTERFACE_NOT_ADMITTED
Development training authorized: true
NO_ADMITTED_COLD_TARGET_UTILITY_CLAIM
```

Use the `drug` conda environment:

```powershell
conda run -n drug python main.py status
conda run -n drug python main.py verify tests
```

## Experimental Contract

- Predict query-ligand affinity for an unseen target after `k=0/1/2/3/5`
  support observations.
- Keep query labels out of model selection.
- Keep measurement modalities separate.
- Train on target episodes or complete panels, not IID activity rows.
- Use an interaction-first BPSF zero-shot endpoint and a positive label-locked
  residual kernel; no ridge, closed-form solve, inner loop, or test-time
  gradient is permitted.
- Require additive, level, ligand-only, SAR-cut, wrong-protein,
  shuffled-protein, foreign-support, and design-nuisance controls.

## Data Boundary

Git tracks source code, tests, compact reports and manifests.
It does not redistribute third-party releases, feature banks, tensors,
checkpoints, prediction dumps, row-level experiment payloads, local caches,
credentials or machine-specific remote configs.

Tracked data is limited to provenance, license, checksum, warning and manifest
files plus compact `RESULT.json` evidence.

## Archive Boundary

Historical theory, retired scripts, result payloads, checkpoints and traces
are archive material, not active runtime authority. Historical summaries and
archive manifests remain for audit context.

## Do Not Reopen

```text
Closed-form or ridge adaptation in the active path
Full MAML or test-time scheduling
Persistent target-specific parameter memory
Atomic 3D claims without legal common-frame coordinates
```
