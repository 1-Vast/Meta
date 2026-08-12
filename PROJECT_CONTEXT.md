# MetaSieve project context

This is the compact root-level context file. `README.md`, `task.md` and
`history.md` remain the detailed entry points.

## Status

```text
Task: unseen-target few-shot drug-target affinity prediction
Stage: R0C_EXACT_DISTANCE_CONFIRMATION_COMPLETE
Status: MARGINAL_OR_SLOT_RECALIBRATION_ONLY
Training authorized: false
NO_VALIDATED_END_TO_END_FEWSHOT_DTA_MODEL
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
```

Use the `drug` conda environment:

```powershell
conda run -n drug python main.py status
conda run -n drug python main.py verify tests
```

## Experimental Contract

- Predict query-ligand affinity for an unseen target after `k=1/2/3/5`
  support affinities.
- Keep query labels out of model selection.
- Keep measurement modalities separate.
- Train on target episodes or complete panels, not IID activity rows.
- Use support-only positive-ridge adaptation with task freedom `d <= 5`.
- Require ligand-only, wrong-partner, foreign-support, permuted-support,
  endpoint-separated and source-separated controls.

No raw pair map or arbitrary neural latent enters `z`; only an independently
confirmed partner-specific biological statistic may be proposed for
`A(F,z)=K(B(z)F(z))`.

## Data Boundary

Git tracks source code, tests, frozen theory, compact reports and manifests.
It does not redistribute third-party releases, feature banks, tensors,
checkpoints, prediction dumps, row-level experiment payloads, local caches,
credentials or machine-specific remote configs.

Tracked data is limited to provenance, license, checksum, warning and manifest
files plus compact `RESULT.json` evidence.

## Archive Boundary

Retired scripts, detailed result payloads, checkpoints and draw-metric traces
are removed from the active Git boundary. Historical report summaries and
archive manifests remain for audit context.

## Do Not Reopen

```text
R1 affinity training from the rejected exact distance residual
Exact residual connection to V1
Centered ridge as a replacement for the uncentered positive ridge
Full MAML/test-time scheduling
Scheduler or fixed-noise baselines as production components
```
