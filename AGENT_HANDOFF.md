# Agent handoff

Start from `PROJECT_SUMMARY.md`, then `task.md`, then `history.md`.

Use the `drug` conda environment for verification and experiments:

```powershell
conda run -n drug python main.py status
conda run -n drug python main.py verify tests
```

Binding state:

```text
MARGINAL_OR_SLOT_RECALIBRATION_ONLY
NO_VALIDATED_END_TO_END_FEWSHOT_DTA_MODEL
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
```

Do not run R1 affinity training, connect the exact distance residual to V1,
restore centered ridge/full MAML/test-time scheduling, or treat scheduler/noise
baselines as production components. Record terminal decisions in `history.md`
and keep duplicate payloads out of Git.
