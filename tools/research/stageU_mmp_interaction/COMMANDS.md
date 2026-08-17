# Stage U — exact command record

Environment: `D:\anaconda\envs\drug\python.exe` (conda env `drug`),
Python 3.11.15, torch 2.6.0+cu124 (CUDA available but no U2 training was
authorized), RDKit 2023.09.6. Working directory `D:\MetaSieve`.

```bash
# 1. Freeze preregistration (done before any U0 statistic was read)
#    SHA-256: fdc0a830aa92882d07b9aea50f22a4c72fc6d93f92c55a3be6bc15cd6a645c11

# 2. U0 measurement-reliability audit (builds U0_PROVENANCE_meta_train.jsonl.gz)
D:/anaconda/envs/drug/python -m tools.research.stageU_mmp_interaction.u0_reliability

# 3. Build the deterministic MMP observation cache
D:/anaconda/envs/drug/python -c "..."   # see observation_cache.py; equivalently
D:/anaconda/envs/drug/python -m tools.research.stageU_mmp_interaction.u0_census

# 4. U0 census, bipartite graph, coverage and frozen admission gate
D:/anaconda/envs/drug/python -m tools.research.stageU_mmp_interaction.u0_census

# 5. Structural verification
RUN_SLOW=1 D:/anaconda/envs/drug/python -m pytest tools/research/stageU_mmp_interaction/tests -q
```

## Outcome

- `U0_CENSUS.json` → `admission.all_pass = false`.
- Frozen degree-concentration gate failed: top-1 target and top-1 component
  each hold **29.63%** of same-panel fit observations against the frozen
  **25%** cap.
- Per the frozen stop rule, U1 (`u1_variance.py`) was not run, U2 was not
  implemented or trained, and no neural model was created.
- The stage is recorded in `U0_DECISION.json` and `REPORT.md`.

No command in this stage reads, opens or evaluates the sealed confirmation
split. No command modifies `model/` or production `scripts/`.
