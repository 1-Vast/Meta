# Stage V — exact commands, in order

Environment: conda env `drug` (`D:\anaconda\envs\drug\python.exe`), Python
3.11.15, torch 2.6.0+cu124, RDKit 2023.09.6, numpy 1.26.4, scipy 1.17.1,
git commit `5bb3736`. Full digests in `ENVIRONMENT.json`.

All commands run from the repository root.

```bash
# Stage 0 - forensic reconciliation of Stage T's core-blind key.
# Reads Stage T's construction code unmodified; writes STAGE0_FORENSICS.json.
python -m tools.research.stageV_core_mmp.stage0_forensics

# V0 census (core-inclusive key) + V0b evaluability + V1 interaction variance.
# Writes V0_V1_RESULT.json. All three gates evaluated in one pass.
python -m tools.research.stageV_core_mmp.v0_census

# Structural and leakage tests (31, including the slow subprocess and
# corpus-mounting cases).
RUN_SLOW=1 python -m pytest tools/research/stageV_core_mmp/tests -q

# Environment and artifact hashes.
python - <<'PY'
import hashlib, json, platform, subprocess, sys
from pathlib import Path
import numpy, torch, rdkit, scipy
# ... see ENVIRONMENT.json for the recorded output
PY

# Maintained suite (unchanged by this stage).
python main.py verify tests
```

## Not run, and why

```bash
# V2 neural arms - NOT RUN.
# The frozen V1 stop rule fired: theta = -0.4059 [-0.6889, -0.0577], resolved
# below zero, so the interaction variance does not exceed supervision noise.
# V0 also failed its concentration caps and the primary evaluation surface has
# 32 rows over 4 protein components (< the inherited 100-row evaluability rule).
# No model was built or trained.
```

## Prior-stage commands reproduced for the forensic comparison

```bash
# Stage T artifacts are read-only inputs to Stage 0; they were not re-run.
# tools/research/stageT_mmp/T1_CENSUS.json, T2_RESULT.json are the originals.
```
