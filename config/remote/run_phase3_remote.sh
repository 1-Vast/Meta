#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/autodl-tmp/MetaSieve
PYTHON=/root/miniconda3/bin/python3.12
EXPECTED_THEORY=3d660448a585662083979c198d42258466cdcca7e0aab197095800cc2d42501e
THEORY=theory/FINAL_FROZEN_THEORY/00_CORE_THEORY/FINAL_THEORY_COMPLETE.md
RUNS=report/phase3_runs
ARCHIVE=report/phase3_archive/remote_rerun_20260804
LOG=report/phase3_remote.log

cd "$ROOT"
exec 9>report/.phase3_remote.lock
if ! flock -n 9; then
  echo "another remote Phase 3 run holds report/.phase3_remote.lock" >&2
  exit 1
fi
echo $$ > report/phase3_remote.pid
trap 'rm -f report/phase3_remote.pid' EXIT

check_theory() {
  local actual
  actual=$(sha256sum "$THEORY" | cut -d' ' -f1)
  if [[ "$actual" != "$EXPECTED_THEORY" ]]; then
    echo "frozen theory hash mismatch: $actual" >&2
    exit 1
  fi
}

mkdir -p "$RUNS" "$ARCHIVE"
for dataset in DAVIS KIBA; do
  for seed in 11 23 37; do
    previous="$RUNS/phase3_${dataset}_s${seed}.json"
    if [[ -f "$previous" ]]; then
      cp -p "$previous" "$ARCHIVE/"
    fi
  done
done

exec > >(tee -a "$LOG") 2>&1
echo "remote Phase 3 started: $(date --iso-8601=seconds)"
"$PYTHON" --version
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
check_theory

for dataset in DAVIS KIBA; do
  for seed in 11 23 37; do
    echo "run started: dataset=$dataset seed=$seed time=$(date --iso-8601=seconds)"
    PYTHONUNBUFFERED=1 CUDA_VISIBLE_DEVICES=0 "$PYTHON" scripts/phase3.py run \
      "$dataset" --seed "$seed" --profile base --out-dir "$RUNS"
    echo "run completed: dataset=$dataset seed=$seed time=$(date --iso-8601=seconds)"
    check_theory
  done
done

"$PYTHON" scripts/phase3.py summarize \
  --datasets DAVIS KIBA --seeds 11 23 37 --out-dir "$RUNS"
"$PYTHON" -m pytest -q
check_theory
echo "remote Phase 3 completed: $(date --iso-8601=seconds)"
