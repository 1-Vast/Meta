#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/autodl-tmp/MetaSieve
PYTHON=/root/miniconda3/bin/python3.12
EXPECTED_THEORY=3d660448a585662083979c198d42258466cdcca7e0aab197095800cc2d42501e
THEORY=theory/FINAL_FROZEN_THEORY/00_CORE_THEORY/FINAL_THEORY_COMPLETE.md
RUNS=report/phase4_runs
LOGS=report/phase4_parallel

cd "$ROOT"
mkdir -p "$RUNS" "$LOGS"
exec 9>report/.phase4_remote.lock
if ! flock -n 9; then
  echo "another remote Phase 4 run holds report/.phase4_remote.lock" >&2
  exit 1
fi
echo $$ > "$LOGS/supervisor.pid"
trap 'rm -f "$LOGS/supervisor.pid"' EXIT

check_theory() {
  local actual
  actual=$(sha256sum "$THEORY" | cut -d' ' -f1)
  if [[ "$actual" != "$EXPECTED_THEORY" ]]; then
    echo "frozen theory hash mismatch: $actual" >&2
    exit 1
  fi
}

run_job() {
  local dataset=$1
  local seed=$2
  local log="$LOGS/${dataset}_s${seed}.log"
  echo "started dataset=$dataset seed=$seed time=$(date --iso-8601=seconds)" > "$log"
  PYTHONUNBUFFERED=1 CUDA_VISIBLE_DEVICES=0 "$PYTHON" \
    research/phase4_experiment.py run "$dataset" --seed "$seed" \
    --out-dir "$RUNS" >> "$log" 2>&1
  check_theory
  echo "completed dataset=$dataset seed=$seed time=$(date --iso-8601=seconds)" >> "$log"
}

run_wave() {
  local dataset=$1
  local pids=()
  local failed=0
  for seed in 11 23 37; do
    run_job "$dataset" "$seed" &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
      failed=1
    fi
  done
  if [[ "$failed" -ne 0 ]]; then
    echo "Phase 4 wave failed: $dataset" >&2
    exit 1
  fi
}

check_theory
"$PYTHON" -m pytest -q tests/test_phase4_interface.py tests/test_phase3_protocol.py
run_wave DAVIS
check_theory
run_wave KIBA
check_theory
"$PYTHON" research/phase4_experiment.py summarize \
  --datasets DAVIS KIBA --seeds 11 23 37 --out-dir "$RUNS"
"$PYTHON" -m pytest -q
check_theory
date --iso-8601=seconds > "$LOGS/completed.txt"
