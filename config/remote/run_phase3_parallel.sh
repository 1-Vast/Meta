#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/autodl-tmp/MetaSieve
PYTHON=/root/miniconda3/bin/python3.12
EXPECTED_THEORY=3d660448a585662083979c198d42258466cdcca7e0aab197095800cc2d42501e
THEORY=theory/FINAL_FROZEN_THEORY/00_CORE_THEORY/FINAL_THEORY_COMPLETE.md
RUNS=report/phase3_runs
LOGS=report/phase3_parallel

cd "$ROOT"
mkdir -p "$RUNS" "$LOGS"
echo $$ > "$LOGS/supervisor.pid"
trap 'rm -f "$LOGS/supervisor.pid"' EXIT

check_theory() {
  [[ "$(sha256sum "$THEORY" | cut -d' ' -f1)" == "$EXPECTED_THEORY" ]]
}

run_job() {
  local dataset=$1
  local seed=$2
  local log="$LOGS/${dataset}_s${seed}.log"
  echo "started dataset=$dataset seed=$seed time=$(date --iso-8601=seconds)" > "$log"
  PYTHONUNBUFFERED=1 CUDA_VISIBLE_DEVICES=0 "$PYTHON" scripts/phase3.py run \
    "$dataset" --seed "$seed" --profile base --out-dir "$RUNS" >> "$log" 2>&1
  check_theory
  echo "completed dataset=$dataset seed=$seed time=$(date --iso-8601=seconds)" >> "$log"
}

check_theory
pids=()
names=()
for spec in DAVIS:23 DAVIS:37 KIBA:11 KIBA:23 KIBA:37; do
  dataset=${spec%%:*}
  seed=${spec##*:}
  run_job "$dataset" "$seed" &
  pids+=("$!")
  names+=("$spec")
done

failed=0
for index in "${!pids[@]}"; do
  if ! wait "${pids[$index]}"; then
    echo "failed: ${names[$index]}" >&2
    failed=1
  fi
done
check_theory
if [[ "$failed" -ne 0 ]]; then
  exit 1
fi
date --iso-8601=seconds > "$LOGS/completed.txt"
