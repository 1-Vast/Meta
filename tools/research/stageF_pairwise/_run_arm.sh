#!/usr/bin/env bash
set -u
cd /d/MetaSieve
CONDA=D:/anaconda/Scripts/conda.exe
run_arm () {
  ARM=$1
  OUT=$2
  STEPS=$3
  echo "[$(date +%T)] $ARM starting"
  $CONDA run -n drug python tools/research/stageF_pairwise/train_stagef.py --arm "$ARM" --steps "$STEPS" --output "$OUT"
  RC=$?
  if [ ! -f "$OUT/RESULT.json" ]; then
    echo "[$(date +%T)] $ARM FAILED (no RESULT.json, rc=$RC)"
    exit 1
  fi
  D:/anaconda/Scripts/conda.exe run -n drug python -c "import json; d=json.load(open('$OUT/RESULT.json')); r=d['report']; assert r['arm']=='$ARM'; assert r['gpu_probe']['torch_cuda_is_available']; assert r['gpu_probe']['batch_device_check']; print('OK best_internal_val_mse_pk=', r['best_internal_val_mse_pk'], 'trainable=', r['trainable_parameters'])" || { echo "[$(date +%T)] $ARM RESULT validation failed"; exit 1; }
  echo "[$(date +%T)] $ARM done"
}
run_arm "$1" "$2" "$3"
