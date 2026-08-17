#!/usr/bin/env bash
set -u
cd /d/MetaSieve
CONDA=D:/anaconda/Scripts/conda.exe

run_std () {
  ARM=$1; OUT=$2; SEED=$3; SCRIPT=$4; ARG=$5
  echo "[$(date +%T)] $ARM seed=$SEED starting"
  $CONDA run -n drug python "$SCRIPT" --arm "$ARG" --seed "$SEED" --steps 1200 --output "$OUT"
  if [ ! -f "$OUT/RESULT.json" ]; then echo "[$(date +%T)] $ARM seed=$SEED FAILED"; exit 1; fi
  echo "[$(date +%T)] $ARM seed=$SEED done"
}

run_std T2 report/meta_fewshot/stageG2_20260817/T2/seed20260816 20260816 tools/research/stageD_level_panel/train_staged.py T2
run_std T2 report/meta_fewshot/stageG2_20260817/T2/seed20260817 20260817 tools/research/stageD_level_panel/train_staged.py T2
run_std G report/meta_fewshot/stageG2_20260817/G/seed20260816 20260816 tools/research/stageG_esm650/train_stageg.py G
run_std G report/meta_fewshot/stageG2_20260817/G/seed20260817 20260817 tools/research/stageG_esm650/train_stageg.py G
echo "[$(date +%T)] G2 training done"
