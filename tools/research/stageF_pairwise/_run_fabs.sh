#!/usr/bin/env bash
set -u
cd /d/MetaSieve
for i in $(seq 1 60); do
  if [ -f report/meta_fewshot/stageF_pairwise_20260817/F/RESULT.json ]; then break; fi
  sleep 15
done
if [ ! -f report/meta_fewshot/stageF_pairwise_20260817/F/RESULT.json ]; then
  echo "[$(date +%T)] F never completed; aborting F-ABS chain"
  exit 1
fi
bash tools/research/stageF_pairwise/_run_arm.sh F-ABS report/meta_fewshot/stageF_pairwise_20260817/F-ABS 1200
