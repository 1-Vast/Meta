#!/usr/bin/env bash
set -u
cd /d/MetaSieve
CONDA=D:/anaconda/Scripts/conda.exe
EVAL=tools/research/stageD_level_panel/evaluate_staged.py

eval_one () {
  NAME=$1; CKPT=$2; OUT=$3; BANK=$4
  echo "[$(date +%T)] eval $NAME"
  if [ -z "$BANK" ]; then
    $CONDA run -n drug python "$EVAL" "$CKPT" --output "$OUT"
  else
    $CONDA run -n drug python "$EVAL" "$CKPT" --protein-bank "$BANK" --output "$OUT"
  fi
}

eval_one T2_s15 report/meta_fewshot/stageD_level_panel_20260817/T2/checkpoint.pt tools/research/stageG_esm650/T2_s15.rows.jsonl ""
eval_one T2_s16 report/meta_fewshot/stageG2_20260817/T2/seed20260816/checkpoint.pt tools/research/stageG_esm650/T2_s16.rows.jsonl ""
eval_one T2_s17 report/meta_fewshot/stageG2_20260817/T2/seed20260817/checkpoint.pt tools/research/stageG_esm650/T2_s17.rows.jsonl ""
eval_one G_s15 report/meta_fewshot/stageG_esm650_20260817/G/checkpoint.pt tools/research/stageG_esm650/G_s15.rows.jsonl tools/runtime/esm2_t33_650M_protein_bank
eval_one G_s16 report/meta_fewshot/stageG2_20260817/G/seed20260816/checkpoint.pt tools/research/stageG_esm650/G_s16.rows.jsonl tools/runtime/esm2_t33_650M_protein_bank
eval_one G_s17 report/meta_fewshot/stageG2_20260817/G/seed20260817/checkpoint.pt tools/research/stageG_esm650/G_s17.rows.jsonl tools/runtime/esm2_t33_650M_protein_bank
echo "[$(date +%T)] G2 evaluation done"
