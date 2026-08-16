# MetaSieve Project Summary

This file is the current status entry point. Numerical claims are controlled by the corresponding
`RESULT.json` and manifest files. Historical theory now lives under
`archive/theory/` and is design provenance, not a frozen neural contract.

## Objective

Learn a shared deep meta-learning model for few-shot affinity prediction on proteins absent from
training. Deployment may use only the protein sequence or legal structure representation, ligand
molecular information, context, and a disjoint support set. Target identifiers, query labels, and
target-specific parameter memory are prohibited.

## Current population and incumbent

The governed **double-cold protocol**
`dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1` is the
development/confirmation population (meta_train 5,643 cells / 346 targets;
meta_val 41 targets / 19 components; meta_test 22 targets / 10 components,
physically sealed and pristine — never opened).

The retained incumbent is the **Stage R3/R4 `similarity_only` grammar**
(three checkpoints, `report/meta_fewshot/stageR3R4_level_shape_20260815/
A0_incumbent_seed*/`, 1200 steps). Double-cold meta_val: k=0 MSE 2.149 /
CI 0.580 / Spearman 0.223 / calibration 1.236 / shape 0.913. It is a
development incumbent, not an admitted model.

## Current status (2026-08-16)

- The R5-R8 relative-transport/gate model family was **closed** under its
  preregistered gates (R6a/R6b/R7 transport mechanisms all measured
  deployment-inert; R8 strong-shape arm ties A0 at k=0 but regresses CI).
  Full cycle summary: `report/meta_fewshot/README.md`.
- The shape-first training method is retained as the project's first
  measured within-target shape source (shape 0.943 -> 0.896; k=5
  activity-cliff sign 0.768, best on record).
- R9 localized the CI regression to the mid-similarity band (0.4-0.6) and
  showed the x4 activity-cliff pair weight is a net negative for ranking;
  the cliff-weight dose response found no single dose passing both
  preregistered advance gates (C1: CI 0.562; C2: k=0 2.119, calib 1.218).
- R10 (in progress) tests the next single variable: `shape_variance_weight`
  1.5 -> 0.5 on the C1 base, three seeds, 1200 steps, via the smoke-first
  stage runner.
- The complete maintained suite passes 393 tests (79 pytest modules).

No Cold Target admission, SOTA, or confirmation claim is authorized.
`meta_test` opens once, only after every preregistered meta_val gate passes.

## Core Repository Surface

1. `docs/PROJECT_FILE_ORGANIZATION.md`: ownership and active-file map.
2. Active models: `model/interaction_grammar.py` + `model/similarity_grammar.py`
   (incumbent family), `model/level_shape.py`, `model/reltransport.py`
   (closed, retained for evidence), `model/qpsmp_meta.py` (control arm).
3. `scripts/qpsmp_data.py`, `scripts/train_qpsmp.py`,
   `scripts/train_reltransport.py`, `scripts/train_grammar_shape.py`,
   `scripts/stageR6_compare_arms.py`, `scripts/stageR9_pair_audit.py`,
   `scripts/stage_smoke.py`, `scripts/run_stage.py`: current data, training,
   evaluation, auditing and orchestration entry points.
4. `report/CURRENT_MODEL_EVIDENCE.md`: consolidated model and experiment ledger.
5. `archive/theory/`: historical mathematical provenance.

Directory-specific indexes now live in `model/`, `scripts/`, `research/`,
`report/`, `dataset/`, `tests/`, `contracts/`, `config/`, `docs/`, and `LLM/`.
Experiment outputs belong only under `report/`; research directories contain
code, not result trees.
