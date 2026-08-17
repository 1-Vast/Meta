# Exact commands for this stage

Environment: the `drug` conda environment. Working directory: `D:\MetaSieve`.
Every command below was run in the order shown.

## Phase 0 — governance repair

```bash
conda run -n drug python -m pytest tools/tests/test_meta_test_seal_contract.py -q
conda run -n drug python -m pytest tools/tests -q
RUN_SLOW=1 RUN_RESEARCH_GATES=1 conda run -n drug python -m pytest tools/tests -q
conda run -n drug python -m scripts.audit_research_record --skip-loading
```

Seal-repair reproduction — the identical R14 audit re-run under the repaired
fail-closed seal, whose 105 numeric A0 fields must be bit-identical to the
recorded ones:

```bash
conda run -n drug python -m scripts.r14_dispersion_audit \
  --arm "A0=report/meta_fewshot/stageR3R4_level_shape_20260815/A0_incumbent_seed20260815/checkpoint.pt" \
  --arm "A0=report/meta_fewshot/stageR3R4_level_shape_20260815/A0_incumbent_seed20260816/checkpoint.pt" \
  --arm "A0=report/meta_fewshot/stageR3R4_level_shape_20260815/A0_incumbent_seed20260817/checkpoint.pt" \
  --split-directory dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1 \
  --output tools/research/a2_readiness_v2/SEAL_REPAIR_REPRODUCTION.json
```

Correct the seven artifacts carrying the false seal claim (numbers untouched;
the script asserts non-seal content is unchanged before writing):

```bash
conda run -n drug python -m tools.research.a2_readiness_v2.correct_seal_claims
```

## Phase 1 — reproduce the diagnosis without leakage

```bash
conda run -n drug python -m tools.research.a2_readiness_v2.branch_ordering_v2 \
  --output tools/research/a2_readiness_v2/BRANCH_ORDERING_V2_meta_val.json
```

Three A0 seeds and ten independent random initialisations, five donor strata,
identical-protein and repeated-forward floors, shuffled-label /
scrambled-protein / foreign-ligand controls, novelty and scaffold strata.

## Phase 2 — internal representations

```bash
conda run -n drug python -m tools.research.a2_readiness_v2.extract_features \
  --split meta_train --output tools/research/a2_readiness_v2/features
conda run -n drug python -m tools.research.a2_readiness_v2.extract_features \
  --split meta_val   --output tools/research/a2_readiness_v2/features
conda run -n drug python -m tools.research.a2_readiness_v2.representation_probe \
  --features tools/research/a2_readiness_v2/features \
  --output tools/research/a2_readiness_v2/REPRESENTATION_PROBE_meta_val.json
```

## Phase 3 — causal attention/readout localization

```bash
conda run -n drug python -m tools.research.a2_readiness_v2.attention_causal_audit \
  --output tools/research/a2_readiness_v2/ATTENTION_CAUSAL_meta_val.json
```

## Stage tests

```bash
conda run -n drug python -m pytest tools/research/a2_readiness_v2/tests -q
```

## Not run, and why

`meta_test` labels were used for no fitting, selection or reported metric; the process-isolation incident remains open. No training was performed: every measurement in
this stage is a forward pass on frozen checkpoints, except the small diagnostic
probes of Phase 2, which are trained on `meta_train` component folds by
ordinary gradient descent and scored once on `meta_val`.

Stage P of `PREREGISTRATION_V2.md` is **not** run and is not authorized by this
package.
