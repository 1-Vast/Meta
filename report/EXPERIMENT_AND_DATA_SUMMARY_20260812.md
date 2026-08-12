# MetaSieve experiment and data summary - 2026-08-12

## Repository structure

```text
model/       Verified operator, encoder and geometry primitives.
contracts/   Data and mechanism schemas.
scripts/     Data acquisition, sealing, audit and preprocessing utilities.
research/    Active research runners and reproducibility scripts.
tests/       Pytest coverage for model, data, research and CLI contracts.
dataset/     Manifest-only provenance in Git; bulk data remains local.
report/      Active experiment JSON reports and concise result evidence.
archive/     Historical reports, failed lines and retired research evidence.
theory/      Frozen theory and development archive.
```

## Active experiment evidence

### Source affinity and SAR-delta

```text
report/source_affinity/chembl_assay_sardelta_gate1_target_smoke/RESULT.json
TERMINAL_VERDICT = CHEMBL_ASSAY_SARDELTA_GATE1_PASS
development pairs = 415
correct MSE = 0.4267805814862487
zero-delta MSE = 0.7864955120880043
component macro reduction = 0.2972134666811634
LCB95 = 0.14637099289452893
V1 integration authorized = false
```

```text
report/crossed_interaction/bindingdb_sardelta_cq_bridge_gate1_smoke/RESULT.json
TERMINAL_VERDICT = BINDINGDB_SARDELTA_CQ_BRIDGE_GATE1_PASS
train pairs = 20423
development pairs = 1033
correct MSE = 0.2337137884759225
zero-delta MSE = 0.5800715132288407
component macro reduction = 0.3767396232892215
LCB95 = 0.27017239997568654
V1 integration authorized = false
```

### Failed crossed-interaction controls retained

The following active evidence remains in `report/crossed_interaction/` because
it records negative controls and failed Gates, not disposable code:

```text
chembl_affinity_teacher_cq_observable_gate1*
sardelta_potential_cq_observable_gate1_smoke
sardelta_edge_cq_observable_gate1_smoke
attention/localizer, kmer, physchem, family/function context, PLM slot,
pocket prototype and structure pocket prior CQ observable sweeps
```

The 2026-08-12 organization audit is
`report/RESEARCH_ORGANIZATION_20260812.md`.

### Meta-fewshot V1

```text
report/meta_fewshot/main_v1_vectorized_reproduction/MAIN_V1_COLD_TARGET_RESULT.json
TERMINAL_VERDICT = COLD_TARGET_FEWSHOT_V1_NOT_YET_GOOD
```

```text
report/meta_fewshot/main_v1_support_only_mlp64_d4_final/MAIN_V1_COLD_TARGET_RESULT.json
TERMINAL_VERDICT = COLD_TARGET_FEWSHOT_V1_NOT_YET_GOOD
RMSE(k=1/2/3/5) = 1.495/1.357/1.290/1.230
CI(k=1/2/3/5) = 0.549/0.553/0.558/0.567
```

```text
report/meta_fewshot/main_v1_ligand_residual_mlp64_d4/MAIN_V1_COLD_TARGET_RESULT.json
TERMINAL_VERDICT = COLD_TARGET_FEWSHOT_V1_NOT_YET_GOOD
pair_feature_mode = ligand_residual
```

Large checkpoints and label-free prediction payloads are local reproducibility
artifacts and are excluded from Git by `.gitignore`; their checksums remain in
the retained JSON reports.

## Dataset evidence

Git tracks dataset provenance, warnings, licenses and manifests only. Raw
third-party releases, embedding banks, tensors and caches remain local.

Representative tracked manifests:

```text
dataset/raw/source_affinity/chembl37_sqlite_v1/release_manifest.json
dataset/processed/source_affinity/energy_pilot_v1/corpus_manifest.json
dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus/manifest.json
dataset/processed/meta_fewshot/bindingdb_ki_v1_development/manifest.json
dataset/processed/correspondence_router/r0c_exact_geometry_v1/manifest.json
dataset/sealed/DAVIS_v1/manifest.json
```

## Cleanup decision

Removed from the Git boundary:

```text
Python bytecode caches
pytest cache/cycle outputs
active report checkpoint directories
compressed label-free prediction dumps
```

Retained:

```text
failed experiment scripts and tests
Gate RESULT.json files
dataset manifest/provenance files
archive manifests and historical negative evidence
```

No admitted end-to-end few-shot DTA model exists yet. The current positive
result is scoped SAR-delta transfer evidence; the current V1 family remains
trainable but scientifically failed closed.
