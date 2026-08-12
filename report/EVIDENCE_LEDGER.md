# Evidence Ledger

| Evidence | Artifact | Scope | Status |
|---|---|---|---|
| X1 | `crossed_interaction/bindingdb_rectangle_interaction_x1_theory_corrected_20260812/RESULT.json` | Governed development panels; observed labels | Diagnostic, not latent-mechanism evidence |
| F163 | `crossed_interaction/bindingdb_rectangle_descriptor_g2_plmslots_ridge1e6_20260812/RESULT.json` | 12 consumed development components | G2 fail-closed; only zero contrast passed |
| QPSMP analytic smoke | `meta_fewshot/qpsmp_core_smoke_separated_20260812/RESULT.json` | Main-v0 meta-test smoke | Fail-closed; comparator only |
| QPSMP neural interface | `../model/qpsmp_meta.py`, `../scripts/train_qpsmp.py`, `../scripts/evaluate_qpsmp.py` | Model and governed CLI contract | Trainability/invariant tests pass; no admission |
| QPSMP evidence-gated smoke | `meta_fewshot/qpsmp_meta_channels_sarcut_20step_20260812/RESULT.json` | Fixed-bank CPU implementation smoke; six test components | SAR-only gain is slightly negative; correct-state specificity unrecognized; no safety or inferential Gate authorized |
| QPSMP repaired k5 seed 20260831 | `meta_fewshot/qpsmp_meta_final_formula_k5_seed20260831_20260812/RESULT.json` | Consumed main-v0 development; 1,000 updates, all eligible targets, five draws | Full beats level; SAR and all state controls positive in this seed only |
| QPSMP repaired k5 seed 20260832 | `meta_fewshot/qpsmp_meta_final_formula_k5_seed20260832_20260812/RESULT.json` | Consumed main-v0 development; matched budget | Full beats level; SAR-only gain negative; no admission |
| QPSMP repaired k5 seed 20260833 | `meta_fewshot/qpsmp_meta_final_formula_k5_seed20260833_20260812/RESULT.json` | Consumed main-v0 development; matched budget | Full beats level; SAR-only and permutation controls fail; no admission |
| QPSMP nested k={1,2,3,5} development | `meta_fewshot/qpsmp_nested_k1235_development_20260812/RESULT.json` | Consumed main-v0; shared checkpoints, nested support, fixed query, 42 targets/6 components | Full-level point estimates positive at every k, but all component LCBs fail; SAR/specificity unrecognized; promotion fail-closed |

Failed or superseded experiments remain retained for auditability. They are not impossibility proofs
for richer architectures or different data designs.
