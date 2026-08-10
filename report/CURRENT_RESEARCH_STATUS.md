# MetaSieve current research status

Updated: 2026-08-10.

## Current verdict

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_AND_PAIR_COMPATIBILITY_IDENTIFIED
FROZEN_ESM2_EXACT_RESIDUE_LOCALISATION_PASS_IN_DEVELOPMENT
TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED
SYNTHETIC_BINARY_ORDINAL_ESTIMATOR_TRAINABLE
REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

## Earliest failed boundary

Phase 2A proved that same-construct scaffold-distinct ligands change the MONN
residue masks beyond the replicate noise floor. S2R then repaired the synthetic
optimization defect with a bounded gauge-free direct matrix and passed a sealed
synthetic seed (`AP_bidir = 0.6620`).

S3R transferred that exact estimator to real structural labels. On 46,818 pairs
from 112 closure components, the candidate scored `0.03588` versus chance
`0.02547`. The gain was positive but only `+0.01041 [LCB +0.00692]`, below the
registered `+0.05` margin. The candidate did not beat frozen B5, foreign-ligand,
context-corrupted or trained permuted-label controls by their margins.

```text
REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
```

The failure is not numerical: module participation, train-only label access,
unit norm, stream equality and bit-exact prediction replay passed. It is scoped
to the current frozen ESM2 residue states plus mean-pooled 41-D ligand atom
features and direct-W ordinal estimand. It does not close all sequence+2D models.

## Current boundary

No active training stage is authorized. A future proposal may change one axis:
replace ligand mean pooling with a frozen graph-aware 2D ligand statistic while
holding the protein branch, estimator, loss, split and Gates fixed. It requires
a new preregistration and cannot reuse heldout-A as confirmation evidence.

Heldout-B, R6, affinity values, few-shot sectioning, biological `z`, CSMO/Band
and the frozen law operator remain unopened and unchanged.

## Canonical evidence

1. `report/s7_l2b_r0r/PHASE2B_S3R_GATE.json`
2. `report/s7_l2b_r0r/PHASE2B_S3R_EVIDENCE_CONSOLIDATION.md`
3. `report/s7_l2b_r0r/PHASE2B_S2R_VERDICT.json`
4. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
5. `history.md`
