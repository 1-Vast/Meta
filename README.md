# MetaSieve-DTA

Mechanism-first few-shot drug-target affinity research with a frozen
probability-law operator.

## Current status

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_AND_PAIR_COMPATIBILITY_IDENTIFIED
EXACT_RESIDUE_LOCALISATION_IDENTIFIED_IN_DEVELOPMENT
TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED
S2R_SYNTHETIC_BINARY_ORDINAL_ESTIMATOR_PASS
S3R_REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

MetaSieve remains a trainable deep-learning bioinformatics research system:
frozen protein and ligand encoders supply biological measurements, and small
trainable interaction heads are admitted only after falsification controls. Its
main mathematical object remains the frozen constrained probability-law
operator `A(F,z)=K(B(z)F(z))`.

The latest real structural experiment did not admit a new biological statistic.
The gauge-free direct-W learner scored `AP_bidir = 0.03588` against chance
`0.02547`; the gain was positive but below the registered practical margin and
did not replicate beyond B5 or shortcut controls. This result is scoped to
frozen ESM2 residue states plus a mean-pooled 41-D ligand atom representation.

## Repository boundaries

- `theory/FINAL_FROZEN_THEORY/`: authoritative mathematics.
- `model/`: passed mathematical, encoder and geometry primitives; no validated
  assembled DTA pipeline.
- `scripts/`: passed data, sealing, structure, geometry and governance tools.
- `research/`: preregistered or executed research stages, including S2R/S3R.
- `report/`: machine Gates, current status and evidence summaries.
- `history.md`: chronological failure and decision ledger.

No active training stage is authorized. Affinity values, heldout-B, R6,
few-shot adaptation, biological `z`, CSMO/Band and the frozen operator remain
untouched.

## Read first

1. `report/CURRENT_RESEARCH_STATUS.md`
2. `report/s7_l2b_r0r/PHASE2B_S3R_EVIDENCE_CONSOLIDATION.md`
3. `report/s7_l2b_r0r/PHASE2B_S3R_GATE.json`
4. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
5. `task.md`
6. `experiment.md`
7. `history.md`

## Verification

```powershell
conda run -n drug python -m pytest -q
```

Large third-party releases, embedding banks and caches are not redistributed;
see `DATA_AVAILABILITY.md`. Current consolidated regression: **134 passed**.
