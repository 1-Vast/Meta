# MetaSieve-DTA

Mechanism-first few-shot drug–target affinity research with a frozen
probability-law operator.

## Current status

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_IDENTIFIED
PAIR_COMPATIBILITY_IDENTIFIED
FROZEN_ESM2_B5_DEVELOPMENT_GATE_PASS_6_OF_6
EXACT_RESIDUE_LOCALISATION_IDENTIFIED_IN_DEVELOPMENT
TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED_IN_LABELS
B5_RESIDUE_MARGINAL_IS_GENERIC_POCKET
B5_LIGAND_DEPENDENCE_CONFINED_TO_THE_COUPLING_TERM
TEACHER_EDGE_COUPLING_NOT_IDENTIFIED
EXACT_RESIDUE_ATOM_COUPLING_NOT_IDENTIFIED
LABEL_SEMANTICS_NOT_AMBIGUOUS
AFFINITY_ENERGETICS_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

Frozen ESM2 resolved the B4 residue-representation deficit and B5 passed all six
registered structural development Gates. The Phase 2A audit then answered the
attribution question: the MONN labels **are** ligand-conditioned at the residue
level (ΔJ `+0.258` [LCB `+0.234`] over a replicate noise floor measured from the
data), B5's residue marginal is **not** (a wrong ligand retains 89% of it), and
edge-level coupling is not identifiable in the teacher or in B5. Affinity value
beyond ligand-only and wrong-protein controls remains untested.

## Repository boundaries

- `theory/FINAL_FROZEN_THEORY/`: authoritative frozen mathematics.
- `model/`: passed mathematical, encoder and P1B geometry primitives; no
  validated assembled DTA pipeline.
- `scripts/`: passed data, sealing, structure, geometry and source-governance
  workflows.
- `research/`: completed S7/L2B development code, the Phase 2A audit
  (`pa0`–`pa5`), and unvalidated follow-on research. The Phase 2B
  preregistration is frozen but not authorized and has no implementation.
- `report/`: current status, split protocol and compact PASS evidence.
- `history.md`: authoritative experimental and failure ledger.

Terminal-negative implementations and duplicate reports were removed after
consolidation.  They remain recoverable from Git history at `3281780`,
`12a2765`, and `608decf`.

## Active stage

S7/L2B Phase 0, Phase 1 and the audit-only Phase 2A are complete. Phase 2A
trained nothing, read no affinity value, and returned
`LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING`, whose single mandated
action — writing and freezing a preregistration for one ligand-conditioned
residue residual head — has been taken
(`research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md`). Phase 2B is
**registered, not authorized**: no Phase 2B code exists and no run has occurred.

No real affinity value, DAVIS/recipient label, confirmation scoring, production
`z`, CSMO/Band change or P2–P4 stage is authorized.

## Read first

1. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
2. `report/CURRENT_RESEARCH_STATUS.md`
3. `report/s7_l2b_r0r/PHASE2A_SYNTHESIS.md`
4. `report/s7_l2b_r0r/PHASE1_EVIDENCE_CONSOLIDATION_AND_PHASE2A_TRIAGE.md`
5. `task.md`
6. `experiment.md`
7. `EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md`
8. `history.md`

## Verification

```powershell
conda run -n drug python -m pytest -q
```

The consolidated suite currently passes 75 tests.  Large third-party releases,
embedding banks and caches are not redistributed; see `DATA_AVAILABILITY.md`.
