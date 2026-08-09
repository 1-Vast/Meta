# MetaSieve-DTA

Mechanism-first few-shot drug–target affinity research with a frozen
probability-law operator.

## Current status

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_IDENTIFIED
PAIR_COMPATIBILITY_IDENTIFIED
PAIR_LOCAL_P1B_MECHANISM_OBSERVABILITY_NOT_TESTED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
EXTERNAL_S5_S9_CLAIMS_NOT_REPRODUCED
```

The current scientific boundary is not whether the encoders see the protein.
P1B already established correct-protein contact/distance geometry.  The open
question is whether its atom-local, residue-local and pair-local observables can
recover a protein-specific structural mechanism that subsequently provides
affinity value beyond ligand-only and wrong-protein controls.

## Repository boundaries

- `theory/FINAL_FROZEN_THEORY/`: authoritative frozen mathematics.
- `model/`: passed mathematical, encoder and P1B geometry primitives; no
  validated assembled DTA pipeline.
- `scripts/`: passed data, sealing, structure, geometry and source-governance
  workflows.
- `research/`: the active S5 preregistration and one conditional metadata-only
  SSL data-feasibility preregistration. Unvalidated code remains here.
- `report/`: current status, split protocol and compact PASS evidence.
- `history.md`: authoritative experimental and failure ledger.

Terminal-negative implementations and duplicate reports were removed after
consolidation.  They remain recoverable from Git history at `3281780`,
`12a2765`, and `608decf`.

## Active stage

`P1R2B-S5_LOCAL_MECHANISM_OBSERVABILITY` tests the actual frozen P1B local
contract.  It begins with chain/mapping and pseudo-teacher audits, then an
observability ladder and synthetic trainability control.  Lightweight GPU
distillation is conditional on those checks.

No real affinity value, DAVIS/recipient label, production `z`, CSMO/Band change
or P2–P4 stage is authorized by S5.

## Read first

1. `report/CURRENT_RESEARCH_STATUS.md`
2. `task.md`
3. `experiment.md`
4. `EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md`
5. `history.md`
6. `report/SSL_DETAILED_REPORT_EVIDENCE_AUDIT.md`
7. `report/SSL_TEST_PLAN_REVIEW_AND_EXECUTION_BOUNDARY.md`

## Verification

```powershell
conda run -n drug python -m pytest -q
```

The consolidated suite currently passes 75 tests.  Large third-party releases,
embedding banks and caches are not redistributed; see `DATA_AVAILABILITY.md`.
