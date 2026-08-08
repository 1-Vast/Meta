# MetaSieve-DTA

Mechanism-first few-shot drug-target affinity research with a frozen convex
law-valued operator.

Updated: 2026-08-07

## Current Verdict

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_IDENTIFIED
PAIR_COMPATIBILITY_IDENTIFIED
FIXED_RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED
TASK_LOCAL_RADIAL_AFFINITY_HEADROOM_OBSERVED
CORRECT_PARTNER_AFFINITY_SECTION_NOT_IDENTIFIED
FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED
CHEMBL_CROSSED_SOURCE_INTERACTION_UNDERDETERMINED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

The implementation is mathematically interface-compatible, but it is not yet
deeply integrated with bioinformatics. `model/` contains the frozen operator
primitives and the P1B-passing local geometry bridge. It intentionally contains
no assembled DTA pipeline because the previous biological state did not provide
registered protein-specific affinity increment.

## Repository Boundaries

- `theory/FINAL_FROZEN_THEORY/`: authoritative read-only mathematics.
- `model/`: verified Band/CSMO/law primitives and P1B geometry components.
- `scripts/`: passed data, structure, geometry, and release-governance workflows.
- `research/e0_identifiability/`: synthetic identifiability evidence, E0R2,
  the negative T-DIR-P0 pilot, and the passing T-BASIS-R0 fixed-radial study;
  no production model.
- `report/`: current split protocol and immutable PASS evidence only.
- `history.md`: complete failure ledger and deleted-artifact record.

Failed experimental implementations and duplicate reports were consolidated and
deleted. Start with `EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md` for the full
theory-to-biology assessment, evidence map, and stop rules.

## Retained Evidence

- P0 canonical data contracts: PASS.
- P1A governed open holo corpus: PASS.
- P1B partner-specific contact/distance geometry: PASS.
- D0-C release-pinned ChEMBL37 Ki/Kd corpus: PASS.
- D1 homology/document closure: PASS.
- E0 real affinity source gate: NOT RUN.
- T-DIR-P0 lightweight pilot: learnability signal not observed.
- T-BASIS-R0 fixed radial basis partner recoverability: PASS in research.
- E-AFF-P0 population-shared radial affinity direction: NOT OBSERVED.
- E-AFF-H0A task-local radial headroom: OBSERVED; partner specificity below Gate.
- E-AFF-H0C support-matched interaction residual: NOT OBSERVED.
- E-AFF-X0 crossed ChEMBL census: INSUFFICIENT INDEPENDENT COMPONENTS; STOP.
- Formal typed-interaction T and P2-P4: FROZEN.
- Recipient-label reads: `0`.

## Read First

1. `EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md`
2. `task.md`
3. `experiment.md`
4. `history.md`
5. `report/CURRENT_DATA_SPLIT_PROTOCOL.md`

## Verification

Environment: `D:\anaconda\envs\drug`, PyTorch `2.6.0+cu124`, CUDA RTX 4060.

```powershell
D:\anaconda\envs\drug\python.exe -m pytest -q
```

Current result: `70 passed in 77.02s`, including E-AFF/L0/X0-B contract tests. T-DIR,
T-BASIS and E-AFF stages remain research-only.

## Data Availability

The Git repository tracks code, frozen theory, provenance manifests, compact
PASS evidence, and registered research artifacts. Large third-party releases,
raw benchmark labels, embedding banks, downloaded tools, and model caches are
not redistributed. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for the
exact boundary and upstream provenance.
