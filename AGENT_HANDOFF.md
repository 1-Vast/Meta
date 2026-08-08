# MetaSieve Agent Handoff

Updated: 2026-08-08

## State

The repository has been consolidated around verified evidence. There is no
validated assembled DTA model.

```text
KINASE_PANEL_COMPONENT_IDENTIFIABILITY_OBSERVED_IN_DEVELOPMENT
F6I_TOTAL_GATE_NOT_ADMISSIBLE
FRESH_ENDPOINT_CONSISTENT_EXTERNAL_ADMISSION_NOT_RUN
WITHIN_TASK_RANKING_DIRECTION_NOT_IDENTIFIED
PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_YET_TESTED
CROSSED_INTERACTION_EXISTENCE_NOT_YET_TESTED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

```text
P0/P1A/P1B/D0-C/D1 PASS
E0R2 SYNTHETIC OBJECTIVE/DESIGN/SOLVER IDENTIFIED
T-DIR-P0 PILOT LEARNABILITY SIGNAL NOT OBSERVED
T-BASIS-R0 FIXED RADIAL RECOVERABILITY IDENTIFIED
E-AFF-P0 SHARED DIRECTION NOT OBSERVED
E-AFF-H0A TASK-LOCAL HEADROOM WITHOUT PARTNER SPECIFICITY
E-AFF-H0C CENTERED RADIAL INTERACTION NOT OBSERVED
E-AFF-X0 STOP SOURCE INTERACTION UNDERDETERMINED
E-AFF-R0 READOUT SCOPE AUDIT COMPLETE
E-AFF-X0-FEAS X0 UNIT REQUIREMENT UNATTAINABLE BY CONSTRUCTION
E-AFF-X0-B CONDITIONAL DESIGN SUPPORT PENDING RHO
E-AFF-L0 NOT RUN: NUMERICAL PRECONDITION FAILED
E-AFF-L0R NOT RUN: POSITIVE CONTROL ABSENT
PKIS/F6I COMPONENT DECOMPOSITION OBSERVED ON CONSUMED DEVELOPMENT PANELS
F6I TOTAL GATE NOT ADMISSIBLE
RESEARCH LAW BRIDGE TESTED; NOT PRODUCTION OPERATOR EQUIVALENCE
X1/X2, ANGULAR BASIS, RFSA, DAVIS, P2-P4 FROZEN
RECIPIENT_LABEL_READS=0
DAVIS_LABEL_READS=0
```

## Two Separate Claims

The affinity question splits into two estimands that must not be conflated.

- **Claim A (E-AFF-L0).** Does the correct protein provide affinity-*location*
  information beyond population, ligand-only and protein-sequence-only
  baselines?
- **Claim B (E-AFF-X0-B, then X1).** Does a non-additive protein-by-ligand
  affinity *interaction* exist?

**L0 and X1 are different estimands.** Neither replaces, implies nor authorizes
the other. An L0 result says nothing about interaction; an X1 result says
nothing about location.

## Scope Of The Historical Negatives

P1C, P1R1, P1R2A, P1R2B0, P1R2B1, E-AFF-P0, H0A and H0C used within-task
concordance or another rank-based metric. Those verdicts and metrics stand
unchanged. E-AFF-R0 established that this readout is **exactly** invariant to
per-task affinity shift and positive rescaling, so those results are evidence
about within-task ranking information only. They were never evidence about
protein-specific affinity location, which has not yet been tested.

No historical FAIL, CLOSED, NOT-RUN or negative result has been relabelled.

## The X0 Stop Is Specification-Induced

X0's verdict `STOP_SOURCE_INTERACTION_UNDERDETERMINED` is retained unchanged.
X0-FEAS subsequently established that X0's registered requirement of 245
effective components per endpoint exceeded the structural ceiling of its own
independence unit at every governed population — best ceiling Ki `97`, Kd `56`.
The stop therefore reflects the unit definition rather than a measurement of how
crossed ChEMBL37 is. This is an annotation of scope, not a reclassification.

X0-B re-registered the unit as the cell-disjoint rectangle with a design-effect
effective sample size, leaving the frozen effect size, alpha, power, 245
requirement and `+0.03` affinity margins unchanged. Its verdict is **conditional
design support**, not evidence that affinity interaction exists.

## Read Order

1. `EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md`
2. `task.md`
3. `experiment.md`
4. `research/e0_identifiability/THEORY_BIOLOGY_INTEGRATION.md`
5. `research/IDENTIFIABILITY_RESOLUTION_INTAKE_AUDIT.md`
6. `history.md`

## Code Boundaries

- `model/` has verified mathematical primitives and P1B geometry only.
- `scripts/` has passed data, geometry, and governance workflows only.
- `research/e0_identifiability/` contains terminal synthetic, structural and
  source-affinity research evidence plus the registered, unrun L0. It is not
  production code.
- `research/pkis_mechanism_pilot/` and `research/section_operator_pilot/`
  contain imported kinase-panel development evidence. Their F6I total Gate is
  not admissible; neither package may be imported by production code.
- `report/` retains only current protocol and PASS evidence.

Do not recreate deleted failed implementations from prose. Do not run X1/X2,
access DAVIS or recipient labels, reuse any consumed evaluation panel as
untouched validation, build angular/reference-state potentials, start RFSA, or
modify the frozen mathematical operator without a new explicit registration.

## Last Verified Facts

- P1B correct contact AUPRC `0.43885`; wrong-protein `0.05149`.
- D0/D1 retained `3,817` tasks, `697` targets, `253` closure components.
- E0R2 deterministic closure train RMSE is `3.188e-8`; synthetic only.
- T-BASIS-R0 test reconstruction/partner gain is `0.5312/0.1561`.
- E-AFF-H0C correct-minus-local/correct-minus-deranged is `-0.00391/+0.00152`.
- X0 found `36` Ki and `12` Kd effective components against a requirement of
  `245`; no affinity values were read.
- X0-FEAS: closure-component universe `245` total against `245` required per
  endpoint; best ceiling Ki `97`, Kd `56`; panels spanning more than one closure
  component `0`.
- X0-B: Ki `11,168` cell-disjoint units in `36` clusters, breakeven
  `rho* = 0.0915`; Kd `1,041` units in `12` clusters, `rho* = 0.0164`; at
  `rho = 1` the design reproduces X0's `36` and `12` exactly.
- R0: within-task concordance invariance deviation `0.0` on all four
  transforms; a perfect task-level predictor scores exactly `0.5000` at every
  variance share up to `0.985`.
- Imported PKIS/F6I research regression: `39 passed` in the `drug` environment.
- Full repository regression after intake: `70 passed`.
