# MetaSieve Current Tasks

Updated: 2026-08-08

## Current Decision

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
P1C/P1R1/P1R2A/P1R2B0/P1R2B1 FAIL OR CLOSED
E0 SYNTHETIC PRE-GATE FAIL-CLOSED
E0R1 PINV WITNESS PASS; E0R2 SYNTHETIC NUMERICAL CLOSURE PASS
T-DIR-P0 LIGHTWEIGHT PILOT: LEARNABILITY SIGNAL NOT OBSERVED
T-BASIS-R0 FIXED RADIAL PARTNER RECOVERABILITY PASS
E-AFF-P0 POPULATION-SHARED DIRECTION NOT OBSERVED
E-AFF-H0A TASK-LOCAL HEADROOM OBSERVED; PARTNER SPECIFICITY BELOW GATE
E-AFF-H0C SUPPORT-MATCHED RADIAL INTERACTION RESIDUAL NOT OBSERVED
E-AFF-X0 CROSSED CENSUS: STOP SOURCE INTERACTION UNDERDETERMINED
E-AFF-R0 READOUT SCOPE AUDIT: RANK METRIC BLIND TO TASK-LEVEL LOCATION
E-AFF-X0-FEAS X0 UNIT REQUIREMENT UNATTAINABLE BY CONSTRUCTION
E-AFF-X0-B CONDITIONAL DESIGN SUPPORT PENDING RHO
E-AFF-L0 NOT RUN: NUMERICAL PRECONDITION FAILED
E-AFF-L0R NOT RUN: POSITIVE CONTROL ABSENT (CLAIM A UNTESTED)
PKIS/KLIFS F6I COMPONENT DECOMPOSITION: DEVELOPMENT SIGNAL OBSERVED
F6I REGISTERED TOTAL GATE: NOT ADMISSIBLE
RESEARCH LAW BRIDGE INVARIANTS PASS; NOT THE FROZEN PRODUCTION OPERATOR
FRESH ENDPOINT-CONSISTENT EXTERNAL ADMISSION: NOT REGISTERED
E-AFF-H0B/X1/X2/RFSA/ANGULAR BASIS/P2-P4 NOT AUTHORIZED OR FROZEN
```

## Two Separate Estimands

**Claim A** — does the correct protein provide affinity-location information
beyond population, ligand-only and protein-sequence-only baselines? Addressed by
`E-AFF-L0`.

**Claim B** — does a non-additive protein-by-ligand affinity interaction exist?
Addressed by `E-AFF-X0-B` and then `X1`.

**L0 and X1 are different estimands.** Neither stage replaces or authorizes the
other. An L0 pass would not be evidence of interaction; an X1 pass would not be
evidence of location.

## Scope Of The Historical Affinity Negatives

The P1C/P1R\*/E-AFF-P0/H0A/H0C evaluations used within-task concordance or
another rank-based metric. Those experiments remain valid and their verdicts are
unchanged. Their negative conclusions are limited to **within-task ranking
information**, because within-task concordance is invariant to task-wise
prediction shifts and positive rescaling. They did not test whether the correct
protein supplies transferable affinity-location information.

## Newly Registered And Recorded

- [x] `E-AFF-R0_READOUT_SCOPE_AUDIT` — readout diagnosis; no affinity label read.
- [x] `E-AFF-X0-FEAS` — unit feasibility audit of the frozen X0 requirement.
- [x] `E-AFF-X0-B` — crossed design re-registration under the cell-disjoint unit.
- [x] `E-AFF-L0` preregistration — Claim A location Gate, registered.
- [x] `E-AFF-L0` execution — ran and failed closed:
  `L0_NOT_RUN_NUMERICAL_PRECONDITION_FAILED`. Claim A remains untested.
- [x] `E-AFF-L0R` — corrected Gate on a fresh 78-component panel; both L0 defects
  fixed, but the registered positive control failed:
  `L0R_NOT_RUN_POSITIVE_CONTROL_ABSENT`. Claim A remains untested.

## Active Research Boundary

R0 established that within-task concordance is exactly invariant to per-task
affinity shift and positive rescaling: maximum deviation `0.0` across prediction
shift, prediction rescale, label shift and label rescale. A simulated predictor
holding a task's affinity level perfectly and nothing else scores exactly
`0.5000` at every variance share tested, up to `0.985`. H0C additionally removed
that channel upstream and shared a correct-protein support-fitted nuisance with
the deranged arm.

X0-FEAS established that the closure-component universe of the governed corpus
is `245` in total while the frozen requirement is `245` per endpoint, that the
best ceiling over every governed population is Ki `97` and Kd `56`, and that no
panel spans more than one closure component, so crossing can never create an X0
unit. X0's verdict is retained unchanged and annotated as specification-induced.

X0-B re-registered the unit as the cell-disjoint rectangle with a design-effect
effective sample size, leaving the effect size, alpha, power, `245` requirement
and `+0.03` affinity margins unchanged. Ki packs `11,168` units in `36`
clusters; Kd packs `1,041` in `12`. At `rho = 1` the design reduces exactly to
X0's `36` and `12`. Breakeven thresholds are `rho* = 0.0915` (Ki) and `0.0164`
(Kd). This is conditional design support, not evidence of interaction.

## Success Criteria For Any Future Admission

Code may return from `research/` to `model/` only after a separately registered
Gate establishes both:

1. correct protein improves over ligand-only by the frozen margin; and
2. correct protein improves over deranged protein with a positive lower bound.

The admitted object must also define explicit biological coordinates for
`z(S,Q,gamma)` without changing the frozen operator.

An L0 pass would establish at most
`PROTEIN_SPECIFIC_AFFINITY_LOCATION_IDENTIFIED_IN_SOURCE`. It would not
authorize biological `z` admission, `model/` promotion, DAVIS or recipient
evaluation, X1, angular or many-body bases, RFSA, CSMO/Band/theory
modification, P2-P4, an end-to-end DTA claim, or a physical binding-free-energy
interpretation. Production or biological-`z` admission still requires a
separately registered independent-source replication and a sealed novel-target
transfer Gate.

## New PKIS/F6I Research Intake

The supplied `pkis_mechanism_pilot` and `section_operator_pilot` packages have
been integrated under `research/` only. Their isolated suite passes `39` tests.
The evidence supports a component-level decomposition on consumed kinase
activity panels: a protein-dependent zero-shot surface plus a bounded
one-dimensional support-location statistic. It does not establish Ki/Kd
affinity energetics, an untouched external result, or production admission.

The registered F6I total verdict remains `F6I_COMPONENTS_NOT_ADMISSIBLE`; it is
not rewritten as PASS. The seven-point categorical `law_bridge.py` is also not
type-equivalent to the frozen 33-point CDF-band operator in `model/`. See
`research/IDENTIFIABILITY_RESOLUTION_INTAKE_AUDIT.md` for the accepted claims,
provenance limitation and next evidence boundary.

## Verification

The imported research package regression is `39 passed`; the full repository
regression is `70 passed` in the `drug` environment. R0, X0-FEAS and X0-B carry
independent audits. R0 read no affinity label; X0-FEAS and X0-B selected
zero affinity value fields. DAVIS and recipient reads remain zero. All evaluated
panels are development evidence and must not be reused as untouched validation.
