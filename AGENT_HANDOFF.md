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
CROSSED_INTERACTION_EXISTENCE_IDENTIFIED_XP1_EXTERNAL_PANELS
INTERACTION_IS_LOW_RANK_r1_TO_r3_XP1
SUPPORT_IDENTIFIED_INTERACTION_SECTION_TRANSFERS_XP1
ZERO_SHOT_PROTEIN_FEATURE_INTERACTION_MAP_NOT_IDENTIFIED_XP1
UNSEEN_LIGAND_LOADING_NOT_TESTED
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
XP1-A CROSSED INTERACTION EXISTENCE IDENTIFIED ON TWO PLATFORMS
XP1-B SUPPORT-IDENTIFIED SECTION PASSES GATE; ZERO-SHOT SURFACE DOES NOT
XP1-C PDSP REPLICATION DIRECTIONAL, EFFECT BELOW NEGLIGIBILITY FLOOR
XP1-D REPRESENTATION FAILURE LOCALISED TO CROSS-GROUP TRANSFER
XP1-E TRUNCATION DESTRUCTIVE CONTROL PASSED
RECIPIENT_LABEL_READS=0
DAVIS_LABEL_READS=0
CHEMBL37_AFFINITY_VALUE_READS_IN_XP1=0
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

## XP1 Crossed-Panel Identification (2026-08-08)

Registered in `research/crossed_panel_identification/PREREG_XP1.md`, reported in
`report/crossed_panel_identification/XP1_RESEARCH_REPORT.md`. Three external
complete crossed panels acquired (Metz 2011, Klaeger 2017, NIMH PDSP). DAVIS,
PKIS2 and Anastassiadis excluded. No ChEMBL37 affinity value read.

```text
XP1 DECISION: OBJECTIVE_OR_PARAMETERIZATION_FAILURE
  secondary REPRESENTATION_FAILURE (zero-shot protein pathway only)
  DATA_IDENTIFIABILITY_FAILURE REJECTED
```

- Protein-by-ligand interaction is `59.6%` of panel affinity variance; `38%` of
  the residual is reproducible; geometry replicates at `r=0.885` across disjoint
  compound halves and `r=0.565` across an independent platform.
- Under strict kinase-group closure, a rank-1..3 section identified from `k=16`
  support labels reaches `R2_gamma = +0.160 [+0.109, +0.195]` and is
  target-specific: `+0.0695 [+0.0524, +0.0874]` against deranged support.
- Every zero-shot protein representation — ESM-2 t30 (production encoder),
  aligned KLIFS pocket, pocket physicochemistry, KLIFS conformational state,
  family, group, homolog kNN — is indistinguishable from a random-feature null at
  group closure, while the aligned pocket reaches `R2 = 0.52` leave-one-protein-out.
  This is `PARTNER_COMPATIBILITY` measured, not affinity direction.
- Destructive control: a synthetic additive panel with identical margins, noise
  and 40% truncation yields `Delta_specific = +0.00004 [-0.0005, +0.0011]`.

Candidate statistic `z_section = <u(L), vhat(S)>` (a bounded scalar, one CSMO
view) passes the registered Gate on the primary panel; PDSP replication is
directional only (`R2_gamma = 0.024`, below the preregistered floor). It is
**not** promoted to production `z`. Named continuation is `XP2`: the transpose,
holding out ligands, testing whether `u(L)` is predictable from chemistry.

## Read Order

1. `report/crossed_panel_identification/XP1_RESEARCH_REPORT.md`
2. `EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md`
3. `task.md`
4. `experiment.md`
5. `report/VERIFIED_EVIDENCE_SUMMARY.md`
6. `history.md`

## Code Boundaries

- `model/` has verified mathematical primitives and P1B geometry only.
- `scripts/` has passed data, geometry, and governance workflows only.
- `research/` contains the policy boundary plus the registered XP1 package in
  `research/crossed_panel_identification/`. Terminal pre-XP1 research
  implementations and artifacts were removed after their conclusions were
  consolidated in `history.md`.
- `model/component_statistic.py` retains only the gauge-separated algebra from
  F6I and is deliberately not exported or connected to the production state.
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
- Removed PKIS/F6I research package: isolated `39 passed`, recoverable at
  `8b7789e`.
- Consolidated repository regression: `73 passed` in the `drug` environment,
  re-verified after XP1.
- XP1 `BLK-METZ-60`: 704 compounds x 82 kinases, 34,764 uncensored cells;
  variance shares ligand `0.283`, protein `0.131`, interaction+noise `0.596`.
- XP1 reproducible interaction sd `0.442` log units, saturating near rank 12-20.
- XP1 group-closure arms: `A2` RMSE `0.8021`, `A4` `0.7352`, oracle `AO1`
  `0.6315`; `A3` best zero-shot `0.7913`.
- XP1 support-size identification curve `R2_gamma`: `0.052 / 0.092 / 0.160 /
  0.221 / 0.258` at `k = 4 / 8 / 16 / 32 / 64`.
- XP1 PDSP measurement noise: per-report `sigma = 0.714` log units over 2,490
  replicated cells.
