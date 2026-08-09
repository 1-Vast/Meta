# MetaSieve current research status

Updated: 2026-08-10.

## Current verdict

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
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

## Phase 2B terminal verdict — the current stopping point

```text
PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED
  -> preregister one repair of the Phase 2B OPTIMIZATION contract
```

The Phase 2A verdict authorized one ligand-conditioned residue residual head.
Its preregistration was audited before use, found to contain eleven design
defects, and superseded by
`PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md` (`5e6688f6…`, committed `b9753db`
**before** any implementation existed). The superseded document
(`ae6d1a01…`) was never executed and is kept byte-identical.

The repaired contract then passed every artifact and numerical precondition —
10,568 parameters exactly, atom-permutation invariance at 0.0, prior
cancellation at `2.05e-15`, projection orthogonality at `6.19e-15`, zero
component and ligand-graph overlap, census matching the registration exactly —
and **stopped at its own synthetic trainability control**: `AP_bidir = 0.3577`
against a preregistered `>= 0.50`.

The real-label run was not executed and gates `R1`–`R6` were not scored, so
**no biological conclusion is drawn**. Diagnostics localise the shortfall to the
registered optimization budget (teacher scores 0.99971 on its own labels;
in-sample 0.3654 vs held-out 0.3577, so no generalisation gap; learned field
correlates with the teacher field at 0.754 median). The a priori `0.50`
threshold may itself have been mis-calibrated; neither it nor the budget may be
adjusted now, because both were frozen and the synthetic holdout has been seen.

## Phase 2A terminal verdict

```text
LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
  -> preregister one ligand-conditioned residue residual head
```

Phase 2A was audit-only over the sealed Phase 1 predictions and the MONN labels.
Nothing was trained, no affinity value was read, no frozen surface was modified,
and nothing was committed. Registration SHA-256 `4e01401d…` with amendments
01–03, each frozen before the phase it governs.

## The correction Phase 2A forced

Phase 1 reported 92.5% wrong-ligand retention of B5's residue AP and concluded
"mostly generic pocket". That was a statement about **B5** and was being carried
as if it were also a statement about the **corpus**. It is not. The Phase 1
control used an arbitrary foreign ligand; Phase 2A used a real alternative
ligand of the same exact construct, against a noise floor measured from the data
(same construct, same ligand, different crystal).

| component-macro Jaccard of residue masks | value |
|---|---:|
| replicate — same ligand, different crystal | 0.6361 |
| alternative scaffold-distinct ligand | 0.4165 |
| **ΔJ, paired over 292 closure components** | **+0.2580 [LCB +0.2344]** |

Registered minimum meaningful effect: 0.05. Chemistry association:
Spearman ρ = **+0.3221 [LCB +0.2987]** with a ligand-permutation control at
median p = 0.03. 80.4% of scaffold-distinct pairs change the mask meaningfully.
About 44% of the alternative-ligand difference is ligand-attributable; the rest
sits at replicate noise level.

## Where B5's ligand dependence lives

| arm | full | residue marg. | atom marg. | additive | coupling |
|---|---:|---:|---:|---:|---:|
| **B5** | **0.06975** | 0.04045 | 0.00514 | 0.03983 | **0.01133** |
| B4 | 0.02323 | 0.01619 | 0.00550 | 0.01510 | 0.00637 |
| BX5 wrong ligand | 0.01969 | 0.03595 | 0.00326 | 0.02769 | 0.00346 |
| BP5 wrong protein | 0.00464 | 0.00461 | 0.00492 | 0.00563 | 0.00355 |
| BL ligand-only | 0.00573 | 0.00330 | 0.00573 | 0.00573 | 0.00305 |

A wrong ligand retains 89% of B5's residue marginal but only 31% of its coupling
term. B5 uses the ligand only through the pair term, and that term is small:
`+0.0060 [LCB +0.0046]` over its degree-preserving rewiring null and
`+0.0079 [+0.0062]` over the wrong-ligand arm, both clearly above zero and both
below the preregistered 0.01 margin. The teacher's own edge coupling also fails
(median `z = +0.413`, threshold 2.0), reproducing I-2 under a stricter rewiring
specification with zero degree-preservation violations.

## The number that sets the priority

The well-posed label-fitted additive ceiling is **0.3889**. B5 reaches
**17.9%** of it. Its residue marginal reaches 19.8% of the true residue-margin
ceiling of 0.2043. **The bottleneck is the residue marginal, not the coupling.**

## Label semantics

Not ambiguous. Water-mediated edges are 8.2% (threshold 20%), and removing them
strengthens the teacher result (ΔJ 0.258 → 0.278). A local dense-distance
comparator was built on 1,909 complexes: 88.1% of PLIP positives lie within
5.0 Å of a ligand heavy atom. A second interaction-annotation tool does not
exist locally, so that one comparison remains **UNRESOLVED**. PU learning and a
soft teacher remain unauthorized.

## Confirmation boundary

Unchanged: MONN cannot supply a time-forward confirmation cohort under the
frozen 2019 cutoff (two qualifying entries; none at 2024). A document-closed
cohort supports development analysis only.

## Next eligible stage

Preregister **one** repair of the Phase 2B **optimization contract**. Before any
real label is read it must fix the budget and sampler caps from a measured
synthetic scaling curve, derive the acceptance threshold from that curve rather
than by assertion, use a fresh synthetic teacher seed (`20260905` has been
observed), and leave the architecture, projection, metric, controls and gates
`R1`–`R6` byte-identical to `PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md`.

Adding capacity, another PLM, attention, a GNN, a geometry branch or a
typed-interaction branch is **not** an admissible response to a
synthetic-precondition failure and is not proposed.

**Not registered, not authorized.** Nothing has been implemented for it.

## Frozen boundary

- real ChEMBL/BindingDB affinity training;
- DAVIS, KIBA and recipient labels;
- independent confirmation scoring;
- new PLM, attention stack, geometry branch, typed-interaction branch, affinity
  head, PU loss, knowledge graph or parallel module;
- few-shot adaptation and section claims;
- biological statistic admission to `z`;
- CSMO, Band, mesh, positive ridge, and frozen theory;
- P2-P4.

## Canonical records

0. `report/s7_l2b_r0r/PHASE2B_REPORT.md`
1. `report/s7_l2b_r0r/PHASE2A_SYNTHESIS.md`
2. `report/s7_l2b_r0r/PHASE2A_VERDICT.json`
3. `report/s7_l2b_r0r/PHASE1_EVIDENCE_CONSOLIDATION_AND_PHASE2A_TRIAGE.md`
4. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
5. `project_state.json`
6. `task.md`
7. `history.md`

Failed and withdrawn stages remain only in the evidence ledger, `history.md`,
and recoverable Git history; they are not current execution instructions.
