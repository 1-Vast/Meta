# Current task

## Current state

```text
PHASE2A_AUDIT_COMPLETE
TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED_IN_LABELS
B5_RESIDUE_MARGINAL_IS_GENERIC_POCKET
B5_LIGAND_DEPENDENCE_CONFINED_TO_THE_COUPLING_TERM
TEACHER_EDGE_COUPLING_NOT_IDENTIFIED
EXACT_RESIDUE_ATOM_COUPLING_NOT_IDENTIFIED
LABEL_SEMANTICS_NOT_AMBIGUOUS
PHASE2B_RESIDUE_RESIDUAL_NOT_RUN_ON_REAL_LABELS
AFFINITY_DIRECTION_NOT_TESTED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
```

## Phase 2B — stopped at a fail-closed precondition

```text
TERMINAL VERDICT   PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED
REAL TRAINING      NOT EXECUTED
GATES R1-R6        NOT SCORED
```

The registered synthetic trainability control returned `AP_bidir = 0.3577`
against a preregistered requirement of `>= 0.50`. The threshold was not lowered,
nothing was tuned against the synthetic holdout, and no second seed was tried.

The contract itself passed everything else. `PHASE2B_CONTRACT_PASS` on all 14
preflight items: 10,568 parameters exactly, `g(L)` atom-permutation invariant to
0.0, ligand-order swap sign-exact to 0.0, prior cancellation `2.05e-15`,
projection orthogonality `6.19e-15`, zero train/held-out component overlap, zero
held-out ligand-graph overlap. Census matched the registration: 226,765 training
and 46,818 held-out A eligible pairs; foreign-pair control coverage 1.000;
derangement with 0 fixed points.

Failure localisation, from gauge-invariant diagnostics:

| diagnostic | value | reading |
|---|---:|---|
| teacher on its own labels | 0.99971 | evaluation code sound |
| in-sample vs held-out AP | 0.3654 / 0.3577 | no generalisation gap |
| field correlation with the teacher | 0.754 median | class is being fitted |
| parameter movement `U` / `V` | 1.808 / 0.426 | the head trained |

So the shortfall is the registered **optimization budget** — and possibly the a
priori `0.50` threshold, which was set without a calibration curve. Neither can
be adjusted now: both were frozen and the synthetic holdout has been seen.

**No biological conclusion is permitted from this run.** It does not show that
the frozen sequence + 2-D representation lacks a ligand-conditioned residue
correction.

## Next eligible work — requires separate authorization

Preregister **one** repair of the Phase 2B **optimization contract**. Before any
real label is read it must:

1. fix the budget and sampler caps from a **measured synthetic scaling curve**;
2. derive the synthetic acceptance threshold from that curve, with the
   derivation stated;
3. use a **fresh synthetic teacher seed** (`20260905` has been observed);
4. leave the architecture, projection, metric, controls and gates `R1`–`R6`
   byte-identical to `PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md`.

Adding capacity, another PLM, attention, a GNN, a geometry branch or a
typed-interaction branch is **not** an admissible response and is not proposed.

## Frozen

- ChEMBL/BindingDB affinity training; DAVIS, KIBA and recipient labels;
- the sealed independent structural confirmation;
- any new PLM, GNN, attention stack, geometry branch, typed-interaction branch,
  affinity head, PU loss, knowledge graph or parallel SSL module;
- ESM2 fine-tuning, hyperparameter search, seed selection;
- few-shot section adaptation and any `k`-shot claim;
- biological `z` admission;
- CSMO, Band, mesh, and frozen theory;
- P2-P4; pushing to GitHub.

## Evidence

- `report/s7_l2b_r0r/PHASE2B_REPORT.md`
- `report/s7_l2b_r0r/PHASE2B_GATE.json`
- `report/s7_l2b_r0r/PHASE2B_SYNTHETIC_AUDIT.json`
- `report/s7_l2b_r0r/PHASE2B_INPUT_MANIFEST.json`
- `report/s7_l2b_r0r/PHASE2B_CONTROL_MANIFEST.json`
- `research/s7_l2b_r0r/PHASE2B_DESIGN_AUDIT.md`
- `report/s7_l2b_r0r/PHASE2A_SYNTHESIS.md`
- `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
- `history.md`
