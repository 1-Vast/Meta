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
AFFINITY_DIRECTION_NOT_TESTED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
```

Phase 2A was audit-only: nothing trained, no affinity value read, no frozen
surface touched, nothing committed.

```text
TERMINAL VERDICT   LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
NEXT ACTION        preregister one ligand-conditioned residue residual head
```

## What Phase 2A established

The labels **are** ligand-conditioned at the residue level, measured against a
noise floor taken from the data itself (same construct, same ligand, different
crystal):

| quantity | value |
|---|---:|
| replicate Jaccard (noise floor) | 0.636 |
| alternative-ligand Jaccard | 0.416 |
| `T1` ΔJ, paired over 292 components | +0.258 [LCB +0.234] |
| `T5` Spearman ρ(mask dissimilarity, chemical distance) | +0.322 [LCB +0.299] |
| `T7` pairs with a meaningful residue change | 80.4% |

B5 is not. Decomposing the sealed logits, a wrong ligand retains **89%** of B5's
residue marginal but only **31%** of its coupling term — B5's ligand dependence
lives entirely in the pair term, and that term is small:

| contrast | Δ | LCB95 | margin | |
|---|---:|---:|---:|---|
| B5 coupling − degree-preserving null | +0.0060 | +0.0046 | 0.01 | FAIL |
| B5 coupling − wrong-ligand coupling | +0.0079 | +0.0062 | 0.01 | FAIL |

Both are clearly above zero and both are below the preregistered practical
margin. The teacher's own edge coupling also fails (median `z = +0.413` against
a threshold of 2.0), reproducing I-2 under a stricter rewiring specification.

The headroom number that sets the priority: B5 reaches **17.9%** of the
well-posed label-fitted additive ceiling (0.0698 of 0.389). The bottleneck is
the residue marginal, not the coupling.

## Next eligible work

`research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md`
(SHA-256 `ae6d1a01…`) is written and frozen. It registers one head:

```text
logit p_r(P, L) = b_r(P) + delta_r(P, L)
```

with `b_r(P)` the frozen B5 residue-marginal prior, `delta_r` one low-rank
bilinear residual (`K <= 8`) over existing frozen states, projected away from
the constant, the pocket prior and ligand-only directions, and supervised only
by the same-protein ligand differential on symmetric-difference residues.

Gates `D1`–`D5` and a fail-closed module-participation audit are frozen in that
document. **It is registered, not authorized.** No Phase 2B code exists and no
run has occurred.

## Frozen

- ChEMBL/BindingDB affinity training;
- DAVIS, KIBA and recipient labels;
- independent confirmation scoring;
- new PLM, attention stack, geometry branch, typed-interaction branch, affinity
  head, PU loss, knowledge graph or parallel module;
- few-shot section adaptation;
- biological `z` admission;
- CSMO, Band, mesh, and frozen theory;
- P2-P4.

## Evidence

- `report/s7_l2b_r0r/PHASE2A_SYNTHESIS.md`
- `report/s7_l2b_r0r/PHASE2A_VERDICT.json`
- `report/s7_l2b_r0r/PHASE2A_TEACHER_CONDITIONALITY.json`
- `report/s7_l2b_r0r/PHASE2A_MARGINAL_COUPLING_AUDIT.json`
- `report/s7_l2b_r0r/PHASE2A_LABEL_SEMANTICS.json`
- `report/s7_l2b_r0r/PHASE2A_DATA_IDENTIFIABILITY_CENSUS.json`
- `report/s7_l2b_r0r/PHASE2A_INPUT_MANIFEST.json`
- `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
- `history.md`
