# MetaSieve current research status

Updated: 2026-08-10.

## Current verdict

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_AND_PAIR_COMPATIBILITY_IDENTIFIED
FROZEN_ESM2_EXACT_RESIDUE_LOCALISATION_PASS_IN_DEVELOPMENT
TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED
SYNTHETIC_BINARY_ORDINAL_ESTIMATOR_TRAINABLE
LIGAND_MEAN_POOLING_COLLAPSE_MEASURED
GRAPH_AWARE_LIGAND_INCREMENT_REAL_BUT_NOT_LIGAND_SPECIFIC
REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

## Earliest failed boundary

Phase 2A proved that same-construct scaffold-distinct ligands change the MONN
residue masks beyond the replicate noise floor. S2R repaired the synthetic
optimization defect and passed a sealed synthetic seed (`AP_bidir = 0.6620`).
S3R transferred that estimator to real labels and failed every Gate, scoring
`0.03588` against chance `0.02547`, with participation and replay intact. The
failure was therefore scoped to the measurement basis rather than the pipeline.

S4R tested that scoping directly. A label-blind audit first confirmed the
mean-pooled 41-D ligand basis is collapsed: pair-difference effective rank
`6.183`, 687 distinct ligand graphs sharing a bit-identical vector, and `85.2%`
of the difference-norm variance explained by heavy-atom count alone. Two
constitutional isomers with identical atom composition and different
connectivity map to the *same* vector. A frozen radius-1 Morgan per-heavy-atom
statistic at `d = 128` raises that effective rank to `20.93` and places `35.5%`
of its difference energy beyond anything the baseline can linearly express.

S4R then swapped only that statistic into the S3R stage. On the same 46,818
pairs and 112 components, with a bit-exact reproduction of the S3R baseline as
an anchor:

```text
candidate                0.046856     baseline41 (= S3R)   0.035880
foreign ligand pair      0.046212     chance               0.025472
R1 candidate - chance   +0.021384 [LCB +0.016064]  needs +0.05   FAIL
R3 candidate - foreign  +0.000644 [LCB -0.009226]  needs +0.03   FAIL
C1 candidate - baseline +0.010976 [LCB +0.004939]  non-gating
```

```text
REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED
```

The representation was a real bottleneck and removing it doubled the
above-chance gain; the candidate also beats the capacity-matched permuted-label
learner, which S3R did not. But the surviving signal is invariant to which
ligands are supplied, so it is a construct-level residue-change prior, not
ligand-conditioned residue selection.

## Current boundary

No active training stage is authorized. The registered stopping rule closes the
pose-free ligand representation repair route, including re-running S4R at a
larger vocabulary or radius. The remaining question — whether any pose-free
sequence-plus-2D estimand can bind a ligand substructure to a residue context
without geometric correspondence — needs a separately governed information
stage, not another repair.

Heldout-B, R6, affinity values, few-shot sectioning, biological `z`, CSMO/Band
and the frozen law operator remain unopened and unchanged. Heldout-B was not
even created by S4R.

## Canonical evidence

1. `report/s7_l2b_r0r/PHASE2B_S4R_GATE.json`
2. `report/s7_l2b_r0r/PHASE2B_S4R_EVIDENCE_CONSOLIDATION.md`
3. `report/s7_l2b_r0r/PHASE2B_S4R_REPRESENTATION_AUDIT.json`
4. `report/s7_l2b_r0r/PHASE2B_S3R_GATE.json`
5. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
6. `history.md`
