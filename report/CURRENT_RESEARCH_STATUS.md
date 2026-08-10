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
LIGAND_STEERING_PRESENT_BUT_BIOLOGICALLY_MISDIRECTED
CONDITIONAL_ESTIMAND_REPAIR_ROUTE_CLOSED
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

## Why, and the falsified explanation

S5D trained nothing and reused the frozen S4R checkpoints to ask why R3 failed.
It registered the obvious mechanism — that the estimator collapses ligand
differences onto one residue direction per protein — and **falsified it**. The
top principal energy fraction of the candidate's residue fields is `0.4793`
against a data-side upper bound of `0.4550`, an excess of `0.0138` where the
rule needed `0.80` and `+0.10`, and the median cosine between a pair's true and
foreign fields is `0.4487`. The estimator does steer on the ligand, and it
steers more than the mean-pooled baseline does.

S5D then measured the symmetric-difference conditional estimand, which cancels
pocket membership exactly by restricting each comparison to residues that
changed. On 40,157 pairs across 107 components it found nothing:

```text
candidate 0.655030   foreign 0.655470   chance 0.643744   baseline41 0.638830
E1 candidate - chance   +0.011285 [LCB -0.007749]  needs +0.05  FAIL
E2 candidate - foreign  -0.000440 [LCB -0.021814]  needs +0.03  FAIL
```

```text
LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED
```

Ligand information is neither lost upstream nor diluted by the metric. It
reaches the residue field, rotates it by a large angle, and the direction it
chooses is unrelated to which residues gained or lost contact.

## Current boundary

No active training stage is authorized. Both pose-free repair routes are closed:
the representation route by S4R, including re-running at a larger vocabulary or
radius, and the estimand route by S5D, which also forbids a fourth estimand
variant on heldout-A. Heldout-A has now been consumed three times and every
number from it is development evidence.

The remaining hypothesis is that the missing ingredient is correspondence —
which ligand substructure sits against which residue — and that a pose-free
sequence-plus-2D estimand has no channel to supply it. Testing that needs a
separately governed information stage about geometry, not another repair.

Heldout-B, R6, affinity values, few-shot sectioning, biological `z`, CSMO/Band
and the frozen law operator remain unopened and unchanged. Heldout-B was not
even created by S4R or S5D.

## Canonical evidence

1. `report/s7_l2b_r0r/PHASE2B_S5D_GATE.json`
2. `report/s7_l2b_r0r/PHASE2B_S5D_EVIDENCE_CONSOLIDATION.md`
3. `report/s7_l2b_r0r/PHASE2B_S4R_GATE.json`
4. `report/s7_l2b_r0r/PHASE2B_S4R_EVIDENCE_CONSOLIDATION.md`
5. `report/s7_l2b_r0r/PHASE2B_S4R_REPRESENTATION_AUDIT.json`
6. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
7. `history.md`
