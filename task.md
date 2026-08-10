# Current task

## Current state

```text
PHASE2A_TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED
S2R_BINARY_ORDINAL_TRAINABILITY_PASS
S3R_REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
S4R_A_LIGAND_MEAN_POOLING_COLLAPSE_MEASURED
S4R_GRAPH_AWARE_INCREMENT_REAL_BUT_NOT_LIGAND_SPECIFIC
POSE_FREE_LIGAND_REPRESENTATION_REPAIR_ROUTE_CLOSED
AFFINITY_DIRECTION_NOT_TESTED
K_SHOT_SECTION_NOT_TESTED
BIOLOGICAL_Z_NOT_ADMITTED
```

S4R executed the single authorized single-axis repair in the `drug`
environment. It replaced only the mean-pooled 41-D ligand statistic with a
frozen radius-1 Morgan per-heavy-atom statistic (`d = 128`, `163,840`
parameters) and held the protein states, gauge, direct-`W` estimator, loss,
sampler, closure split, 210-update stream, seeds, control maps and R1-R5
margins byte-identical to S3R.

The representation was not the excuse and the change was not cosmetic:

```text
candidate                0.046856
baseline41 (= S3R)       0.035880
chance                   0.025472
C1 candidate - baseline  +0.010976  [LCB +0.004939]
R5 candidate - permuted  +0.010563  [LCB +0.003880]   S3R had -0.001245
```

But the gain is not ligand-conditioned. Foreign ligand pairs cost only
`+0.000644 [LCB -0.009226]` and the within-construct chemistry shuffle scores
`0.051322`, above the candidate. R1 needs `+0.05` and observed `+0.021384`.

Terminal verdict: `REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED`.

## Next research decision

No experiment is authorized, and no repair of this estimand is eligible. The
registered stopping rule closes the pose-free representation route, including
re-running S4R at `d = 256` or `d = 512`.

Two stages of evidence now bound the question. Phase 2A proved the labels carry
same-construct ligand conditionality; S2R proved the estimator trains; S4R
proved that removing the ligand-representation bottleneck doubles the
above-chance signal yet leaves it invariant to which ligands are supplied. The
open question is no longer representational richness but whether any pose-free
sequence-plus-2D estimand can bind a ligand substructure to a residue context
without geometric correspondence. Answering it requires a separately governed
information stage with its own preregistration, not another repair here.

## Frozen

- heldout-B and R6 amplitude/B5 integration;
- ChEMBL/BindingDB affinity, DAVIS, KIBA and recipient labels;
- larger PLM, second protein encoder, attention stack, parallel branch,
  geometry/pose branch, typed interaction branch, KG, PU loss or affinity head;
- larger ligand vocabulary or radius as a rescue of S4R;
- few-shot sectioning and biological `z` admission;
- CSMO, Band, mesh and `A(F,z)=K(B(z)F(z))`.

## Read first

1. `report/s7_l2b_r0r/PHASE2B_S4R_EVIDENCE_CONSOLIDATION.md`
2. `report/s7_l2b_r0r/PHASE2B_S4R_GATE.json`
3. `report/s7_l2b_r0r/PHASE2B_S4R_REPRESENTATION_AUDIT.md`
4. `report/CURRENT_RESEARCH_STATUS.md`
5. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
