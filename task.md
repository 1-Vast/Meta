# Current task

## Current state

```text
PHASE2A_TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED
S2R_BINARY_ORDINAL_TRAINABILITY_PASS
S3R_REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
S4R_A_LIGAND_MEAN_POOLING_COLLAPSE_MEASURED
S4R_GRAPH_AWARE_INCREMENT_REAL_BUT_NOT_LIGAND_SPECIFIC
POSE_FREE_LIGAND_REPRESENTATION_REPAIR_ROUTE_CLOSED
S5D_LIGAND_STEERING_PRESENT_BUT_BIOLOGICALLY_MISDIRECTED
CONDITIONAL_ESTIMAND_REPAIR_ROUTE_CLOSED
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

## S5D — the follow-up diagnostic, and its falsification

S5D trained nothing and reused the frozen S4R checkpoints. It registered one
mechanism for R3 — that the estimator collapses every ligand difference onto
about one residue direction per protein — and **that mechanism is falsified**:

```text
rho_dg  (data-side bound)   0.4550      true-vs-foreign field cosine  0.4487
rho_graph (candidate)       0.4793      excess over data              0.0138
rho_base  (baseline41)      0.5758      rule needed >= 0.80 and +0.10
```

The estimator does steer on the ligand: a foreign pair moves the residue field
by a large angle. The graph statistic also produces *more* diverse fields than
the mean-pooled baseline, exactly as the audit predicted.

D2 then measured the symmetric-difference conditional estimand, which cancels
pocket membership exactly and non-parametrically. On 40,157 eligible pairs
across 107 components:

```text
candidate 0.655030   foreign 0.655470   chance 0.643744
baseline41 0.638830  permuted 0.628586
E1 candidate - chance    +0.011285 [LCB -0.007749]  needs +0.05  FAIL
E2 candidate - foreign   -0.000440 [LCB -0.021814]  needs +0.03  FAIL
E3 candidate - permuted  +0.026444 [LCB -0.002977]  needs +0.03  FAIL
```

Terminal verdict: `LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED`. The S4R verdict is
unchanged.

## Next research decision

No experiment is authorized, and no repair of this estimand is eligible. S4R's
stopping rule closed the representation route; S5D closed the conditional
estimand route and forbids a fourth estimand variant on heldout-A, which has
now been consumed three times.

The question has changed shape. Ligand information is not lost upstream and is
not diluted by the metric — it reaches the residue field at roughly half the
field direction and points somewhere biologically wrong. The natural reading is
that the missing ingredient is **correspondence**, which ligand substructure
sits against which residue, and that a pose-free sequence-plus-2D estimand has
no channel to supply it. Testing that is a separately governed information
stage about geometry, with its own preregistration. It is not authorized here.

## Frozen

- heldout-B and R6 amplitude/B5 integration;
- ChEMBL/BindingDB affinity, DAVIS, KIBA and recipient labels;
- larger PLM, second protein encoder, attention stack, parallel branch,
  geometry/pose branch, typed interaction branch, KG, PU loss or affinity head;
- larger ligand vocabulary or radius as a rescue of S4R;
- a fourth estimand variant on heldout-A;
- few-shot sectioning and biological `z` admission;
- CSMO, Band, mesh and `A(F,z)=K(B(z)F(z))`.

## Read first

1. `report/s7_l2b_r0r/PHASE2B_S5D_EVIDENCE_CONSOLIDATION.md`
2. `report/s7_l2b_r0r/PHASE2B_S5D_GATE.json`
3. `report/s7_l2b_r0r/PHASE2B_S4R_EVIDENCE_CONSOLIDATION.md`
4. `report/s7_l2b_r0r/PHASE2B_S4R_REPRESENTATION_AUDIT.md`
5. `report/CURRENT_RESEARCH_STATUS.md`
6. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
