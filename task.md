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
UNTOUCHED_CORRESPONDENCE_CORPUS_IDENTIFIABLE
EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER
WITHIN_SLOT_DECONVOLUTION_SATURATED_BY_ADDITIVE_MARGINALS
X1A_AMENDED_AUDIT_HISTORICALLY_PASSED
X1A_ICC_PRECONDITION_NOT_ESTABLISHED
X1A_R_DIRECT_DD_RECTANGLE_MANIFEST_MATERIALIZED
CROSSED_INTERACTION_EXISTENCE_NOT_YET_TESTED
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

## C0/C1 — the correspondence question, answered on an untouched corpus

S5D pointed at **correspondence** as the missing ingredient. C0/C1 tested it
directly, audit-only, on 1,862 systems built from raw RCSB mmCIF coordinates
that no MetaSieve stage had ever touched — 24,874 exposed PDB ids excluded.

C0 passed every admissibility Gate: 496 independent inference components,
largest fraction `0.0811`, minimum detectable effect `0.00453` against a `0.05`
requirement. The union closure produced only 89 components and blew the giant-
component cap, so the registered DataSAIL-style fallback was used, exactly as
the instruction anticipated.

C1 then failed:

```text
within-slot AP, empirical              0.985611
within-slot AP, fixed-degree rewire    0.953959
C1a empirical - rewire  +0.031652 [LCB +0.029690]  needs +0.05  FAIL
```

Terminal verdict: `EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER`.

The decisive number is not the Gate but the ceiling. The maximum possible gain
above the fixed-degree null is `1 - 0.953959 = 0.046041`, so the `+0.05` margin
is unreachable in principle. The value `1 - 0.985611 = 0.014389` is only the
empirical residual to perfect AP. At the
frozen `6.0 A` threshold a slot holds about three sequence-adjacent — hence
spatially adjacent — residues, so there is very little left to deconvolve.

The C2 router was **not preregistered and not trained**, as the stopping rule
requires.

## X1A — crossed affinity dependence precondition

X1A tested whether the ChEMBL37 crossed design permits the estimand
`DD = y(P1,La) - y(P1,Lb) - y(P2,La) + y(P2,Lb)` to be tested at all. It
trained nothing and read ChEMBL37 only after its preregistration was committed.

```text
G1 Ki UCB95(rho) 3.88e-07 < 0.0915   G3 capped share Ki 0.0387 Kd 0.2066 <= 0.25
G2 Kd UCB95(rho) 1.33e-03 < 0.0164   G4 effective n Ki 827.0 Kd 604.3 >= 245
X1_ICC_PRECONDITION_PASSED
```

Two of this stage's own instruments were defective and were corrected before
the verdict: the registered ICC estimator was degenerate (a per-panel intercept
forces `var(cluster)` to zero for any data; its void first run was a pass), and
G3/G4 used measurement counts instead of the registered X0-B DD unit.

Three caveats travel with the historical pass. `rho` is biased toward zero by additive
over-parameterization and that bias favours passing; `var(panel)` truncated to
zero so the decomposition is only partly identified; and replicate noise is
over 99.9% of adjusted variance, giving detectable interaction RMS of `0.309`
(Ki) and `0.930` (Kd) log units at the frozen 0.5 ratio. Subsequent independent
review found that the global fit absorbs cluster-exclusive targets and that its
signed residual ICC does not identify dependence of `q=DD^2-v_noise`; the final
artifact also used 2,000 rather than 10,000 registered bootstrap draws.

## Next research decision

**X1A-R direct-DD dependence, and nothing else.** The label-blind rectangle
manifest now reproduces Ki `11,168/36` and Kd `1,041/12`, with frozen caps
selecting 827 and 605 units. X1A-R must estimate dependence of the same statistic
that X1B would test, `q=DD^2-v_noise`, without target/ligand nuisance fitting.
Current verdict is `X1A_ICC_PRECONDITION_NOT_ESTABLISHED`; X1B is not authorized.

X2 and every trainable component stay unauthorized until X1A-R and X1B pass
under separate preregistrations. Three routes remain closed by preregistered Gates:
representation (S4R), estimand (S5D) and geometry-gated correspondence (C1).

## Frozen

- heldout-B and R6 amplitude/B5 integration;
- additional ChEMBL fields, BindingDB affinity, DAVIS, KIBA and recipient
  labels; X1A's scoped 63,859 pChEMBL reads are recorded, with zero training reads;
- larger PLM, second protein encoder, attention stack, parallel branch,
  geometry/pose branch, typed interaction branch, KG, PU loss or affinity head;
- larger ligand vocabulary or radius as a rescue of S4R;
- a fourth estimand variant on heldout-A;
- the C2 geometry-gated correspondence router at the `6.0 A` / 128-slot contract;
- widening the correspondence corpus or relaxing its threshold to chase C1a;
- few-shot sectioning and biological `z` admission;
- CSMO, Band, mesh and `A(F,z)=K(B(z)F(z))`.

## Read first

1. `report/crossed_interaction/X1A_REPAIR_AND_CURRENT_BOUNDARY.md`
2. `report/correspondence_router/C0_C1_EVIDENCE_CONSOLIDATION.md`
3. `report/correspondence_router/C1_INFORMATION_AUDIT.json`
4. `report/s7_l2b_r0r/PHASE2B_S5D_EVIDENCE_CONSOLIDATION.md`
5. `report/s7_l2b_r0r/PHASE2B_S4R_EVIDENCE_CONSOLIDATION.md`
6. `report/CURRENT_RESEARCH_STATUS.md`
7. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
