# Current task

## Current state

```text
PHASE2A_TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED
S2R_BINARY_ORDINAL_TRAINABILITY_PASS
S3R_REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
CURRENT_ESM2_PLUS_POOLED_41D_LIGAND_BASIS_NOT_ADMITTED
AFFINITY_DIRECTION_NOT_TESTED
K_SHOT_SECTION_NOT_TESTED
BIOLOGICAL_Z_NOT_ADMITTED
```

S3R completed the authorized real structural transfer in the `drug` environment.
The candidate scored `AP_bidir = 0.03588` versus chance `0.02547`; its gain
`+0.01041 [LCB +0.00692]` failed the registered `+0.05` margin. It also failed
to beat B5, foreign-ligand, context-corrupted and trained permuted-label controls
by their registered margins.

The model trained and replayed correctly. The current failure is therefore the
real structural identifiability of the **specific** frozen ESM2 plus mean-pooled
41-D ligand basis, not an optimizer failure and not proof that all sequence+2D
representations are insufficient.

## Next research decision

No experiment is automatically authorized. The highest-value eligible proposal
is a separately preregistered, single-axis representation audit that replaces
only the ligand mean with a frozen graph-aware 2D ligand statistic while keeping
the protein states, direct-W ordinal estimator, split, loss, stream and R1-R5
unchanged. It must not be executed until its teacher ceiling, information gain,
parameter budget and controls are frozen.

## Frozen

- heldout-B and R6 amplitude/B5 integration;
- ChEMBL/BindingDB affinity, DAVIS, KIBA and recipient labels;
- larger PLM, new protein encoder, attention stack, geometry/pose branch, typed
  interaction branch, KG, PU loss, affinity head or parallel module;
- few-shot sectioning and biological `z` admission;
- CSMO, Band, mesh and `A(F,z)=K(B(z)F(z))`.

## Read first

1. `report/s7_l2b_r0r/PHASE2B_S3R_EVIDENCE_CONSOLIDATION.md`
2. `report/s7_l2b_r0r/PHASE2B_S3R_GATE.json`
3. `report/CURRENT_RESEARCH_STATUS.md`
4. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
