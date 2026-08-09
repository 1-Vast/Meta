# S3R preregistration amendment 01 - pre-execution correction

Stage: `P1R2B-PHASE2B-S3R_REAL_STRUCTURAL_DIRECT_W`

Written and committed before S3R implementation and before any S3R score.
This amendment is binding wherever it is more specific than
`PREREG_PHASE2B_S3R_REAL_STRUCTURAL_DIRECT_W.md`. The parent registration is
retained byte-identical and marked
`SUPERSEDED_BEFORE_EXECUTION_ABSOLUTE_SCALE_AND_ORIGIN_NOT_IDENTIFIED`.

## A1. Train/held-out label firewall

Data preparation must materialize separate hashed train and held-out residue-mask
views. The candidate and R5 training processes may open only the train label
view. The candidate checkpoint and training trace must be frozen before a
separate scoring process opens held-out A labels. Held-out A is already a
development panel, not an independent confirmation, but implementation leakage
is still prohibited.

## A2. R4 isolates residue context

For R4, use the correct protein's frozen nuisance basis `Q_P`. Shuffle only the
contextual ESM2 residue states within amino-acid type and within the same exact
sequence, then compute:

```text
d_R4 = (I - Q_P Q_P^T) H_shuffled W [g(La)-g(Lb)]
```

Do not recompute `b^P` or `Q_P` from shuffled states. This keeps the control
attributable to contextual residue assignment rather than a changed projector.

## A3. Real stream and common masks

The 210-update stream is newly generated from the real eligible train pairs;
the S0R/S2R metadata-only stream is forbidden. Candidate, same-seed repeat and
R5 use one identical serialized stream and semantic hash.

Before every paired Gate, assert exact equality of pair IDs, construct IDs and
component IDs between arms. Intersections or arm-specific row dropping are
forbidden. A nonfinite prediction is fail-closed. A finite pair with
`g(La)=g(Lb)` or near-zero raw RMS is retained as a tied/chance contribution and
counted; it is not removed after labels are known.

The aggregation invariant is duplication of the complete observation panel.
No invariance is claimed for duplicating only one pair inside a construct.

## A4. R6 is removed from S3R adjudication

S3R does not compute or inspect the old R6 held-out integration metric. The
pairwise objective identifies only differences:

```text
Pi_P H_P W [g(La)-g(Lb)]
```

Unit norm selects one parameter representative but does not identify an output
amplitude relative to B5 logits. Pairwise labels also do not identify the
absolute ligand-feature origin or W directions that annihilate all observed
training ligand differences. Those directions cancel from R1-R5 but can change
a single-ligand `delta(P,L)`. Therefore `G_B5 + delta_raw` with coefficient one
is not an identified estimand.

If R1-R5 pass, the only authorized next action is a separately frozen
`S3R-I0_POSITIVE_SCALAR_B5_INTEGRATION` stage. I0 must use training ligands to
freeze a ligand-feature origin and the training difference-span projector,
freeze S3R W, and estimate one `alpha >= 0` from training components only with a
predeclared one-dimensional convex objective and positive ridge. Only after a
positive training alpha is frozen may held-out R6 be opened once. I0 is not
executed or scored in S3R.

Removing R6 is not a lowered Gate: its prerequisites are absent in the ordinal
contract, so it is deferred rather than passed.

## A5. Exact claim boundary

An R1-R5 development PASS may be described only as an annotation-defined,
ligand-conditioned residue ranking statistic on MONN development components.
It is not a low-rank claim: the primary learner is full W. The rank-8 SVD is
secondary and cannot rescue the candidate. It is not exact residue-atom
coupling, physical interaction energy, affinity, selectivity, few-shot
adaptation or biological `z` admission.

## A6. Revised terminal verdicts

Exactly one verdict is written in earliest-failure order:

```text
S3R_CONTRACT_OR_LABEL_FIREWALL_FAIL_CLOSED
S3R_DIRECT_W_TRAINING_OR_PARTICIPATION_FAILED
REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
NO_INCREMENT_BEYOND_FROZEN_B5_DIFFERENTIAL
REAL_RESIDUE_STATISTIC_SHORTCUT_DEPENDENT
STRUCTURAL_LIGAND_CONDITIONED_RESIDUE_DIRECTION_IDENTIFIED_IN_DEVELOPMENT
```

R1 failure selects the third verdict; R1 PASS and R2 failure selects the fourth;
R1-R2 PASS and any R3-R5 or module-participation failure selects the fifth; only
R1-R5 plus participation PASS selects the final verdict. The final verdict
authorizes I0 preregistration and one independent structural confirmation
design, but neither is run automatically and no affinity or production boundary
is opened.
