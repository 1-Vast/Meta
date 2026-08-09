# Preregistration — P1R2B-PHASE2B-S1R

## All-residue bidirectional ranking objective repair

Stage identifier: `P1R2B-PHASE2B-S1R_PAIRWISE_RANK_OBJECTIVE_REPAIR`

Written: 2026-08-10. Base evidence: S0R terminal verdict
`SURROGATE_AP_MISALIGNMENT_FULL_PANEL`. This document is frozen before any
S1R teacher, calibration or verification metric is computed.

## 1. Authorization and single change

S0R showed on all 44,746 held-out pairs / 81 closure components that the
registered group-balanced BCE decreases while bidirectional AP falls, even
from the balanced train-only ray optimum. S1R therefore changes exactly one
object: the surrogate loss. It does not change the head, rank, sampler,
optimizer, learning rate, weight decay, update budget, projection, panel,
threshold or model inputs.

The old BCE and its failed S0/S0R results remain historical evidence. They are
not overwritten.

## 2. Frozen data and firewall

Use the S0R metadata-only artifacts byte-for-byte:

```text
metadata_only_records.jsonl
  dcc2a0cd0cd958640f6548ed7cd3e5076e6c0ae9f6e455e35ad3e82f534dc856
train_pairs.jsonl
  847a0770856e7e26903e56bc40d7249bb4bc21082392aa7edf1e3166014c1195
heldoutA_pairs.jsonl
  29925ed5139b1d054aa3b2a26d5a0b281336e717fc801dcfae553b5e6cc340ae
```

The frozen 210-update stream from S0R is reused by semantic hash. No real
MONN residue-edge label, affinity value, ChEMBL, BindingDB, DAVIS, KIBA,
recipient or metaval value may be opened.

## 3. Fixed model and optimizer

The candidate remains the single 10,568-parameter rank-8 bilinear residue
residual head:

```text
delta(P, La, Lb) = (I - Q_P Q_P^T) H_P U^T V [g(La)-g(Lb)].
```

Frozen optimizer contract:

```text
AdamW; lr=1e-3; weight_decay=1e-4; gradient clip=5.0;
210 updates; identical hierarchical stream; parameter seed=20260901.
```

## 4. The only repaired object

For one synthetic ligand pair, let `G` be the teacher top-8 gain residues,
`L` the bottom-8 loss residues and `d_r` the predicted differential score.
Define

```text
L_gain = mean_{g in G, j notin G} softplus(-(d_g - d_j))
L_loss = mean_{l in L, j notin L} softplus(-((-d_l) - (-d_j)))
L_rank = 0.5 * (L_gain + L_loss).
```

All valid residues participate. There is no margin, temperature, class weight,
hard-negative mining, subsampling or point BCE. This is the standard smooth
pairwise ordering surrogate for the same gain and loss rankings measured by
`AP_bidir`; it removes the incompatible `unchanged target = 0.5` constraint.

Loss aggregation remains pair mean -> construct mean -> component mean ->
batch mean. Ordered `(a,b)` and `(b,a)` are not duplicated.

## 5. Seed isolation

```text
burned historical seed: 20260905 (never used in S1R)
calibration teacher seeds: 20260921, 20260922, 20260923
sealed verification teacher seed: 20260998
sampler seed: 20260902
bootstrap seed: 20260903
```

The sealed teacher is not instantiated until every calibration rule is
evaluated and a machine artifact records calibration PASS.

## 6. Calibration checks

For every calibration seed, run two trajectories on the complete panels.

### A. Teacher-start alignment certificate

Initialise `Head=(U*,V*)`, optimize the repaired loss for 100 updates, and
evaluate all 44,746 held-out pairs / 81 components at updates 0 and 100.
PASS iff:

```text
AP_100 >= AP_0 - 0.05
and UCB95(AP_100 - AP_0) >= -0.05.
```

The repaired loss must decrease on the complete train panel. Any seed failure
is `PAIRWISE_SURROGATE_STILL_MISALIGNED`; stop before scoring a student.

### B. Random-start recoverability

From parameter seed 20260901, train 210 updates. PASS iff every calibration
seed has complete-held-out component-macro `AP_bidir >= 0.50`. Report chance
AP and oracle-normalized recovery; do not lower 0.50.

If alignment passes but any student fails, terminal verdict is
`PAIRWISE_OBJECTIVE_ALIGNED_BUDGET_NOT_IDENTIFIED`. It authorizes only a
separate one-axis update-budget preregistration.

## 7. Sealed verification

Only after all three calibration seeds pass A and B, instantiate seed
20260998 once and run the identical teacher-start certificate and random-start
210-update student. No setting may change.

PASS requires the alignment rule and student `AP_bidir >= 0.50`. A sealed
failure is `PAIRWISE_OBJECTIVE_VERIFICATION_FAILED`; no retry or seed
replacement is allowed.

## 8. Artifact contract

Before scoring, freeze and hash the reused stream. For every evaluated
trajectory save checkpoint, optimizer state, complete held-out residue-score
table and pair/construct/component metrics. Reload each artifact and reproduce
AP/loss. Bootstrap units are closure components, never pair rows.

The program writes terminal and all downstream `NOT_RUN` artifacts itself.

## 9. Mutually exclusive terminal verdicts

```text
S1R_CONTRACT_INVALID
PAIRWISE_SURROGATE_STILL_MISALIGNED
PAIRWISE_OBJECTIVE_ALIGNED_BUDGET_NOT_IDENTIFIED
PAIRWISE_OBJECTIVE_VERIFICATION_FAILED
SYNTHETIC_IDENTIFIABILITY_REPAIRED
```

Only `SYNTHETIC_IDENTIFIABILITY_REPAIRED` authorizes a separately frozen,
single-run real structural Phase 2B contract. It does not itself authorize
real labels, affinity, structural confirmation, few-shot adaptation or z.

## 10. Frozen mathematical boundary

S1R is an upstream structural-control repair. It neither verifies nor modifies

```text
A(F, z) = K(B(z)F(z)).
```

No raw residue score or pair map enters `z`. Biology, affinity direction and
support-identifiable sectioning remain unclaimed.
