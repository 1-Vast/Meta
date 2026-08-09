# Preregistration — P1R2B-PHASE2B-S2R

## Gauge-free bounded binary ordinal witness

Stage identifier: `P1R2B-PHASE2B-S2R_GAUGE_FREE_DIRECT_W`

Written: 2026-08-10, before any S2R calibration or sealed metric. Historical
S0R and S1R remain failed development evidence. A post-S1R audit found that
S1R used float32 rather than the parent-contract float64 projection; S2R fixes
that contract defect explicitly and does not reinterpret S1R as a clean
loss-only experiment.

## 1. Question

Can the existing frozen residue/ligand states identify a transferable binary
ordinal residue statistic when both non-identifiable scales are removed?

The failed factor head has

```text
W = U^T V,  U -> A U, V -> A^{-T} V,
```

so Adam follows different paths for the same function. Separable logistic loss
also decreases along `cW, c>1` without changing AP. S2R removes both freedoms.

## 2. Frozen inputs and firewall

Reuse byte-for-byte the S0R metadata-only records, pair panels and 210-update
stream. No MONN residue edge, affinity value, ChEMBL, BindingDB, DAVIS, KIBA,
recipient or metaval value may be read.

All biology remains frozen: ESM2 residue states, 41-D ligand atom mean,
protein nuisance basis `Q_P`, exact sequence/scaffold closure and component
aggregation.

## 3. Only candidate

One matrix, no bias:

```text
W in R^(1280 x 41), 52,480 trainable parameters
d_raw = (I-Q_P Q_P^T) H_P W [g(La)-g(Lb)]
d_dir = d_raw / sqrt(mean_r(d_raw^2) + 1e-12)
```

The projection is computed in float64 exactly as in the parent Phase 2B
contract. `d_dir` is used only by the loss; raw and normalized scores have the
same ranking. After every optimizer update, project `W <- W / ||W||F`.

This is not an added parallel module. It replaces the gauge-dependent
factorization for this witness.

## 4. Objective and optimizer

Use the S1R all-residue bidirectional pairwise logistic loss on `d_dir`:

```text
L = 0.5 * [
 mean_{g in G,j notin G} softplus(-(d_g-d_j))
 + mean_{l in L,j notin L} softplus(-((-d_l)-(-d_j)))
]
```

Aggregation remains pair -> construct -> component -> batch. Frozen optimizer:

```text
Adam, lr=1e-3, no weight decay, gradient clip=5.0, 210 updates,
parameter seed=20260901, sampler seed=20260902.
```

No LR, budget, sampler, rank, margin, threshold or loss search is permitted.

## 5. Seed isolation

```text
development seed already burned: 20260921
fresh calibration seeds: 20260931, 20260932, 20260933
sealed verification seed: 20260997
bootstrap seed: 20260903
```

The sealed teacher may be instantiated once only after a persisted calibration
PASS artifact exists.

## 6. Numerical preconditions

Before calibration:

1. positive rescaling of any pair score changes AP by at most `1e-12`;
2. `d_dir(c d) = d_dir(d)` within `1e-8` for positive `c`;
3. antisymmetry and identical-ligand zero pass existing tolerances;
4. float64 projection matches the parent path within `1e-12`;
5. every checkpoint reloads, `||W||F` is within `1e-5` of 1 and predictions
   reproduce exactly.

Failure is `S2R_CONTRACT_INVALID`.

## 7. Calibration and sealed Gate

For every calibration seed, train once from seed 20260901. Complete-train and
complete-held-out component-macro `AP_bidir` must both be at least `0.50`.
All three seeds must pass. No mean-over-seeds rescue is allowed.

Only after calibration PASS, run sealed seed 20260997 once. It must also have
train and held-out `AP_bidir >=0.50`.

Report the singular spectrum of `W` and the AP of its rank-8 SVD truncation as
a non-gating compression diagnostic. Do not select rank from sealed results.

## 8. Terminal verdicts

Exactly one:

```text
S2R_CONTRACT_INVALID
BINARY_ORDINAL_TRAIN_FIT_NOT_IDENTIFIED
FINITE_DESIGN_GENERALIZATION_NOT_IDENTIFIED
GAUGE_FREE_BINARY_ORDINAL_VERIFICATION_FAILED
BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED
```

The earliest failed boundary controls the verdict. Only
`BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED` authorizes a separately frozen,
single-run real structural Phase 2B experiment. It does not authorize that run
automatically.

## 9. Mathematical and biological boundary

S2R identifies at most a bounded binary residue-ranking statistic. It does not
identify physical interaction energy, affinity direction, selectivity,
few-shot adaptation or a biological `z`. It does not modify

```text
A(F,z) = K(B(z)F(z)).
```

Any future biological admission still requires correct-ligand, wrong-ligand,
wrong-protein, independent structural and affinity-incremental Gates.
