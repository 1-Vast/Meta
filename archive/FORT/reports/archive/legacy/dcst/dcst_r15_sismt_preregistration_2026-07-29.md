# DCST-R15 SISMT preregistration

Date: 2026-07-29  
Status: frozen before implementation and execution

## Claim

`Support-Intersection Spectral Mechanism Transport` tests whether R6 contains
a small intersection of directions that are both source-mechanism certified
and variable in the ChEMBL-train covariate domain. It does not align the full
domains and does not retrain either encoder.

## Inputs and representation

- PLINDER: the 2,106 exact-target, firewalled R6 TRAIN rows;
- ChEMBL: TRAIN covariates only while constructing the projector;
- frozen R6-Priv and R6-NoPriv teachers;
- matched uniform, random-frozen, wrong-target, and wrong-ligand controls.

For every pair, use the `32 x 8` R6 mechanism probability, residualized
against the same fixed 16 ChEMBL-train ligand anchors:

```text
z(t,d) = P_R6(t,d) - mean_a P_R6(t,a).
```

No ChEMBL affinity column may be requested or loaded. Confirmation and sealed
features remain excluded.

## Common spectral support

Estimate source and target covariance on CUDA and solve

```text
Sigma_T v = lambda (Sigma_S + epsilon I) v
```

through a symmetric source-whitened eigendecomposition. The ridge is
`1e-4 * trace(Sigma_S) / 256`. A direction is support-compatible when
`lambda` is in `[0.25, 4.0]`.

R6 source `theta` is decomposed into its eight rank-one mechanism directions.
Only its four directions belonging to the two held-source active certificate
bands carry privileged credit. A common eigendirection must have squared
projection of at least `0.05` onto that active basis. NoPriv and random
teachers receive no mechanism credit unless their own preregistered source
certificate is active.

Twenty target-block bootstrap repetitions re-estimate the common subspace.
Every retained direction must have median squared projection of at least
`0.50` into the corresponding bootstrap support subspace. At most 16
directions may be retained.

## Partial transport selector

Partial entropic transport is computed only after projection. Equal-weight
protein, ligand, and retained-mechanism views form the cost. Protein and
ligand views are reduced on source plus ChEMBL TRAIN covariates only; transport
never edits or synthesizes a representation.

For transported-mass fractions `0.10, 0.20, 0.40, 0.60`, report:

- achieved real-to-real mass;
- source exact-target ESS;
- ChEMBL target ESS;
- target and ligand coverage;
- maximum and 99th-percentile source weights.

The transport implementation uses an augmented dummy-source/dummy-target
Sinkhorn problem; dummy-to-dummy transport is prohibited. A numerical result
is valid only when both marginal maximum error and transported-mass error are
below `1e-3`.

## Frozen admission gate

SISMT enters affinity fitting only if all hold:

1. Priv retains between 1 and 16 stable directions;
2. NoPriv, uniform, and matched random retain fewer directions;
3. wrong-target and wrong-ligand each remove at least 50% of the Priv retained
   mechanism score or transportable mass;
4. a transported mass of at least `0.20` has both source-target ESS and
   ChEMBL-target ESS at least `88`;
5. at that mass, at least 20% of ChEMBL targets and 20% of ChEMBL ligands
   receive at least half their uniform transported mass.

Failure returns `STOP_SISMT_LABEL_BLIND_SUPPORT_GATE`. It does not permit
tuning thresholds or reading affinity.

## Conditional Stage 2

If admitted, freeze B0 and the projector. The only new affinity parameter is
an exact-null linear coefficient:

```text
y_hat = B0(d) + theta' Q' z(t,d),  H0: theta'=0.
```

Priv, NoPriv, uniform, random, wrong-target, and wrong-ligand arms use the
same optimizer, initialization, and steps. The existing `0.0586` MDE,
positive grouped-bootstrap LCB95, RMSE safety, and 70% destruction-removal
requirements remain binding.

