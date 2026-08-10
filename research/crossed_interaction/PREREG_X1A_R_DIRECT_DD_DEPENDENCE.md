# Preregistration — X1A-R direct-DD dependence repair

Stage: `E-AFF-X1A-R_DIRECT_DD_DEPENDENCE`

Status: **registered, not executed**. This supersedes X1A Amendment 01 for
authorization purposes without deleting or rewriting its historical result.

## Why repair is required

The amended X1A global target/ligand fit does not identify the dependence that
controls the planned interaction-variance test:

- every target is confined to one dependency cluster, so target fixed effects
  absorb the cluster indicator;
- singleton ligand effects exactly fit their cell means;
- its residual ICC is measurement-level, while X1B would test the rectangle
  statistic `q = DD^2 - v_noise`;
- the final artifact used 2,000 rather than the registered 10,000 bootstrap
  draws.

Therefore `X1_ICC_PRECONDITION_PASSED` remains a historical amended-audit
result but no longer authorizes X1B execution. Current status is
`X1A_ICC_PRECONDITION_NOT_ESTABLISHED`.

## Frozen unit and estimand

First materialize the exact deterministic X0-B cell-disjoint rectangles using
the original greedy packing algorithm and byte-verified X0 cells. It must
reproduce exactly:

```text
Ki  11,168 rectangles, 36 clusters, cap 32
Kd   1,041 rectangles, 12 clusters, cap 125
```

For every rectangle, use the four registered cell means directly:

```text
DD = y(P1,La) - y(P1,Lb) - y(P2,La) + y(P2,Lb)
q  = DD^2 - v_noise
```

No target, ligand or panel nuisance model is fitted: their additive effects
cancel algebraically in DD. Ki and Kd remain separate.

`v_noise` is the sum of the four cell-mean measurement variances. Report two
frozen strata:

1. exact-assay replicate-supported rectangles, where all four variances are
   directly estimable;
2. pooled-noise rectangles, using an endpoint-level variance estimated only
   from exact-assay replicate cells and fixed before q is inspected.

The primary dependence quantity is the intra-cluster correlation of q, not of
signed affinity residuals. Estimate it with a one-way random-intercept profile
likelihood and a one-sided 95% profile-likelihood upper bound. Kd's 12 clusters
must additionally receive a conservative small-G sensitivity interval. No
ordinary pair/rectangle bootstrap is permitted.

Use the same conservative upper bound for both the rho threshold and effective
sample-size calculation. Per-cluster cap membership is fixed label-blind by
the lexicographic rectangle id in the materialized manifest.

## Verdicts

Exactly one:

```text
X1A_R_CONTRACT_INVALID
X1A_R_DEPENDENCE_PRECONDITION_FAILED
X1A_R_DEPENDENCE_PRECONDITION_PASSED
```

X1B is conditionally preregistered before any direct-DD value is opened. Only
the last verdict authorizes **execution** of X1B for the passing endpoint(s).
X2 and every trainable model remain unauthorized.

## Boundaries

No BindingDB, DAVIS, KIBA, PDBbind or recipient values. No model training. No
change to `model/`, production `scripts/`, `theory/`, z, CSMO, Band, mesh or
`A(F,z)=K(B(z)F(z))`.
