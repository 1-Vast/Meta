# Biological gauge and partner-information audit

Date frozen: 2026-08-11

## Scope

Run a source/meta-validation audit only. Do not train or modify an affinity
model, open main-v0 test, add an anchor, unfreeze the frontend, or migrate code
to production. The audit distinguishes:

```text
H1  frozen representation lacks partner-selective affinity information
H2  information exists but the section is locally gauge-equivalent
H3  information exists before the 288D aggregation but is lost by aggregation
```

Stages are sequential. A later stage may run only if its data-independence Gate
passes. No stage can produce a confirmation claim.

## A0: gauge audit

Use only the five cluster-balanced v0 development checkpoints and the sealed
meta-validation episodes. For each `(seed,target,draw)` form correct/wrong
support matrices `M_c,M_w` and held-out query matrices `Q_c,Q_w`.

Fit the local orthogonal map on support only:

```text
R* = argmin_R ||M_w R - M_c||_F, R^T R = I
```

Without refitting, report support and query Procrustes residuals, support/query
Frobenius norms, Gram traces/eigenvalues, normalized and unnormalized Gram
distortion, and effective ridge scale `lambda/(trace(G)/k)`.

For all four support/query combinations report

```text
H_XY = Q_Y M_X^T (M_X M_X^T + lambda I)^-1
```

and decompose full prediction into `mu_q`, support residual and `H_XY r_X`.
Pair recomputed quantities with the existing correct/correct, correct/wrong,
wrong/correct and wrong/wrong prediction rows by
`(seed,target,draw,cell_id)`.

Controls: exact orthogonal rotation, sign/permutation rotations, non-orthogonal
shear, scale `0.1/10`, rank-deficient support, and ligand-row shuffle. Exact
orthogonal controls must preserve `H` and correction to numerical tolerance;
scale and shear must not be declared gauge-equivalent. Aggregate queries to
episode, episodes to target, and targets to CD-HIT40 cluster before inference.

A0 supports only **local gauge-like equivalence** when support-fitted maps also
fit held-out queries and low `H_CC-H_WW` distortion tracks low full-prediction
distortion. A leave-one-cluster-out global map is reported separately; failure
of a global map forbids calling the effect a single parameter gauge.

## A1: measured selectivity dependency Gate

Use only source measured-partner groups. Before reading a response in any probe,
aggregate repeated targets within `(panel,ligand,CDHIT40-family)` by median.
One `(panel,ligand)` contributes one normalized group loss, regardless of its
number of targets or family pairs.

Construct label-blind dependency closure. Groups sharing any document, exact
ligand, Murcko scaffold, CD-HIT40 family, or verified protein local identity
`>=0.40` belong to one component. Report components, largest component share,
family/ligand/scaffold/document depth and component-unit MDE.

Probe training is authorized only with at least 18 eligible components,
component MDE `<=0.600`, and scoreable groups in every outer fold. Otherwise
stop A1 at `SELECTIVITY_PROBE_NOT_IDENTIFIABLE` and do not run A2.

If authorized, outer folds are complete components. Standardization, optional
dimension reduction and ridge selection occur inside outer training only.
Compare capacity-matched nuisance/additive, correct-pair and rewired-coupling
null probes. Use group-macro then component-macro loss. Primary fixed-sequence
contrasts require one-sided 95% component-bootstrap LCB above zero and at least
999 block-label refits with p<=0.05. A planted positive control must pass.

## A2: information localization

A2 is authorized only after an identifiable A1 probe. Under exactly the same
components, folds, capacities and nulls, compare frozen representation levels:

```text
protein length/composition
ESM pooled protein
288D calibrated T-BASIS
calibration-input interaction aggregate
atom-residue contact/distance-derived local summary
```

These are frozen representation probes, not frontend training. If A1 cannot
support component-level inference, A2 remains unrun rather than reported from
dependent pair rows.

