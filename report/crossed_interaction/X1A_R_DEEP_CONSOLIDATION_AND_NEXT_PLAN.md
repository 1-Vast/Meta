# X1A-R deep consolidation and next plan

> **SUPERSEDED PLAN.** X1A-R has executed and failed for both endpoints. X1B
> was not run. Current evidence is in `X1A_R_FINAL_SYNTHESIS.md`; the planning
> record below is retained for chronology.

Updated: 2026-08-10.

## Executive state

```text
C1 exact 6 A binary within-slot route ........ CLOSED
X0 original component unit ................... UNDERDETERMINED
X0-B cell-disjoint rectangle design .......... CONDITIONAL, REPRODUCED
X1A amended residual ICC ..................... HISTORICAL PASS, NOT VALID FOR AUTHORIZATION
X1A-R direct-DD dependence ................... REGISTERED, NOT EXECUTED
X1B interaction existence .................... NOT AUTHORIZED
X2 trainable q_theta .......................... NOT AUTHORIZED
GPU training ................................. NOT AUTHORIZED
```

MetaSieve still has no identified affinity-directed protein-ligand statistic.
The next failure boundary is data/estimand dependence, not neural capacity.

## What was repaired

### Statistical authorization

The amended X1A result is no longer treated as evidence that X1B can run. Its
fixed effects absorb cluster-exclusive targets and singleton ligands, and its
signed residual ICC does not measure dependence of the intended squared
double-difference statistic. This is an estimand failure, not a request for a
different optimizer.

### Statistical unit

X0-B previously preserved only counts and an algorithm description. The exact
rectangle rows are now materialized, hashed and fixed before direct DD values
are computed. Counts reproduce exactly, including frozen cap selection.

### Evidence governance

- ChEMBL37 X1A pChEMBL reads: 63,859.
- ChEMBL affinity training reads: 0.
- BindingDB/DAVIS/KIBA/PDBbind/recipient reads: 0.
- X1B execution: not authorized.
- Historical X1A artifacts: retained, not overwritten.

### C1 numerical language

`0.014389` is empirical-to-perfect residual. The maximum possible improvement
over the fixed-degree null is `0.046041`; that is the correct mathematical
reason a `+0.05` Gate cannot pass.

## Sole next stage: X1A-R

The stage must operate on the fixed rectangles and direct contrast. It must not
fit target or ligand nuisance functions.

For scale consistency with the X0-B cell interaction ratio:

```text
D   = DD / 2
v_D = (v11 + v12 + v21 + v22) / 4
Z   = D^2 - v_D
```

Primary rectangles require that, within each target, its two ligand values have
a common exact assay id. A label-blind rule chooses among multiple common
assays. Rectangles without this property do not silently enter the primary
stratum. After this filter, the frozen cap, effective sample size and largest
component weight must be recomputed.

Dependence and inference use closure components, never rectangle rows. Kd has
only 12 components and requires a conservative small-G analysis. X1A-R must
finish before X1B interaction existence is preregistered.

## Conditional X2 model

Only an endpoint passing X1A-R and X1B may train. X2 is limited to one
low-dimensional trainable `q_theta(P,L)` with at most four channels and rank at
most eight over existing frozen local protein, ligand and geometry states. No
new PLM, encoder, attention stack, parallel 3D branch, typed head, KG or support
adapter is allowed.

The primary objective is crossed difference only. A raw point-affinity loss is
forbidden because it reopens ligand and assay shortcuts.

## GPU interpretation

X1A-R/X1B are statistical stages and need not use GPU. For X2, performance is
measured by unique cells/s, rectangles/s, edge tokens/s, end-to-end versus
compute-only throughput and wall-clock completion. GPU utilization percentage
is diagnostic only. A small model must not be enlarged to raise `nvidia-smi`.

Before X2, cache all frozen embeddings and geometry, deduplicate `(P,L)` cells
within batches, bucket by token count, use pinned/nonblocking transfers and
mixed precision where numerically equivalent, and freeze the best label-blind
batch configuration before scientific training.

## Stop tree

```text
rectangle/assay manifest mismatch
  -> contract failure; no interaction claim

assay-aligned effective design < 245 or component dominance > 0.25
  -> data underdetermined; no model training

X1A-R dependence UCB exceeds frozen threshold
  -> dependence precondition failed; no X1B

X1B noise-corrected interaction absent or below margin
  -> source/estimand unsupported; no X2

X2 synthetic or module-participation failure
  -> implementation failure; no biological conclusion

X2 does not beat ligand-only and wrong-protein controls
  -> shortcut or representation failure; no section/z

all development Gates pass
  -> freeze q_theta and register independent confirmation only
```

The frozen operator `A(F,z)=K(B(z)F(z))` remains unchanged throughout.
