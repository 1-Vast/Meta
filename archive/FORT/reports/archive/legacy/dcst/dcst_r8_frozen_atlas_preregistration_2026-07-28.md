# DCST-R8 frozen protein-language atlas preregistration

Date: 2026-07-28  
Status: frozen before atlas construction, implementation, and training

## Hypothesis and innovation

R8 replaces R7's learned router with a **Frozen Protein-Language Atlas
(FPLA)**. Eight role centroids are fitted once from the L2-normalized frozen
ESM segment embeddings of unique, firewalled PLINDER source-train exact
targets only. No affinity, structural annotation, source-development target,
or downstream embedding participates in atlas fitting.

Atlas construction is fixed:

- `MiniBatchKMeans`, eight clusters;
- random seed 1729, `n_init=20`, batch size 2048, maximum 200 iterations,
  `reassignment_ratio=0`;
- each unique train target contributes exactly 32 segments once;
- output centroids are L2-normalized and stored with the source-feature hash,
  target count, segment count, cluster occupancy, and zero-label provenance.

At model time, each target segment is assigned one-hot to its nearest frozen
centroid by cosine similarity. There is no trainable router and no
temperature:

```text
role(segment) = argmax_r cosine(ESM_segment, centroid_r)
P(role,type) = sum_segment P(segment,type) 1[role(segment)=role]
g = <Theta_atlas, P(role,type)>.
```

The atlas is an immutable, content-defined cross-protein coordinate. Stage-1
structure and affinity gradients can change the predicted segment interaction
measure but cannot rotate, exchange, or erase role identity.

## Frozen comparison

R8 changes only R7's learned role assignment into fixed hard atlas assignment.
It retains:

- eight roles and an `8 × 8` exact-null affinity matrix;
- exact UniProt/SIFTS source data and cross-source firewall;
- 32 ESM segments and active Morgan environments;
- segment-level structural losses and all weights;
- seed, optimizers, 4,000-step source budgets;
- matched FPLA-NoPriv architecture and affinity exposure;
- existing two-direction spectral certificate and source gates.

The role count is inherited from the preregistered R7 and is not selected
using R7 outcomes.

## Source gate

Before another ChEMBL training run:

1. the original segment-level privileged mechanism probe passes;
2. FPLA has at least one held-source certified role×type spectral band;
3. FPLA certifies strictly more bands than FPLA-NoPriv.

Failure stops R8 before another downstream label load.

## Conditional Stage 2

On source pass, the same centroid artifact is embedded in every Stage-2 arm
and frozen transfer bridge. The previously verified ChEMBL 32-segment ESM
cache is routed by the source-train atlas without refitting.

All R6 Stage-2 models, paired target bootstrap, MDE, destruction, negative
transfer, RMSE, confirmation, and sealed policies remain unchanged. R8 must
pass the full development gate before confirmation is authorized.

