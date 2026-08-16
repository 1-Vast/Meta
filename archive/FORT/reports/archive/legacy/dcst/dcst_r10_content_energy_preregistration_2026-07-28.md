# DCST-R10 continuous content-energy preregistration

Date: 2026-07-28  
Status: frozen before implementation and training

## Hypothesis and innovation

R10 defines the **Privileged Measure × Content Energy (PMCE)** interface.
It loads and freezes the accepted R6 segment interaction representation. For
each target segment, the same R6 target projection produces a 32-dimensional
content token without adding the learned segment-position embedding:

```text
h_s = GELU(R6_target_projection(ESM_segment_s))
P_s,t = R6_privileged_interaction_measure(segment=s, type=t)
g(target, ligand) =
    sum_(s,t,d) P_s,t h_s,d Theta_d,t.
```

`Theta` is a `32 × 8` continuous-content × interaction-type energy matrix,
reset to the exact null and trained for 4,000 source steps. All upstream R6
parameters are frozen.

Unlike R6, the first axis of `Theta` is not an absolute sequence position.
Unlike R7-R9, it is not an arbitrary or clustered role. A spectral direction
is a continuous ESM content pattern × physical interaction-type pattern. The
pair-specific privileged measure decides where and for which ligand that
energy is expressed.

## Matched attribution control

PMCE-NoPriv loads the R6 no-privileged representation from the same checkpoint,
freezes it, resets the identical `32 × 8` matrix, and receives the same source
affinity steps. The only difference is upstream privileged structural
supervision.

The copied state-key audit, R6 checkpoint hash, exact-null preservation, and
trainable parameter count must be reported. No adapter, router, centroid, or
extra content MLP is allowed.

## Spectral certificate

For `Theta = U diag(s) V^T`, rank-direction `k` contributes

```text
s_k sum_(segment,type)
    P(segment,type)
    <h_segment, U_k>
    V_k(type).
```

The existing target and ligand destruction operators recompute or derange the
pair-specific measure. Two-direction bands, utility rule, scale, and held
source target split are unchanged.

## Source gate

Before another ChEMBL training run:

1. the frozen privileged R6 representation reproduces the full segment
   mechanism pass;
2. PMCE certifies at least one held-source continuous-content spectral band;
3. PMCE certifies strictly more bands than PMCE-NoPriv.

Failure stops R10.

## Conditional Stage 2

On source pass, the certified PMCE representation is frozen in the transfer
bridge. The downstream PMCE residual follows the already registered R6 arm
and control protocol using ChEMBL train/development only. MDE, paired
bootstrap, destruction removal, RMSE, negative-transfer, confirmation, and
sealed policies are unchanged.

The novelty claim is specifically structure-privileged interaction-measure
pretraining followed by continuous content-energy transfer with destructive
spectral certification, not generic cross-attention or encoder fine-tuning.

