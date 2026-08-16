# DCST-R13 privileged-prompted CFRI preregistration

Date: 2026-07-28  
Status: frozen before implementation and ChEMBL training

## Hypothesis and architecture

R13 tests whether the R6 structural measure is useful as privileged context
for a strong zero-shot Stage-2 interaction learner rather than as a complete
representation.

The trainable direct branch is CFRI-compatible:

- 32 ESM target segments projected to width 64;
- a Morgan/descriptor ligand encoder and two ligand-to-target cross-attention
  blocks;
- a 64-dimensional direct pair feature.

The frozen R6 teacher supplies the R12 `32 × 8 -> 256` structural interaction
moment. Concatenate the direct pair feature and structural moment, then use

```text
LayerNorm(320) -> Linear(320,128) -> GELU -> Linear(128,1).
```

The final layer starts at exact zero. The direct branch and fusion head train
on ChEMBL; the Stage-1 teacher remains frozen.

This **Privileged-Prompted CFRI (PP-CFRI)** gives Stage 2 its own
dataset-specific interaction capacity while exposing the high-quality
Stage-1 measure as a pair-specific inductive bias. It does not transfer source
affinity energies or read target support labels.

## Matched controls

All arms share the exact direct-branch and fusion-head initialization,
episodes, seed, optimizer, base residual, and 4,000 steps:

- `PP-CFRI-Priv`: R6 privileged moment;
- `PP-CFRI-NoPriv`: R6 no-privileged moment;
- `PP-CFRI-Uniform`: uniform structural measure; the direct CFRI feature
  remains fully trainable and therefore defines the strong no-prompt control;
- B0.

Parameter counts outside the frozen teacher are identical. Uniform removes
only pair-specific Stage-1 information.

## Frozen gate

Using ChEMBL strict dual-cold development and the existing MDE `0.0586`:

1. PP-CFRI-Priv minus B0 reaches the MDE;
2. its grouped-bootstrap LCB95 is positive;
3. it beats both PP-CFRI-NoPriv and PP-CFRI-Uniform with positive LCB95;
4. RMSE is no more than 5% worse than B0;
5. target and ligand destruction each remove at least 70% of the primary
   positive effect.

R6 source admission, ChEMBL train/development-only access, confirmation, and
sealed policies remain unchanged. Failure stops before confirmation.

