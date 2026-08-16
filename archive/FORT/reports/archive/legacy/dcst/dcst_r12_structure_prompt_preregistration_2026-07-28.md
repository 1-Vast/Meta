# DCST-R12 structural interaction-moment prompt preregistration

Date: 2026-07-28  
Status: frozen before implementation and ChEMBL training

## Hypothesis and innovation

R12 transfers the accepted R6 Stage-1 representation but no source affinity
matrix. For a target–ligand pair, the frozen privileged teacher predicts the
32-segment × 8-type measure `P(segment,type)`. The same frozen R6 target
projection produces a position-free 32-dimensional content token `h_segment`.
Their interaction moment is

```text
M[d,type] = sum_segment P(segment,type) h_segment[d].
```

The `32 × 8` moment is flattened and passed to a downstream-only head:

```text
LayerNorm(256) -> Linear(256,128) -> GELU -> Linear(128,1).
```

The final linear layer starts at exact zero, so the complete interaction
residual has an exact regular null. Stage 1 is frozen; only the 33k-scale
observation head trains on ChEMBL.

This **Structural Interaction-Moment Prompt (SIMP)** preserves the
pair-specific privileged interaction measure without assuming that source and
downstream affinities share an energy matrix. It is not generic encoder
fine-tuning: the only downstream input is a sufficient-moment candidate
constructed from the Stage-1 structural measure and position-free protein
content.

## Frozen Stage-1 admission

The exact R6 checkpoint and source report are fixed inputs. R12 requires:

- R6 privileged segment mechanism pass;
- R6 privileged 2/4 versus no-privileged 0/4 source certificate;
- exact checkpoint and target-feature hashes already recorded.

R12 does not retrain or reselect the source teacher.

## Matched Stage-2 controls

All heads use the same initialization, train episodes, base residuals,
optimizer, seed 1729, and 4,000 steps:

- `SIMP-Priv`: frozen R6 privileged measure;
- `SIMP-NoPriv`: frozen R6 no-privileged measure;
- `SIMP-Uniform`: privileged teacher content tokens but a uniform
  segment×type measure, removing pair-specific structural information;
- `B0`: the unchanged ligand-only base.

SIMP-Uniform has no ligand-varying moment within a target and is the exact
information-null control, not a lower-capacity network.

## Development gate

On ChEMBL-37 strict dual-cold development:

1. `SIMP-Priv - B0` mean target-Spearman gain is at least the frozen MDE
   `0.0586`;
2. its 95% grouped-bootstrap lower bound is above zero;
3. it beats both SIMP-NoPriv and SIMP-Uniform with lower bounds above zero;
4. RMSE is no more than 5% worse than B0;
5. target and ligand destruction each remove at least 70% of the positive
   primary effect.

The target and ligand destruction operators, bootstrap units, confirmation,
and sealed policies are unchanged. Failure stops before confirmation.

## Claim boundary

A pass supports a two-stage, structure-privileged dual-cold prediction claim
for the registered train/development substrate. It does not establish sealed
generalization until separately authorized.

