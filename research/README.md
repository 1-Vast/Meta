# Research boundary

## Current executable path

`s7_l2b_r0r/` contains the completed MONN provenance/mapping workflow, frozen
ESM2 B5 exact-residue localizer, Phase 2A attribution audit, synthetic S2R
direct-W witness and real structural S3R transfer.

```text
Phase 2A  LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
S2R       BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED
S3R       REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
```

S3R is the current stopping point. It used one bounded direct matrix over frozen
ESM2 residue states and the existing mean-pooled 41-D ligand atom features.
Training and deterministic replay passed, but the candidate failed R1-R5. The
result does not authorize another optimizer attempt or a larger parallel model.

## Eligible proposal, not authorized

The only evidence-aligned representation proposal is a single-axis ligand
information audit: replace the global ligand mean by one frozen graph-aware 2D
ligand statistic while preserving the protein states, direct-W estimator,
ordinal objective, sampler, closure split and R1-R5. It must be separately
preregistered and must include ligand-only, foreign-ligand, context and
permuted-label controls.

Do not simultaneously add attention, a new PLM, geometry, pose, typed channels,
affinity supervision or KG features. A failed single-axis audit would close that
specific pose-free repair route and justify either a separately governed 3D
information stage or stopping.

## Promotion rule

Research code may enter `model/` or `scripts/` only after its own frozen Gate,
independent structural confirmation, source-affinity increment over ligand-only
and wrong-protein controls, and sealed transfer. No current research statistic
is admitted to production `z`.

Historical failed and superseded stages are evidence in `history.md` and Git,
not current execution instructions.
