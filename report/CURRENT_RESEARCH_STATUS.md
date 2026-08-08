# MetaSieve current research status

Updated: 2026-08-08.

## Proven facts

1. The frozen mathematical operator and its software contracts pass regression.
2. P1A/P1B established a reproducible open holo corpus and correct-protein
   contact/distance geometry.
3. ChEMBL37 D0/D1 established an immutable Ki/Kd source corpus with zero
   homology/document leakage across folds.
4. Historical MIF, nonlinear MIF, interaction-residual and shared radial
   affinity readouts did not establish a transferable correct-protein affinity
   increment over ligand-only and wrong-protein controls.
5. XP1/XP2 found interaction in consumed kinase panels, but the `k<=5`,
   double-held-out and external-replication requirements were not met.
6. XP3/XP4 quantified the public-data boundary: low-noise panels have too few
   independent protein groups, while broader BindingDB panels have assay noise
   larger than the estimated interaction component.
7. S2 produced deterministic six-channel structural pseudo-labels.  S4 showed
   that a mean-pooled ESM + ECFP Ridge probe predicts aggregate H-bond and
   hydrophobic totals, but that signal is ligand-dominated.

## Interpretation correction

S4 is not an upper bound on all sequence+2D models.  It omitted the actual P1B
atom-local states, residue-local states and atom-by-slot contact/distance
tensors.  The current status is therefore:

```text
GEOMETRY_IDENTIFIED
PAIR_COMPATIBILITY_IDENTIFIED
AGGREGATE_ESM_ECFP_PROBE_NOT_PROTEIN_SPECIFIC
PAIR_LOCAL_P1B_OBSERVABILITY_NOT_TESTED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

## Active research stage

The only active research stage is
`P1R2B-S5_LOCAL_MECHANISM_OBSERVABILITY`, registered at
`research/ssl_b2_structural_observability/PREREG_S5_LOCAL_MECHANISM_OBSERVABILITY.md`.

It performs, in order:

1. CCD atom and canonical residue mapping;
2. single-chain versus oligomer-interface deployment audit;
3. chemistry-faithful pseudo-teacher and 128-slot information ceiling;
4. actual frozen-P1B local observability ladder;
5. synthetic trainability control;
6. conditional lightweight pair-local GPU distillation.

The already opened 1,118-complex S4 set is development evidence.  A new
score-blind RCSB confirmation block is required.

## Promotion boundary

A structural S5 PASS only authorizes a separately registered source-affinity
Gate.  Production admission still requires source closure-OOF
`correct-ligand>=0.03` and `correct-deranged>=0.03`, both with 95% LCB above
zero, followed by a sealed transfer Gate.  Only then may a compact named
mechanism basis enter a `k<=5` row-space-constrained support section and the
unchanged `K(B(z)F(z))` operator.

## Repository policy

Terminal-negative implementations and duplicate reports were deleted after
their metrics and verdicts were consolidated into `history.md`.  They remain
recoverable from Git commits `3281780`, `12a2765`, and `608decf`.  Production
directories contain only passed core/data/geometry workflows; unvalidated S5
work remains isolated under `research/`.
