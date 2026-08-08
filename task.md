# Current task

## Objective

Determine whether the frozen P1B atom-local, residue-local and contact/distance
geometry contains a correct-protein structural mechanism that can later support
transferable affinity calibration and a `k<=5` identifiable meta-section.

## Active stage

```text
P1R2B-S5_LOCAL_MECHANISM_OBSERVABILITY
```

Registered at
`research/ssl_b2_structural_observability/PREREG_S5_LOCAL_MECHANISM_OBSERVABILITY.md`.

### Required order

1. Build deterministic CCD-atom and canonical-residue mappings.
2. Separate single-chain pockets from oligomer/interface complexes.
3. Audit a chemistry-faithful deterministic structural pseudo-teacher.
4. Quantify exact-residue to 128-slot information retention.
5. Run the actual frozen-P1B observability ladder with ligand-only, random,
   deranged-protein and pair-shuffle controls.
6. Run a synthetic trainability control.
7. Only if steps 1–6 pass, train a small research-only pair-local GPU head.
8. Freeze passing structural channels before any source-affinity stage.

## Frozen surfaces

- `theory/FINAL_FROZEN_THEORY/`
- production CSMO, Band, mesh and `K`
- DAVIS, recipient and other protected labels
- real ChEMBL/BindingDB affinity training until separately authorized
- production `z` admission
- P2–P4

## Promotion rule

A structural S5 PASS does not prove affinity.  Production admission additionally
requires source closure-OOF `correct-ligand>=0.03` and
`correct-deranged>=0.03`, both with 95% LCB above zero, then a sealed transfer
Gate.  Only afterward may a compact named mechanism basis enter a support
row-space-constrained section and the unchanged probability-law operator.

## Current verdict

```text
GEOMETRY_IDENTIFIED
AGGREGATE_ESM_ECFP_PROBE_NOT_PROTEIN_SPECIFIC
PAIR_LOCAL_P1B_OBSERVABILITY_NOT_TESTED
POSE_FREE_CLASS_NOT_CLOSED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
```
