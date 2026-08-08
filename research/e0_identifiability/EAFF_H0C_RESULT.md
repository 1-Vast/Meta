# E-AFF-H0C Fixed Radial Interaction Residual Result

## Verdict

```text
FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED
```

After matching the 20-shot ligand nuisance and removing both fixed-tensor
marginals, the remaining radial chemistry-distance residual did not improve
held-out affinity ranking and did not distinguish the correct protein.

## Design

- 54 new tasks from 54 closure components; all 107 H0A tasks excluded.
- 20 support and 20 test ligands per task, with zero within-task Murcko
  scaffold overlap (`2,160` rows total).
- Frozen closure-OOF global ligand prior plus a support-matched task-local
  nuisance on the existing frozen 128D pooled ligand state.
- The interaction head was fitted only on support cross-fitted nuisance
  residual differences using the fixed algebraic tensor
  `psi=(phi-phi_null)/total`.
- Correct and deranged arms used exactly the same fitted direction. Wrong
  proteins never entered training.
- No DAVIS or recipient label reads.

## Primary Results

| Arm/contrast | CI or delta | 95% component-bootstrap CI |
|---|---:|---:|
| Global ligand OOF | 0.55487 | -- |
| Support-matched Local-L | 0.59635 | -- |
| Local-L + Interaction-C | 0.59244 | -- |
| Local-L + Interaction-D | 0.59092 | -- |
| Correct - Local-L | -0.00391 | [-0.02040, +0.01191] |
| Correct - Deranged | +0.00152 | [-0.01024, +0.01424] |

The Local-L gain over the global prior was `+0.04147`, confirming that the
support labels contain transferable series-local SAR. The centered radial
interaction residual did not add to it. Its tiny correct-minus-deranged
contrast crosses zero and is far below the frozen `+0.03` margin.

Ki (`48` components) had correct-minus-Local-L `-0.00646` and
correct-minus-deranged `+0.00024`. Kd had only six components; its positive
descriptive deltas cannot support an endpoint claim.

## Interpretation

H0C closes the proposed shortcut-removal rescue for the current fixed radial
basis. H0A's large task-local headroom is explained primarily by ligand/series
adaptation, not by a recoverable chemistry-distance coupling tied to the
correct protein.

This does not prove that source affinity lacks a protein-by-ligand interaction,
and it does not prove that orientation is required. The next discriminating
question is whether real same-document target-by-ligand double differences are
strong and reproducible. Only if that interaction exists but the fixed radial
basis cannot recover it would a new biological basis be justified.

## Audit

The independent audit passed every check. It reproduced all task, component
and bootstrap metrics exactly; verified 20/20 scaffold-disjoint partitions,
zero H0A task overlap, one-to-one `<0.40` derangement, artifact hashes, and
centered marginal error `2.29e-17`.

The first launch session was orphaned by the command timeout and produced only
selection/derangement files; no features or metrics were generated there. That
runtime residue was removed during repository consolidation under `history.md`
F-87. The completed immutable result remains `artifacts/eaff_h0c_v1_run2/`.
