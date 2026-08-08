# E-AFF-H0A Task-Local Radial Headroom Result

## Verdict

```text
TASK_LOCAL_RADIAL_HEADROOM_WITHOUT_PARTNER_SPECIFICITY
```

The fixed 288D radial basis contains substantial task-local held-out affinity
headroom, but the improvement is not sufficiently specific to the correct
protein under the frozen partner margin.

## Design

- 107 deep tasks from 107 closure components.
- One score-blind task per component and 40 score-blind ligand states per task.
- Per task: 20 fit ligands and 20 untouched test ligands.
- One fixed-alpha task-local Ridge direction fitted only on correct-protein
  residual differences.
- Correct, deranged and coupling-null test arms shared exactly the same task
  direction.
- No DAVIS or recipient label reads.

## Primary Test Results

| Arm/contrast | CI or delta | 95% component-bootstrap CI |
|---|---:|---:|
| Ligand OOF | 0.55404 | -- |
| Correct task-local | 0.64226 | -- |
| Deranged | 0.63362 | -- |
| Coupling null | 0.63817 | -- |
| Correct - ligand | +0.08822 | [+0.06761, +0.10998] |
| Correct - deranged | +0.00864 | [+0.00338, +0.01462] |
| Correct - null | +0.00408 | [+0.00061, +0.00747] |

The correct-minus-ligand headroom condition passed strongly. The partner
contrast is reproducibly positive but did not reach the preregistered `+0.03`
effect-size requirement. Therefore the result cannot be called target-specific
radial coefficient heterogeneity.

Ki and Kd both showed positive task-local headroom (`+0.0840` and `+0.1139`).
The Kd coupling-null contrast was negative, while only 15 Kd components were
available, so no endpoint-general coupling claim is supported.

## Interpretation

The shared-direction failure was not caused by a complete absence of affinity
signal in the radial tensor. A task can use the tensor to improve ranking of
unseen ligands. However, nearly the same learned direction transfers to a
wrong protein, so most of this headroom can still be explained by ligand-side
chemistry encoded in the pair tensor, assay/series-local structure, or a weak
protein marginal rather than a strongly correct-partner-specific mechanism.

Under the preregistered decision rule H0-B cross-assay target transport is not
authorized by this result. The next repair must first isolate a stronger
partner-conditioned component; it must not be described as a successful RFSA
or biological section.
