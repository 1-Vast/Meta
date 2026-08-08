# E-AFF-P0 Fixed-Radial Source Feasibility Result

## Verdict

```text
SHARED_DIRECTION_NOT_OBSERVED_H0_DATA_SUPPORTED
```

This research-only pilot read release-pinned ChEMBL37 Ki/Kd source labels. It
did not read DAVIS or recipient labels. It does not establish that the fixed
radial basis lacks affinity information.

## Panel And Controls

- 245 tasks from 245 closure components, one score-blind task per component.
- 20 distinct ligand states per task; 4,900 evaluated observations.
- Five fixed D1 outer folds; task-balanced residual-difference training.
- A closure-OOF 128D pooled-ligand Ridge prior trained on all 147,761 unique
  task-ligand observations derived from 152,737 governed measurements.
- One correct-only 288D shared radial direction per outer fold.
- Correct, score-blind deranged and marginal-preserving coupling-null arms used
  the same frozen direction; neither control entered training.
- 245 one-to-one wrong proteins, zero reuse, maximum local identity 0.39394.

## Primary Component-Macro Results

| Arm/contrast | CI or delta | 95% component-bootstrap CI |
|---|---:|---:|
| Ligand OOF | 0.55225 | -- |
| Correct | 0.54209 | -- |
| Deranged | 0.54445 | -- |
| Coupling null | 0.54210 | -- |
| Correct - ligand | -0.01016 | [-0.02069, -0.00018] |
| Correct - deranged | -0.00236 | [-0.01301, 0.00805] |
| Correct - null | -0.00001 | [-0.00558, 0.00530] |

All preregistered feasibility conditions failed. Ki and Kd were directionally
consistent for correct-minus-ligand (`-0.01108` and `-0.00680`). No endpoint
showed the required positive increment.

## Interpretation

The frozen radial tensor is structurally recoverable by T-BASIS-R0, but one
population-shared coefficient tensor did not provide source affinity ranking
value beyond the ligand prior on this panel. The result rejects the tested
shared-direction hypothesis; it does not reject target/task-specific radial
directions.

The label-blind H0 census found 891 tasks with at least 40 ligands across 107
closure components and 374 target-endpoint groups with at least two tasks and
documents. A separately registered task-local headroom diagnostic is therefore
data-supported. Cross-assay target transport remains a distinct subsequent
question.

## Audit

The independent post-run audit reproduced all component-macro values and
bootstrap intervals, verified all core output hashes, verified exact feature
shapes and finite values, and confirmed that the coupling null preserved both
marginals within `1.7e-7` maximum error.
