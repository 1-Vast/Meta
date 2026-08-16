# Contrast Objective Diagnosis

Date: 2026-07-30  
Evidence: `reports/active/adaptcontrast.v1.json`  
Scope: 64 TRAIN episodes and 64 development episodes; diagnostic only.

## Decision

Adding a unit-weight antisymmetric query contrast loss did not pass the
predefined 64/64 mechanism gate. It is stopped and is not the default training
objective.

| Metric | Joint posterior | Ligand null | Previous joint run |
| --- | ---: | ---: | ---: |
| RMSE | 1.188 | 1.190 | 1.167 |
| MAE | 0.973 | 0.977 | 0.958 |
| Spearman | 0.123 | 0.117 | 0.165 |
| Pairwise accuracy | 0.541 | 0.543 | 0.562 |
| NLL | 1.664 | 1.676 | 1.626 |

The component-bootstrap gain over the ligand null crossed zero for RMSE,
MAE, Spearman, and pairwise accuracy. Correct protein also did not beat the
deterministic wrong-protein control. Across 2,493 queries, their increments
relative to the ligand null remained correlated at 0.932.

The run used CUDA throughout the numerical path. Wall time was 54.0 seconds;
mean/peak utilization was 30.2/53%, mean/peak board power was 10.22/22.88 W,
and peak NVIDIA/Torch memory was 5,072/2,246 MiB.

## Consequence

Do not expand this objective to 480/133, add epochs, or tune its weight on the
development split. The next isolated mechanism is a TRAIN-only
sequence-conditioned pair prior with a correct-versus-wrong-protein gate,
followed by the same exact k=5 support posterior.
