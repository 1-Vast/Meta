# AnchorDelta P0 single-seed precheck

This is a larger diagnostic run, not the final multi-seed gate: 64 fit targets,
30 gate targets, 23 independent homology components, 32 pairs per fit target,
and two epochs.

| Arm | RMSE | Spearman | Pairwise |
| --- | ---: | ---: | ---: |
| AnchorDelta | 1.0886 | 0.1467 | 0.5440 |
| Calibration | 1.2402 | -0.0026 | 0.4872 |
| Wrong-target support labels | 1.5300 | 0.1467 | 0.5440 |
| Wrong protein | 1.0888 | 0.1430 | 0.5419 |

The anchor output improves the absolute RMSE over calibration, and replacing
the support labels worsens RMSE. However, the wrong-protein arm is essentially
identical to the correct-protein arm, and the ranking metrics are unchanged by
wrong support labels. This indicates that the current gain is an anchor/offset
effect rather than evidence of protein-conditioned relative affinity.

Decision: `NO_GO_FOR_P1`; do not add graph/context/uncertainty modules. A formal
multi-seed P0 may still be run to close the hypothesis, but it must pass the
wrong-protein and similarity-bin gates before any architecture expansion.
