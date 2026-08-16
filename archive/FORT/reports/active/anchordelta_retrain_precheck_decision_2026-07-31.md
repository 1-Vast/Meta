# AnchorDelta trainable-encoder retraining decision

This run retrained the interaction encoder and antisymmetric comparator jointly
from the existing protein-conditioned checkpoint. It used 64 fit targets, 30
gate targets, 23 independent homology components, 32 balanced pairs per target,
two epochs, and one seed.

| Arm | RMSE | Spearman | Pairwise |
| --- | ---: | ---: | ---: |
| Retrained AnchorDelta | 1.0775 | 0.2521 | 0.5886 |
| Calibration | 1.2402 | -0.0026 | 0.4872 |
| Wrong-target support labels | 1.5173 | 0.2521 | 0.5886 |
| Wrong protein | 1.0775 | 0.2522 | 0.5887 |

The retraining improves ranking over calibration, but the wrong-protein arm is
indistinguishable from the correct-protein arm. Component-bootstrap intervals
for the correct-minus-wrong-protein RMSE, Spearman, and pairwise gains all cross
zero. Wrong-target labels change absolute RMSE but not ranking, confirming that
the ranking improvement is not protein-specific.

Decision: `NO_GO_FOR_PROTEIN_CONDITIONED_ANCHORDELTA`. Do not report this as a
protein-conditioned improvement and do not add graph/context/uncertainty
modules. A future improvement requires a new validated protein-interaction
signal or an explicitly different data protocol, not more epochs or a larger
head.
