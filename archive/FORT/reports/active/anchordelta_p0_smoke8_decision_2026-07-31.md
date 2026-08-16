# AnchorDelta P0 smoke-8 diagnostic

This run is deliberately non-decisional: eight fit targets, eight gate targets,
two pairs per fit target, one epoch, and one seed. It covered six independent
homology components.

| Arm | RMSE | Spearman | Pairwise |
| --- | ---: | ---: | ---: |
| AnchorDelta | 1.1431 | -0.0400 | 0.4609 |
| Calibration | 1.4726 | 0.0009 | 0.4993 |
| Permuted support labels | 1.1431 | -0.0400 | 0.4609 |
| Wrong-target support labels | 1.4749 | -0.0400 | 0.4609 |
| Wrong protein | 1.1430 | -0.0406 | 0.4603 |

The exact match between correct and permuted support labels is expected under
uniform aggregation: a pure column permutation cannot change a mean. It is
therefore only an invariance check, not a label-use test. Wrong-target labels do
degrade absolute RMSE, so the anchor coordinate is being used. However, correct
and wrong protein are effectively identical and ranking is below 0.5. The
protein-conditioned transfer mechanism is therefore not established. No graph
encoder, context transformer, uncertainty head, or task scheduler is justified
until the full component/seed gates pass.
