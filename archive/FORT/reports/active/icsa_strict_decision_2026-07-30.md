# ICSA-DTA Strict D0 Decision

Date: 2026-07-30. Environment: `D:\anaconda\envs\drug`, CUDA, RTX 4060.
Result: `STOP_A1_BELOW_CALIBRATION`.

## Registered Question

Does an identifiable, protein-free, rank-2 contrast posterior improve over an
affine calibration null for k=5 unseen-target adaptation? The posterior uses
leave-one-out calibration residuals and projects both residuals and features
onto the exact orthogonal complement of `[1, B0]` before the Bayesian update.

## Protocol

- 64 TRAIN episodes, one epoch, seed 1729.
- 30 strict pKi development episodes in 29 homology components.
- Support and query are disjoint in chemical component, document token, and
  assay token.
- Exact evidence gate is the primary arm; joint and soft arms are diagnostics.
- Protein conditioning is disabled. The complete model has 299,164 parameters,
  only 22 more than the frozen ligand model.

## Result

| Arm | RMSE | Spearman | Pairwise |
| --- | ---: | ---: | ---: |
| Calibration | 1.35076 | 0.04969 | 0.52177 |
| Ligand posterior | 1.40082 | 0.04535 | 0.51532 |
| Gradient baseline | 1.62505 | 0.01191 | 0.50561 |
| ICSA exact | 1.35170 | 0.05103 | 0.52129 |

Component-paired ICSA-minus-calibration gains:

| Metric | Mean gain | 95% bootstrap interval |
| --- | ---: | ---: |
| RMSE | -0.00102 | [-0.00896, 0.00491] |
| Spearman | +0.00139 | [-0.01178, 0.01973] |
| Pairwise accuracy | -0.00050 | [-0.00543, 0.00517] |

All admission lower bounds fail. The ranking inclusion probability remains near
its 0.05 prior (mean 0.0521), so the support contrasts provide little evidence
for the added rank-2 function.

## Decision

Do not run the protein-conditioned increment, add backbone capacity, or tune
epochs against this development roster. The strict result localizes the current
limit to the identifying information in k=5 support observations, not posterior
expressiveness. The next credible protein-specific route requires new
interaction-identifying supervision or a denser, provenance-closed factorial
panel. The implemented contrast posterior remains a valid, lightweight
diagnostic and exact nested null, but it is not admitted as a performance or
paper innovation claim.

Machine-readable result: `reports/active/icsa_strict.v1.json`.
