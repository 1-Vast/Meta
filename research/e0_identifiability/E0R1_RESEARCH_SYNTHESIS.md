# E0R1 Three-Agent Research Synthesis

Updated: 2026-08-07

Three independent read-only agents reviewed `E0R1_PROPOSAL.md`, both external
analyses, E0S/E0R0 evidence, implementation and governance. None modified
files, trained a model, read real affinity or accessed DAVIS.

## Consensus

1. E0R0's first confirmed defect is objective semantics, not typed-tensor
   expressivity and not yet AdamW alone. Training ranked residual-only scores
   `Xw` against total labels `y=b+r`, while evaluation ranked `b+Xw`.
2. Aligned total-score logistic ranking is descriptive but is not an exact
   teacher-identification objective: finite perfect margins still have nonzero
   logistic gradient. The corrected synthetic objective is point residual
   Huber plus all-within-task residual-difference Huber.
3. `StandardScaler + Ridge(alpha=10)` is not an exact witness. It changes the
   objective and regularization geometry. E0R1-C uses raw centered float64
   Moore-Penrose `X+ r`, with no alpha or model selection.
4. Rank deficiency means parameter nonuniqueness, not automatically prediction
   nonidentifiability. E0R1-B must measure holdout row-space transport and the
   teacher-specific orthogonal contribution.
5. E0R1 runs A/B/C. Corrected deterministic D is conditionally authorized only
   if C reconstructs train residuals and passes the unchanged holdout Gate.
6. The historical 8-task holdout remains the primary paired diagnostic. It is
   not evidence for architecture-wide generalization. A new derangement map is
   selected without scores and must verify local identity `<0.40` pair by pair.

## Independent Findings

The objective-math audit found a 6.09% residual-order versus total-order
conflict, teacher-point old-rank gradient around `2.82e-3`, and near-zero point
residual/difference gradients. It found numerical rank `225/240`, effective
rank about `37.2`, condition about `1.14e8`, and high holdout row-space coverage.

The code audit found that minibatch pair coverage, feature conditioning and CP
initialization are real secondary issues, but none explains Full-240 before the
objective mismatch is removed. E0R0 trace gradients are minibatch gradients and
cannot certify full-risk convergence. AdamW weight decay contributes only a
small shrink and is not equivalent to Ridge.

The governance audit requires A/B/C plus conditional D, an exact `<40%`
derangement map, unchanged historical CI thresholds as continuity checks, and
continued freeze of affinity, DAVIS, PLIP/T, production, CSMO/Band and P2-P4.

This synthesis authorizes only
`P1R2B-E0R1_OBJECTIVE_DESIGN_SOLVER_AUDIT` under the registered contract.

