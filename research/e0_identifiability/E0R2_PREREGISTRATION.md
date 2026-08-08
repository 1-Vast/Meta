# E0R2 Proposal And Numerical Closure Preregistration

Registered: 2026-08-07

## Scope

This synthetic-only research stage tests the first executable claim in the
directional-potential proposal: the corrected residual objective must be
numerically realizable before any new biological representation is introduced.

It does not read real ChEMBL affinity values, DAVIS, recipient labels, or typed
interaction labels. It does not run PLIP, construct a reference-state potential,
train a few-shot adapter, modify the frozen theory, or change `model/` and normal
`scripts/`.

## Frozen Inputs

- the E0R1 `synthetic_design.npz`;
- the original 32-train-task / 8-holdout-task synthetic selection;
- the corrected score-blind `<40%` derangement features materialized by E0R1;
- the 240-dimensional centered typed statistic;
- the residual target and residual-difference target;
- `rcond = 1e-10`, with no model or tolerance selection.

The existing eight-task holdout is a historical development diagnostic, not an
untouched estimate of architecture-wide generalization.

## Deterministic Solve

For train design `X`, residual `r`, and all same-task pair differences `D`, solve

```text
min_w mean((Xw-r)^2) + mean((DXw-Dr)^2)
```

by one float64 SVD least-squares solve of the normalized augmented system. At
residuals below the Huber transition this is the exact quadratic region of the
registered point-Huber plus difference-Huber objective.

## Gate

- train RMSE `<= 1e-6`;
- corrected objective `<= 1e-12`;
- corrected full-gradient L2 `<= 1e-8`;
- correct holdout CI `>= 0.80`;
- correct-minus-ligand CI `>= 0.10`;
- correct-minus-deranged CI `>= 0.10`;
- permutation error remains `<= 1e-6` by the frozen statistic contract.

PASS means only `SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED`. It does not mean
that affinity energetics, directionality, a reference-state potential, or a
production biological statistic has been identified.

