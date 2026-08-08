# E0R1 Objective, Design And Solver Audit

Decision: `NOT_RUN_NUMERICAL_PRECONDITION_FAILED`.

## Objective

Residual-order versus total-order conflict is
`6.0855%`. At the analytic
teacher, old rank full-gradient L2 is
`0.00282`, while the
residual-difference gradient is
`3.58e-10`.

## Design

The 640 x 240 train design has rank `225`,
effective rank `37.20` and identified
condition number `1.14e+08`. Mean
holdout row-space coverage is `0.999816`;
teacher-specific unseen relative L2 is
`0.01113`.

## Exact Witness

Primary Moore-Penrose train RMSE is `3.18e-08`. Holdout
correct/deranged CI is `0.99737/0.61447`,
with partner delta `+0.38289`. Gate:
`True`.

## Corrected Deterministic Solve

Authorized: `True`. Converged:
`False`.
Correct/deranged CI:
`0.98421/0.61118`.
Gate pass: `False`.

The deterministic solver's prediction metrics would satisfy the historical CI
checks (`0.98421/0.61118`, partner delta `+0.37303`), and its relative gradient
is `2.60e-9`. However, train RMSE is `7.92e-4`, above the preregistered `1e-6`
numerical precondition after 526 closure calls. Its Gate was therefore not
evaluated. The high CI cannot be substituted for numerical convergence.

The formal findings are:

```text
OBJECTIVE_SEMANTICS_DEFECT_CONFIRMED;
TRAIN_DESIGN_PREDICTION_IDENTIFIABLE;
EXACT_LINEAR_WITNESS_PASS;
CORRECTED_ITERATIVE_SOLVER_NOT_RUN_NUMERICAL_PRECONDITION_FAILED.
```

No real affinity, DAVIS, PLIP/T, production, CSMO/Band or downstream work was
executed or authorized.
