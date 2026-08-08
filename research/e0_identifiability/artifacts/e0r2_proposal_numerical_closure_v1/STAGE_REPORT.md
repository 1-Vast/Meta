# E0R2 Proposal And Numerical Closure

Decision: `SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED`.

The corrected synthetic objective was solved once in float64 with a frozen SVD
augmented least-squares solve. Train RMSE is `3.19e-08`, full
gradient L2 is `6.15e-17`, and corrected objective is
`1.57e-15`.

On the historical eight-task development diagnostic, correct/deranged CI is
`0.99737/0.61447`; correct-minus-ligand
is `+0.51283` and correct-minus-deranged is
`+0.38289`.

This closes the synthetic objective/design/solver boundary only. The holdout is
not untouched and the inherited derangement control reuses wrong proteins. No
real affinity, directionality, reference-state potential, few-shot adapter, or
production biological statistic was tested. No code is authorized for promotion
to `model/` or normal `scripts/`.
