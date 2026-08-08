# F6I preregistration: protein-independent task location

Date frozen: 2026-08-08, after the single residual failure of F6P and before
evaluating this orthogonality correction.  KCGS numeric outcomes remain unread.

F6P described `tau` as a nuisance/task-location coordinate but computed it from
`r-s(P,L)`, allowing a wrong protein curve to change its own intercept.  F6I
enforces the stated product decomposition:

\[
 \tau(S)=\Pi_{[-0.5,0.5]}\frac{\sum_i r_i}{k+\lambda_\tau},
 \qquad \hat r(P,L)=s(P,L)+\tau(S).
\]

`tau` is now independent of protein features.  Therefore correct, zero, and
nearest-protein controls share the same location when given the same support
labels, and protein necessity is tested only through the interaction curve.
Wrong-target support receives its own (wrong) location.  Label scrambling and
support-pair reordering remain exact invariants.

All source selection, data, support/query rules, seeds, bounds, controls, and
componentwise gate are otherwise identical to F6P.  This is the only change.
