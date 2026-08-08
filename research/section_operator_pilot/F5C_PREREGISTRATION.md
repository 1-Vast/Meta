# F5C preregistration: counterfactual Fisher support design

Date frozen: 2026-08-08, after F4O and before evaluating this support policy.
KCGS numeric outcomes remain unread.

## Failure addressed

Ordinary D-optimal design conditions coefficients *within* the supplied protein
model.  It does not distinguish that protein from a plausible alternative.
F5C instead chooses supports to maximize local identifiability against the
nearest non-self KLIFS pocket.

For every candidate ligand, form the label-free vector

\[
 q_L=[1,s_0(P,L),\Delta s_0(L),\Delta s_1(L),\Delta s_2(L)],
 \quad \Delta s_j=s_j(P,L)-s_j(P',L),
\]

where `P'` is the nearest non-self pocket and `(s0,s1,s2)` are the three atlas
coordinates.  Columns are scaled by candidate-set standard deviations.  One
representative per generic Murcko scaffold is retained, then five supports are
chosen greedily to maximize `log det(1e-3 I + sum q q^T)`.  No activity outcome
enters support selection.

The predictive law is the bounded orthogonal F4O posterior and still adapts
only `(m,a)`, so `d_adapt=2<=5`; the extra contrast coordinates design the
experiment but are not fitted task parameters.

## Frozen protocol and gate

Penalties are reselected on PKIS1 simultaneous scaffold-cold and kinase-group-
cold folds using this exact counterfactual design.  Evaluation uses one
deterministic support set per target and a target-cluster bootstrap.  Query
compounds sharing any support scaffold are excluded for every arm.

Raw PKIS2 prediction must beat support-free, location-only, uniform protein,
nearest protein, wrong-target support, and permuted support with positive lower
95% bounds; all point estimates must be positive on Anastassiadis2011.
Centred interaction must beat nearest protein with the same criteria.

Support adaptation is not required to improve the already biological
support-free interaction curve.  Instead it has a fixed non-inferiority margin
of `1e-4` MSE versus support-free interaction on both panels.  This endpoint is
appropriate because amplitude calibration targets raw affinity while the
protein-identification claim is the correct-versus-wrong interaction contrast.

Passing licenses the statistic and its counterfactual observation policy for
later law-valued integration; it does not license a scalar bypass.
