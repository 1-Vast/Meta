# Shift and Generalization Audit

## Mandatory assumptions

Historical tasks influence the current decision only under a declared bridge.
At least one of the following is required:

- IID sampling from a common population law;
- a specified joint exchangeable model or ambiguity class including the current
  task, not merely symmetry among historical tasks;
- a declared transport class linking historical and current laws, such as a
  density-ratio or total-variation constraint.

For conditioning on noisy current observations, a stochastic likelihood is
required for a single posterior. With only the frozen bounded-noise support, the
honest result is a posterior ambiguity class over all compatible likelihoods.

Historical decision-functional estimation requires coverage or censoring rules
for each historical member. More historical members reduce sampling error, not
systematic censoring width.

If tasks carry contexts or the future query/context differs from historical
ones, the bridge must be conditional: conditional exchangeability/stationarity,
support overlap or positivity, and a declared conditional transport class. Phase
7's unconditional law on functions does not silently provide this.

## Robust decision

Given a nonempty weakly compact ambiguity class of current laws, a compact action
space, and suitable bounded lower-semicontinuous loss, robust expected-risk
minimization is well-defined. Weakening assumptions enlarges the ambiguity class.
With every law on the identified state set admitted, the criterion reduces to
set-robust minimization.

## Failure without a bridge

If the current population law may be any law supported on members of both
rankings, historical frequencies impose no restriction. An adversary can reverse
the historical majority. History then has zero robust value and a committal
historical-majority rule can be strictly worse than randomization or abstention.

## Audit conclusion

Phase 7 correctly identifies exchangeability/transport and likelihood as
declarations. It must not call exchangeability alone a unique predictive law,
and its generic learning claim needs an explicit population-law complexity
condition outside scalar/ranking regimes.
