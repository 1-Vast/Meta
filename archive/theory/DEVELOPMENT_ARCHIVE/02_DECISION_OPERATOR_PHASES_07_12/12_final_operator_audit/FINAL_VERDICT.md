# Phase 12 Final Operator Learnability Audit

## Consolidated verdict

`OPERATOR_LEARNABILITY_INCOMPLETE`

The Phase-11 repair makes substantial, valid progress. It excludes the earlier
unbounded-outcome counterexample through a declared finite-outcome Route A,
proves a uniform Hoffman transfer for the probability-set coordinate, replaces
the pseudometric with a genuine product metric, assigns all target coordinates,
proves estimator coherence under pullback closure, and adds context allocation
and a vanishing confidence schedule.

The stopping claim `OPERATOR_LEARNABILITY_CLOSED` is still too strong. An
explicit theorem is false at its stated type, the finite-N consistency theorem
misses a zero-fiber rung case, and the engineering theorem remains about a
canonical empirical estimator rather than a parameterized trainable operator.

## Audit results

| Requirement | Judgment |
|---|---|
| Phase 0-7 freeze | `THEORY_FREEZE_CONFIRMED` |
| Phase 8-11 additive extensions | Confirmed by package isolation and chronology |
| Route-A bounded outcome/event complexity | Pass |
| Effective finite constraint geometry | Pass |
| Uniform Hoffman constant | Pass |
| Probability-coordinate metric transfer | Pass |
| OM-2 transfer in full `d_M` | Fail |
| Product metric defined on all coordinates | Pass |
| Complete, non-degenerate metric | Pass on Route-A finite-outcome class |
| Full target typing | Pass |
| OM-5-pre finite-history rung alignment | Fail for unseen positive-mass contexts |
| Full-index canonical asymptotic consistency | Supported under S1-S5 and `rho=0` |
| Finite-N full-metric bound as stated | Incomplete |
| Hidden latent/current/future information | None detected |
| Trainable-model handoff without new theorem choices | Fail |

## Decisive counterexample

OM-2 claims that Route A gives

`d_M(M_1,M_2) <= (H_bar/2) d_desc(M_1,M_2)`.

Let the index set and outcome space each contain one point. Give `M_1` and
`M_2` identical probability sets and identical rungs, but confidence
coordinates zero and one. Then

`d_desc(M_1,M_2)=0`

while

`d_M(M_1,M_2)=1`.

All Route-A dimension bounds hold. Hence bounded query outcome complexity does
not by itself transfer description convergence to the complete operator metric.
The Hoffman proof transfers only `d_H^TV` on probability sets. Confidence and
rung must be handled separately.

OM-5-pre and S5 do handle confidence separately for the canonical estimator,
which is why the repair is incomplete rather than fundamentally invalid. They
do not repair the universal statement of OM-2.

## Finite-sample consistency defect

The inherited zero-fiber contract emits a vacuous rung-1 value when `N_c=0`.
For a positive-mass context whose target assumption stack assigns rung 2 or 3,
the estimator and target rung coordinates differ by one until the first task in
that context is observed. OM-5-pre nevertheless cancels the rung coordinate at
every index, and OM-7 takes a supremum over all contexts while OM-6 exempts zero
fibers.

The almost-sure limit can still hold when the context set is finite and
`N_min -> infinity`. The displayed finite-N bound needs an explicit
all-relevant-fibers-observed condition, a probability term for missing fibers,
or a trivial bound during that event.

## Operator and target space

OM-3 successfully fixes the earlier pseudometric:

- probability objects use Hausdorff-TV distance;
- confidence lies in complete `[0,1]`;
- rungs use a discrete metric; and
- the sup product metric is non-degenerate.

OM-4 also fully types
`M^dagger=(K^dagger,1,r_decl)`. These parts pass on Route A.

The claim that every probability simplex is TV compact is not valid for
continuous outcome spaces. Thus Route A is a real restriction of the earlier
scalar/ranking interface. Continuous scalar output requires Route B or another
verified compact/stable value class; it is not covered automatically.

## Existence, identification, and learning

The three levels are correctly separated:

- ideal marked objects exist only relative to a declared marked law;
- observable data identify the outer operator `M^dagger`, leaving lift
  ambiguity explicit; and
- finite histories estimate that outer operator under the declared sampling,
  VC, atlas, confidence-schedule, geometry, and transport conditions.

No latent mark, hidden task identity, future response, or unannounced current
query label enters inference. The query is an explicit argument and current
observations enter the population arm only through `kappa`.

The theorem's learned output is `M_hat_N=A_phi(H_N)`, converging to
`M^dagger`. It does not estimate a separately defined population map also
called `A_phi`; that notation must be aligned if the required claim is
`A_hat_phi -> A_phi`.

## Engineering judgment

The contract now identifies the observable input, query-indexed population
object, current-support adaptation, decision output, uncertainty certificate,
and loss-typed abstention/failure semantics.

It still does not let an independent researcher claim learnability of a
trainable parameterized model without additional theorem choices. OM-7 proves
the canonical empirical estimator. LC-17 leaves approximation consistency of a
constrained parameterized family to instantiation, and the metric closure does
not add that theorem. Continuous outputs also require Route-B stability rather
than the default finite-outcome Route A.

## Process and model verdict

PROCESS_VERDICT: `THEORY_FREEZE_CONFIRMED`

OPERATOR_VERDICT: `OPERATOR_LEARNABILITY_INCOMPLETE`
