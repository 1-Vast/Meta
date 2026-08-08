# Model Compilation Readiness

## Decision

C.

THEORY_EXTENSION_REQUIRED

## Why explicit assumptions alone are not enough

The program correctly identifies the assumptions needed for population use, but
three gaps are internal to the mathematical interface rather than choices a
model builder can safely fill:

1. The information floor and fallback are not typed to arbitrary `(A,L)`. A
   theorem must replace the universal half-diameter comparison by
   `inf_a sup_v L(a,v)` and recover half-diameter only as the scalar absolute-loss
   corollary.
2. The final operator lacks a declared general tie-break and its undominated
   selection proof imports an undeclared measure under insufficient regularity.
   The output must be set-valued or a valid selection theorem must be supplied.
3. Learnability is proved for scalar CDF/BV and Bernoulli ranking objects, not for
   the generic multivariate joint law and arbitrary loss class advertised by the
   interface. A uniform convergence theorem or a formal scope restriction is
   required.

## Central pipeline

The pipeline

$$
\text{history}\to\text{population information}
\to\text{decision under ambiguity}\to\text{current action}
$$

is mathematically well-defined in restricted contexts when the following are
declared: current identified set, finite query pushforward, action space, loss,
population bridge, likelihood or likelihood ambiguity, population-law ambiguity
class, transport class under shift, coverage/censoring mechanism, and tie-break
or set-valued output.

That restricted existence does not establish the generic Phase-7 handoff claim.
The exact missing theorem package is the type-correct robust-floor/fallback
theorem, the honest selection theorem, and the joint-object uniform learnability
theorem.
