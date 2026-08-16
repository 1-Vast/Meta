# Decision Extension Audit

## Correct results

Phase 7 correctly proves the central separation:

- the identified set determines pointwise dominance of loss profiles;
- dominance-incomparable actions cannot be ordered from the set alone;
- a non-minimax choice requires declared decision information;
- a stochastic noise support does not define a likelihood;
- historical frequencies may change preference but not identification;
- Bayes, robust Bayes, ranking, and asymmetric decisions consume assumptions not
  present in Phase 0-6;
- minimax is recovered when the ambiguity class contains every law on the
  identified state set.

It also explicitly rejects a canonical uniform prior. No implicit stationarity,
exchangeability, or transport assumption is allowed in its stated logic.

## Hidden or invalid additions

Three defects block the claim that the generic decision interface is complete.

### 1. Loss-type mismatch

The ledger and honesty axiom compare every unconditional guarantee to
`rho_id=diameter/2`. That quantity is the frozen lower bound only for the same
scalar target under absolute error. It is not type-compatible with squared loss,
ranking 0-1 loss, abstention cost, or an arbitrary action space.

The missing object is the context-specific robust floor

$$
R_{\rm set}(J;\mathcal A,L)
=\inf_{a\in\mathcal A}\sup_{v\in J}L(a,v).
$$

### 2. Hidden tie-breaking measure

`DE-O2` introduces a full-support measure `mu_0` to select an undominated robust
minimizer. That measure is not included in `Delta`. It is decision information,
and its proof is invalid under the stated regularity: strict improvement at one
state need not have positive `mu_0` mass without state-continuity or an atom.

The operator must either return the entire argmin set, declare the tie-break, or
prove existence of an undominated selector under stronger conditions.

### 3. Generic learnability overclaim

The scalar DKW/BV and Bernoulli ranking bounds are valid in their stated scalar
regimes. They do not prove learnability of a general law on `R^m` or uniform risk
convergence for arbitrary action/loss classes. A finite-dimensional query vector
does not make its probability law finite-dimensional.

A general theorem needs a declared law class, metric, and uniform complexity
condition, or the interface must be restricted to the proved scalar/ranking
regimes.

## Additional corrections

- Exchangeability alone defines a class of predictive laws, not a unique
  predictive law. A specific exchangeable model or ambiguity class must be
  declared.
- The empty-declaration fallback for arbitrary `(A,L)` is set-robust minimization,
  not always the scalar frozen midpoint.
- Phase 7's statement that means over a bounded nonclosed scalar set sweep its
  closed convex hull is too strong at unattained endpoints; this is nonblocking
  but should be corrected.

## Decision-extension result

The separation and restricted robust/Bayes constructions are valid. The generic
operator is not yet a complete mathematical contract.
