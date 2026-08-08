# Target Typing Audit

## Result

`PASS, WITH A FINITE-SAMPLE ALIGNMENT DEFECT`

OM-4 assigns every target coordinate:

`M^dagger(iota)=(K^dagger_iota,1,r_decl(iota))`.

- `K^dagger` is the observable-law forcing/compatibility outer polytope.
- Confidence one is explicitly declared for the deterministic population
  functional.
- The rung is the deterministic rung supplied by the declared assumption stack.
- Zero-mass contexts receive a vacuous rung-1 value.

On the Route-A finite-outcome class, the target is therefore a complete element
of `M`, not merely a constraint description. The population pullback argument
is sufficient for projective coherence.

OM-5-pre nevertheless overstates finite-sample alignment. It says target and
estimator rungs cancel at every index. The inherited zero-fiber rule requires an
unobserved context to emit rung 1, while a positive-mass target context may have
declared rung 2 or 3. Until that context appears in history, the rung distance is
one. Consequently the displayed equality

`d_M(M_hat,M^dagger)=max(sup_iota d_H^TV(K_hat,K^dagger),delta)`

requires the additional event `N_c>=1` for every positive-mass relevant
context. It is not true for every finite history as stated.

This does not destroy asymptotic typing when the finite context set satisfies
`N_min -> infinity`, but it matters to the claimed finite-N bound.

