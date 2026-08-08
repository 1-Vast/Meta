# Conditional Population Object Audit

## What the declared fiber assumptions do establish

For the actual joint data-generating law, the declaration

\[
g(f)\perp \mathrm{record}\mid\kappa(\mathrm{record})
\]

is exactly a decision-sufficiency assumption. Together with conditional IID within the finite \(\kappa\)-fiber, it implies

\[
P(g(f_\beta)\in\cdot\mid O,Q)
=P(g(f_\beta)\in\cdot\mid\kappa(O),Q),
\]

provided the declaration is indexed by the query-dependent \(g,Q\). Historical tasks in the same fiber can estimate this conditional law using the repaired IID/union-bound theorem. This is a legitimate likelihood-free route, but it is a strong declared assumption, not a consequence of bounded-noise identification.

## Why DC-C2 is false as written

DC-C2 assumes `SUFF-kappa` “in the population model” and then asserts its conditional identity **for every admissible noise kernel** \(\lambda\in\Lambda_{\rm adm}\). Sufficiency involves the joint law of \((g,\mathrm{record})\), so it depends on the noise kernel. A declaration for one joint law does not imply uniform sufficiency across all bounded-support kernels.

Counterexample: let \(g\in\{0,1\}\) be fair, choose a bounded-noise tolerance large enough to admit both kernels below, and let \(\kappa\) be constant.

- Under \(\lambda_0\), the record is constant and independent of \(g\); `SUFF-kappa` holds and \(P(g\mid O)=\operatorname{Ber}(1/2)\).
- Under \(\lambda_1\), the record reveals \(g\); it remains support-admissible, but \(P(g\mid O)\) is degenerate and `SUFF-kappa` fails.

The theorem becomes valid only for the actual declared joint law, or if `SUFF-kappa` is explicitly required uniformly for every \(\lambda\) over which the conclusion quantifies.

The notation “fiber conditional restricted to \(g(I(O))\)” also needs precise typing. Under genuine sufficiency, the true fiber law is already supported on the exact identified image almost surely. For an estimated outer class, intersecting with the support constraint is a conservative class operation; it is not an additional likelihood-free reweighting identity.

## Posterior rung error

DC-C3 says that adding a declared likelihood produces a **single posterior**. A known likelihood does not collapse finite-history uncertainty about the population law. If \(\widehat{\mathcal Q}\) contains two distinct priors \(P_1,P_2\), the same likelihood generally produces two distinct posteriors. The valid rung is a posterior ambiguity class unless a singleton prior/population law is also declared or identified.

## Finite-fiber issue

The theorem uses \(\eta_{n_{\kappa(O)}}\) but does not define the finite-sample output when the current fiber has no historical members. Valid closure requires an explicit vacuous-class fallback for \(n_{\kappa(O)}=0\), or a declared positive-count condition. Conditional IID alone does not prevent an empty historical fiber.

## Verdict

`FAIL`

The fiber-sufficiency route can establish \(P(g\mid O,Q)\) under a correctly indexed actual-law declaration, but the theorem and rung ladder overclaim what their written assumptions prove.
