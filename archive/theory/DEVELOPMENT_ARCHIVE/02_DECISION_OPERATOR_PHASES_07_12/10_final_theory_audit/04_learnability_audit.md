# Learnability Audit

## 1. Existence

`CONDITIONAL_PASS`

An argmin set exists under the stated compactness and lower-semicontinuity assumptions. Without those assumptions, the minimax infimum need not be attained. The repaired floor file nevertheless calls \(G_{\rm cert}\) attained without restating the existence hypotheses.

For realizable approximate computation, an \(\eta\)-argmin action has outer risk at most \(G_{\rm cert}+\eta\); its reported guarantee must include that slack or use its directly evaluated outer risk.

## 2. Identification

`PASS`

Given the exact \(J_Q(O)\), declared population ambiguity object, action set, loss, criterion, tolerance, and tie-break, the decision object is mathematically determined as an argmin set or an explicit \(\tau\)-selection. At an undeclared tie, a canonical strict singleton is intentionally not identified.

The ranking object is correctly the jointly realizable order set \(\Sigma(J_Q(O))\), not independent marginals or independently compatible pairwise signs.

## 3. Finite-history learning

`PARTIAL_PASS`

DR-L3-R proves a genuine finite-sample result under IID or conditional IID: the forced/compatible empirical endpoints plus Hoeffding and a finite union bound simultaneously cover all declared events. A separately declared concentration condition can validly replace IID if its constants are supplied. Bare exchangeability is no longer used for rates.

The result does not yet prove learning of the full decision object in every allowed case:

- simultaneous event intervals are not automatically a probability law on a strict subset of orders or on overlapping pairwise events;
- a law on the decision target alone does not determine its conditional law after observing the current support;
- the theorem for the population component does not prove that a neural \(I_\theta\) satisfies exact outer containment; the frozen Phase-5 approximation theorem supplies this only for declared tame, stability-certified classes with one-sided rounding and its other stated scope conditions.

## Separation check

The repaired documents name existence, identification, and estimation separately, so they are not textually conflated. The terminal stopping criterion nevertheless overstates closure by treating conditional validity predicates \(V_I\) and \(V_M\) as if DR-L3-R established every required trainable realization.

## Learnability verdict

`THEOREM_CHAIN_INCOMPLETE`
