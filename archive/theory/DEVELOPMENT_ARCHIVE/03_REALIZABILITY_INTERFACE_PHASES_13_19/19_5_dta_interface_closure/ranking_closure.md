# Ranking Closure (Part 4)

> **Status:** Phase-19.5, 2026-08-03. Closes audit gap 4. The determination is made, with proof; the mandated fallback is adopted. New results **DT-7–DT-8**, tagged **[proved] / [declared]**.

---

## The determination

**Theorem DT-7 (continuous affinity $\to$ ranking is impossible under the current objects). [proved]**
The Route-B scalar interface carries, per item, a marginal CDF-band class on the value interval. Order events are functionals of the **joint** law of the affinity vector, and marginals do not determine them: the frozen witnesses lift verbatim — two joint laws with identical scalar marginals and opposite order behavior (diagonal vs anti-diagonal coupling: comonotone coupling gives a deterministic tie/fixed order, antithetic coupling gives maximal order uncertainty, with every scalar marginal equal); at the law level, the equal-marginal pair of DC-R4 makes even Bayes-optimal orders differ. Hence no map from the current per-item objects (marginal band classes, however tight) to an order law is well-defined — the prohibition already typed in the earlier phases is not a caution but a theorem, and constructing the missing **joint continuous-vector band object** (multivariate CDF-band classes with coherence, stability, statistics, and a parameterization) would be a genuine theory extension, out of the mandate's scope. $\square$

## The adopted closure: ranking as a separate supervised objective

**Definition DT-8 (the ranking channel). [declared, built from existing proved objects]**
Ranking is served by the **Route-A finite order object**, with its own supervision and its own head of the operator — never derived from the scalar channel:
- **Identification tier:** for the current task, order information comes from the identified **joint** object of the queried items (the frozen $\Sigma$ machinery on $J_Q(O)$): Tier-1 sign-identified orders are certified; ambiguous orders carry the standing trichotomy and abstention semantics.
- **Population tier:** historical tasks whose own data identify (or interval-identify) order events among the deployed item-pair/list indices contribute forced/compatible order frequencies — the DC-R1 constraint polytopes on the full order space $S_m$, with the proved simultaneous-coverage statistics and LP decision brackets. This is order-*supervised* learning: the supervision is the per-task identified order set, not a transformation of scalar affinities.
- **Typed prohibition (carried, now load-bearing):** no arrow from the scalar band channel to the order channel exists in the interface; an implementation wishing to exploit cross-item dependence must first supply the missing joint object as new mathematics — the contract makes silently faking it unrepresentable.
- **Degenerate but honest case:** if the deployment's historical tasks never identify order events (no task observed the relevant items jointly), the ranking head has a vacuous population component and degrades to current-task Tier-1 identification plus abstention — the correct, certified behavior, stated in advance. $\square$

**Consequence for the verdict.** A DTA system demanding *coherent ranking derived from its continuous affinity predictions* is outside the closed scope — by theorem, not by omission. A system accepting ranking as a separately supervised objective (with the order supervision the frozen theory can actually certify) is inside it. This is the precise scope limit the final verdict must carry.
