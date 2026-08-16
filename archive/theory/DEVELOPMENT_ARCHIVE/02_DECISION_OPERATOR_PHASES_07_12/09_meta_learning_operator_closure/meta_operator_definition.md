# The Meta-Learning Operator (§2)

> **Status:** Phase-9 closure, 2026-08-03. New results carry **MC-** numbers, tagged **[proved] / [conditional] / [declared]**. No architectures, losses, or implementation vocabulary.

---

## 1. Codomain first: the transferable operator space $\mathbb M$

**Definition MC-6.** A **transferable operator object** is a map
$$M:\ C_\kappa\times\Gamma\ \longrightarrow\ \mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung},\qquad (c,\gamma)\mapsto\big(\widehat{\mathcal Q}_{c,\gamma},\,1-\delta,\,r\big),$$
assigning to every (context, decision-specification) pair a **decision-information object**: a rung-tagged, confidence-tagged constraint class of laws on the decision-relevant outcome space — full-space ranking polytopes (DC-R1), interval classes for scalar pushforwards — equivalently, by the Phase-7 minimal-primitive theorem (DE-P), a declared generator of monotone completions of the dominance order. $\mathbb M$ is the set of such maps satisfying: rung-consistency (tags match the assumptions used), projective consistency across pushforward quotients, and the zero-fiber convention (MC-5).

**Theorem MC-7 (what $M$ must not be — three exclusions, each derived). [proved]**
(i) **Not an embedding / latent vector.** If a vector $z\in\mathbb R^p$ were the interface object: either $z$ is decodable by a fixed map into an $\mathbb M$-element actually consumed downstream — then the composite is the operator and $z$ is contractually void; or some coordinate of $z$ carries meaning not mediated by a validity predicate — refuted by the frozen theory: the family's gauge/parametrization is unidentifiable (CP §4.5) and continuous-summary dimension laws (F19/F20, CP-3) obstruct faithful finite coordinatizations. Unconstrained coordinates either mean nothing certifiable or fabricate.
(ii) **Not a task ID.** An identifier is decision-relevant only through a lookup table; the table *is* the object, and the pair (ID, table) reduces to case (i)'s decodable branch. Moreover an ID scheme is non-transferable by construction: it assigns nothing to unseen contexts, violating the totality of $\mathbb M$-elements (every $(c,\gamma)$ receives at least the rung-1 vacuous value).
(iii) **Positive characterization.** The interface object must carry decision-relevant information *with validity predicates attached* — and the typed classes of MC-6 are exactly the objects for which the Phase-7/8 validity theorems exist. Hence $\mathbb M$-typing is forced, not stylistic. $\square$

## 2. The operator

**Definition MC-8 (meta-learning operator).**
$$A_\phi:\ \bigcup_{N\ge0}\mathbb T^N\ \longrightarrow\ \mathbb M,\qquad H_N\mapsto M_{H_N},$$
$\phi$ ranging over an unspecified approximation family (the approximable slot; never instantiated here).

**Validity conditions $V_A$:**
1. **(typed sample use) [proved, MC-3/4]** Under declared task-exchangeability $A_\phi$ factors through $\mathrm{mult}(H_N)$; order-use only under declared drift; $A_\phi$ never factors through $\mathrm{set}(H_N)$.
2. **(rung-typed coverage) [conditional]** For each $(c,\gamma)$ with $N_c\ge1$: under the declared rung-$r$ assumption stack (§4–§5), the value covers the true rung-$r$ object with probability $\ge1-\delta$, simultaneously over the declared finite event family, at multiset count $N_c$; for $N_c=0$: the vacuous rung-1 value.
3. **(ceiling) [proved, frozen]** No value of $M_{H_N}$ reaches the identification channel (DE-H2/H3).
4. **(auditability) [declared]** Every value carries its assumption tags; downstream emissions echo them (H6).

**Canonical construction (the existence witness for the definition).** The forced/compatible interval polytope estimator: for each event $E$ in the declared family and fiber $c$, bounds $\big[\tfrac1{N_c}\sum l_i(E),\ \tfrac1{N_c}\sum u_i(E)\big]\pm\eta_{N_c}$ from each historical task's own identified object, assembled into the DC-R1 constraint class. This is a concrete, assumption-tagged element of $\mathbb M$ for every $H_N$ — so $A_\phi$'s codomain is inhabited constructively, not just axiomatically. [proved, from the repaired 8.1/9 learning theorems]
