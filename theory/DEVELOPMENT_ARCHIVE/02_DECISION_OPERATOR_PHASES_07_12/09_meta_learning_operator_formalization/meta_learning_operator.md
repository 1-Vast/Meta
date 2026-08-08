# The Meta-Learning Operator (Part II)

> **Status:** Phase-9, 2026-08-03. Phases 0–7 frozen and cited. New results carry **ML-A** numbers, tagged **[proved] / [conditional] / [declared]**. No architectures; no implementation vocabulary. This file defines the object between historical tasks and the decision operator.

---

## 1. The codomain: transferable operator objects, not embeddings

**Definition ML-A1 (the space $\mathcal M$).** A **transferable operator object** is a map
$$M:\ C_\kappa\times\Gamma\ \longrightarrow\ \mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung},\qquad (c,\gamma)\ \mapsto\ \big(\widehat{\mathcal Q}_{c,\gamma},\ 1-\delta,\ r\big),$$
assigning to each (context, decision-specification) pair a rung-tagged, confidence-tagged **decision-information object**: a constraint class of laws on the decision-relevant outcome space (DC-R1 polytopes for ranking; interval classes for scalar pushforwards) — equivalently, by the Phase-7 primitive theorem (DE-P), a declared generator of monotone completions of the dominance order. $\mathcal M$ is the space of all such maps satisfying: rung-consistency (a value's tag matches the assumptions it was built under), the zero-fiber convention (vacuous rung-1 value where unsupported, ML-T4), and — when several $g$'s are indexed — projective consistency across pushforwards (8.1: coarser quotients receive the pushed-forward class).

$M$ is *operational*: its values plug directly into $D_\psi$'s population slot. It is evaluated, not decoded.

**Theorem ML-A2 (why the output must not be a task embedding — derived, not stylistic). [proved]**
Suppose instead $A_\phi$ output a point $z\in\mathbb R^p$ in an unconstrained representation space, to be consumed downstream. Two exhaustive cases:
(i) *$z$ is decodable*: some fixed map $\beta$ turns $z$ into an $\mathcal M$-element actually consumed by $D_\psi$. Then $\beta\circ A_\phi$ is the meta-learning operator and $z$ was an internal variable of no contractual status — the interface object is still $\mathcal M$-valued; nothing was gained by naming $z$.
(ii) *$z$ is consumed as coordinates* (some coordinate carries meaning not mediated by a validity predicate). Then the interface attaches semantics to a representation of the task/family — refuted by the frozen theory: the parametrization/gauge of the family is unidentifiable (CP §4.5: two families with identical windows differ in member sets), and continuous-summary dimension laws (F19/F20, CP-3; frozen handoff §8 item 2) bound and obstruct any faithful finite coordinateization; coordinates therefore either carry no certifiable meaning or fabricate one.
Hence a mathematically valid interface object is necessarily of type $\mathcal M$ (an operator/decision-information object with validity predicates); "no embeddings" is a theorem of the frozen corpus, not a preference. $\square$

---

## 2. The operator

**Definition ML-A3 (meta-learning operator).**
$$A_\phi:\ \bigcup_{N\ge0}\mathcal T^N\ \longrightarrow\ \mathcal M,\qquad H_N\ \mapsto\ M_{H_N},$$
with $\phi$ ranging over an unspecified approximation family (the approximable slot — never instantiated in this program).

**Validity conditions $V_A$. [declared/proved as marked]**
1. **(typed sample use) [proved from ML-T3]** Under declared task-EXCH, $A_\phi$ factors through $\mathrm{mult}(H_N)$; order-dependence is legal only under a declared drift/transport structure. $A_\phi$ never factors through $\mathrm{set}(H_N)$ — that quotient belongs to the feasibility channel alone (ML-T4).
2. **(rung-typed coverage) [conditional]** For every $(c,\gamma)$ with $N_c\ge1$: under the declared assumption stack for rung $r$ (`meta_learning_theory.md`, `conditional_population_repair.md`), $\Pr\big(\text{true rung-}r\text{ object}\in\widehat{\mathcal Q}_{c,\gamma}\big)\ge1-\delta$, with **simultaneous** coverage over the declared finite event family (union-bound constants consuming the *multiset* count $N_c$). For $N_c=0$: the vacuous rung-1 value (ML-T4) — the zero-fiber fallback is part of the definition, not an error state.
3. **(ceiling) [proved, frozen]** No value of $M_{H_N}$ ever feeds the identification channel; population output weights the inside of identified sets only (DE-H2/H3).
4. **(auditability) [declared]** $M_{H_N}$ carries its assumption tags; every downstream emission echoes them (H6).

**Remark (what $A_\phi$ is, in frozen vocabulary).** Phase 3 proved: *meta-learning is identification of the window system; adaptation is sectioning it.* Phase 9 adds the decision layer's parallel: **meta-learning of decision information is estimation of the context-indexed population assignment; adaptation is evaluation at the current context plus support restriction.** The two layers share the archive but consume disjoint quotients of it (ML-T4) — the same separation, one level up.
