# The Final Decision Operator (Parts VII + synthesis)

> **Status:** Phase-7, 2026-08-03. Frozen corpus cited, not modified. New results carry **DE-O** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. This file assembles the phase into a single mathematical interface: a decision operator that selects actions under residual ambiguity **without ever claiming the ambiguity has disappeared**. No architectures, no parameterizations, no training losses (mandate XIII).

---

## 1. Signature

$$\mathbb D:\ \underbrace{\big(J_Q(O),\ \rho_{\mathrm{id}},\ \text{flags}\big)}_{\text{frozen identification output}}\ \times\ \underbrace{\Delta}_{\text{declared decision information}}\ \times\ \underbrace{(\mathcal A,L)}_{\text{decision context}}\ \longrightarrow\ \big(a^*,\ \mathrm{Ledger}\big)$$

with $\Delta=\big(\mathcal Q\ \text{[ambiguity class of laws / single law / }\emptyset\text{]},\ \lambda\text{-declaration or }\Lambda_{\mathrm{adm}},\ \text{transport class},\ \text{criterion},\ \text{tie convention}\big)$ — every slot **declared, auditable, and possibly empty**. The first argument is the frozen operator's output (Phase-6 conditional/union form included: any nonempty admissible set plugs in, per the Part-I header). The Ledger is the four-row object of `decision_uncertainty_ledger.md`.

**Composition (definition DE-O0).**
$$\mathbb D=\text{(frozen }I\text{, untouched)}\ \to\ \text{pushforward }e_Q\ \to\ \text{condition/restrict }\mathcal Q\ \text{to}\ \Delta(J_Q(O))\ \text{(DE-H2)}\ \to\ \text{criterion }\rho_{\mathcal Q}\ \to\ \text{undominated selection}\ \to\ \text{Ledger emission}.$$

---

## 2. Honesty axioms

- **H1 (separation).** The identification component of the output is a function of $O$ alone: for all $\Delta,\Delta'$, the reported $(J_Q(O),\rho_{\mathrm{id}},\text{flags})$ coincide. Identified set and decision preference remain **separate output objects** — the mandate's Part VII requirement, made an axiom.
- **H2 (ceiling).** Every unconditional guarantee emitted is $\ge\rho_{\mathrm{id}}$; any tighter risk statement carries the tags of the axioms that bought it (DE-L5).
- **H3 (support).** Every weighting used is supported in $I(O)$ (automatic under the (LIK) support condition, DE-H2; enforced as a validity check otherwise).
- **H4 (fallback).** $\Delta=\emptyset\ \Rightarrow\ \mathbb D=$ the frozen canonical operator (Chebyshev center + radius); and continuously, $\mathcal Q\uparrow\Delta(J_Q(O))\Rightarrow\mathbb D\to$ the frozen minimax rule (DE-T4(iii)). The extension is **conservative over the frozen theory**.
- **H5 (admissibility).** $a^*$ is undominated on $J_Q(O)$ (never contradicts the identified set).
- **H6 (ledger completeness).** The output always contains all four ledger rows with epistemic tags, the ranking tier where applicable, and **$\Delta$ itself** (the declared information is part of the answer, so the answer is auditable and the ambiguity attributable).

---

## 3. Theorems

**Theorem DE-O1 (reduction at empty declaration). [proved]**
Under H4's definition, with $\Delta=\emptyset$ the only available monotone completions are the set-canonical ones (DE-P6), among which the frozen minimax selection is the unique certificate-bearing choice (Theorem 1); $\mathbb D$ then outputs exactly the frozen pair (center, radius). Nothing in Phase 7 alters any Phase-0–6 output. $\square$

**Theorem DE-O2 (existence and admissible selection). [proved]**
Let $\mathcal A$ be compact metric, $L$ bounded, l.s.c. in $a$, $\mathcal Q$ nonempty convex weak-*-compact (conditioned per H3). Then:
(i) the robust risk $\rho_{\mathcal Q}(a)=\sup_{\mu\in\mathcal Q}\int L(a,\cdot)\,d\mu$ attains its minimum (DE-T4(i)); the argmin set $M$ is nonempty and closed;
(ii) an **undominated** minimizer exists: fix any full-support $\mu_0$ on $J_Q(O)$ (or a countable-support surrogate on a dense subset) and take $a^*\in\arg\min_{a\in M}\int L(a,\cdot)\,d\mu_0$ — if some $a'$ dominated $a^*$ strictly somewhere, then $\rho_{\mathcal Q}(a')\le\rho_{\mathcal Q}(a^*)$ (so $a'\in M$) while $\int L(a')\,d\mu_0<\int L(a^*)\,d\mu_0$, a contradiction. So H5 is satisfiable whenever the criterion is well-posed; with strictly convex risks the argmin is a singleton and the step is vacuous. $\square$

**Theorem DE-O3 (selection is epistemically inert). [proved]**
By construction (H1), $\mathbb D$'s selection step reads, and does not write, the identification layer: for every $\Delta$, post-decision $I(O)$, $J_Q(O)$, $\rho_{\mathrm{id}}$, and all partiality flags are unchanged, and remain part of the output. Consequently **an action may be selected while the member — even the sign, even the value — remains unidentified**, with the residual ambiguity stated rather than resolved; and conversely no selection, however confident its conditional risk, ever narrows what is identified. The mandated meaningful post-decision uncertainty is exactly the Ledger: the unconditional row ($\rho_{\mathrm{id}}$) because it survives all declarations; the conditional rows ($\sigma_{\mathrm{prob}},r_{\mathrm{dec}},u_{\mathrm{rank}}$-tier) because their tags say *which* assumptions would have to hold for them to bind. $\square$

**Theorem DE-O4 (the operator never claims disappearance of ambiguity — formal version). [proved]**
Suppose an output claimed residual uncertainty $<\rho_{\mathrm{id}}$ unconditionally. By H2 this violates the axioms; by DE-L5(i) it also violates Theorem 1 — so within this interface, honesty is not a policy but a consistency condition: **a dishonest output is mathematically inconsistent with the frozen corpus, not merely disallowed.** $\square$

---

## 4. Output schema (the interface, complete)

$$\mathbb D(\cdots)\;=\;\Big(\ a^*\ ;\ \underbrace{J_Q(O)\ \text{or its envelope},\ \rho_{\mathrm{id}},\ \text{flags}}_{\text{identification, verbatim (H1)}}\ ;\ \underbrace{\sigma_{\mathrm{prob}}\,[\text{tags}],\ r_{\mathrm{dec}}(a^*)\,[\text{tags}],\ u_{\mathrm{rank}}\text{-tier}}_{\text{conditional rows}}\ ;\ \underbrace{\Delta}_{\text{declared inputs, echoed}}\Big)$$

Degenerate and boundary behavior, fixed by the axioms: empty $\Delta$ → frozen operator (DE-O1); unbounded sections / off-coverage queries → $\rho_{\mathrm{id}}=+\infty$ is reported even when a declared prior yields finite $\sigma_{\mathrm{prob}}$ (DE-U2 is the case this schema exists to survive); ranking at Tier 3 → the strict output is emitted **only if** a tie convention was declared, else abstention (DE-R3(iii)); violated declarations detected via empty conditioned sections → the frozen misspecification flag fires (unrealizability semantics, unchanged).

**Relation to the frozen interface.** Phase 5's operator emitted (center, radius, flags); $\mathbb D$ is the unique extension of that interface which (a) accepts declared decision information, (b) is conservative at the empty declaration, (c) preserves the ceiling, and (d) emits its own audit trail. Uniqueness here means: any interface satisfying H1–H6 with the same slots agrees with $\mathbb D$ up to the choice of completions inside $\Delta$ — which is exactly the freedom Part II proved irreducible.
