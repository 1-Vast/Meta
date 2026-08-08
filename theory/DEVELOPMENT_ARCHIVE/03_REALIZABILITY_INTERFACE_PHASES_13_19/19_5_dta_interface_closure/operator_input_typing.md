# Operator Input Typing (Part 1)

> **Status:** Phase-19.5, 2026-08-03. Closes audit gap 1 (history/support mismatch) by choosing **exactly one** input type and proving its consistency. New results **DT-1–DT-2**, tagged **[proved] / [declared]**.

---

## The choice: **B** — history compressed into a declared state

$$\boxed{\ A_\theta\big(z_H,\ S,\ Q,\ \gamma\big)\ }$$

with the **deployment state** $z_H$ produced once, at the end of meta-training, and frozen:
$$z_H\ =\ \Big(\ \hat\theta\ \text{(trained parameters)},\ \ \big(N_c\big)_{c\in C_\kappa}\ \text{(fiber counts)},\ \ \big(b^{\mathrm{pop}}_c\big)_{c}\ \text{(population bands per context, with margins)},\ \ \text{tags}\ \Big),$$
and the current-task arguments: support $S$ (finite set, MI-2), from which the operator computes the **identification object** — the exact identified interval/order data $I(S)$ **as an explicit typed component** (DT-0.5) and the context $\kappa(S)$ — plus the index $(Q,\gamma)$.

**Theorem DT-1 (type consistency — why B, and why B is forced). [proved]**
(i) *B is sufficient:* the deployed operator's dependence on the raw history $H_N$ factors exactly through $z_H$: the Phase-17 family's outputs depend on $H$ only through the trained parameters (which fix $\lambda,b_j$), the fiber counts (which fix confidence margins and rungs), and the per-context population bands (the canonical $\lambda{=}0$ component at the current context). Nothing else in $H$ is ever read at deployment — so $z_H$ is the sufficiency quotient of $H$ for the deployed map, and option **A** ($H$ explicit at deployment) carries provably unused type weight (and would re-open the multiset-transport question at every call).
(ii) *C is the degenerate case, not an alternative:* with empty history, $z_H$ is the vacuous state (zero counts) and the operator degrades by the standing fallback — rung-1 population component, frozen minimax semantics; option **C** is B at $z_H=\varnothing$, so choosing C would discard the population channel that the whole program exists to type.
(iii) *The audited mismatch is resolved by splitting the conflated symbol:* the **population band** $b^{\mathrm{pop}}_{\kappa(S)}$ (from $z_H$ — historical counts, historical margins) and the **identification band/object** $b_{\mathrm{id}}(S), I(S)$ (from the current support — no counts, no margins, exact) are different objects with different types and different guarantees; the deployed assembly is
$$b_\theta\ =\ (1-\lambda)\,b^{\mathrm{pop}}_{\kappa(S)}\ +\ \lambda\sum_j\varphi_j\big(z(S;Q,\gamma)\big)\,b_j,\qquad \text{output}\ =\ \Big(K(b_\theta)\big|_{\,\mathrm{supp}\,I(S)},\ \ \mathrm{confidence}(z_H),\ \ \mathrm{rung}(z_H),\ \ \text{certificate row from }I(S)\Big),$$
with support restriction reading $I(S)$ directly (likelihood-free, frozen DE-H2) and the certificate row computed from the current identification alone — $\theta$- and $z_H$-invariant, as the honesty axioms require. Confidence and rung are functions of $z_H$'s counts — well-typed now, because the counts *exist* in the input (the audit's precise objection). $\square$

**Definition DT-2 (the four previously-ambiguous alternatives, decided). [declared]**
Historical counts: retained in $z_H$ as learned-state metadata — alternative (i) of the audit, adopted. The operator does not receive raw $(H_N,\dots)$ at deployment (alternative (ii) rejected as non-minimal). $b_{\mathrm{can}}(S_T)$-as-population-band (alternative (iii)) is retracted — the current support never manufactures population information. Zero-history default (alternative (iv)) is the $z_H=\varnothing$ fallback, inherited, with its flag. One type, one semantics, no inference left to the implementer.
