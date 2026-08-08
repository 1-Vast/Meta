# Conditional Population Information (§4)

> **Status:** Phase-9 closure, 2026-08-03. This file proves **Route B** (sufficient statistic) in full, with the kernel-indexed quantifier that the `../11_final_closure_audit/` demanded; Routes A and C are stated as declared alternatives and shown to normalize to B. **No marginal-to-conditional leap occurs anywhere**: every conditional claim is either proved from a declared conditional-independence assumption or forbidden by the rung typing. New results **MC-11–MC-15**, tagged **[proved] / [conditional] / [declared] / [impossible]**.

---

## 1. The probability space (explicit, so the quantifiers are checkable)

Within a context fiber $c$: latent member $f\sim\pi(\cdot\mid c)$; observable record $R$ generated from $f$ by a noise kernel $\lambda$ with the frozen bounded-support property; decision target $g(f)$ for the query-indexed pushforward $g$ of $\gamma$. The joint law is $\Pi_c\otimes\lambda$ on (member, record). Declared objects:

- **($\kappa$-DESIGN) [declared]** measurable $\kappa:\ \text{records}\to C_\kappa$ finite, computable from any task's record (historical or current).
- **(C-IID-$\kappa$) [declared]** conditional on $\kappa=c$, the (record, member) pairs — historical and current — are IID under the actual generating law. *Task-level only; the within-task noise stays adversarial-bounded (no within-task stochastic assumption).*
- **(SUFF-$\kappa$-$\lambda$) [declared, kernel-indexed]** under $\Pi_c\otimes\lambda$: $\ g(f)\ \perp\ R\ \mid\ \kappa(R)$. Declared either for the *actual declared kernel* $\lambda^*$, or *uniformly* over $\Lambda_{\mathrm{adm}}$ (SUFF-$\kappa$-U). The index is load-bearing (audit witness, carried): with constant $\kappa$, an uninformative admissible kernel satisfies the independence while a revealing one violates it — sufficiency never transfers across kernels for free.

## 2. Route B, proved

**Theorem MC-11 (conditioning identity — Route B). [proved from the declared assumptions]**
Under ($\kappa$-DESIGN) + (C-IID-$\kappa$) + (SUFF-$\kappa$-$\lambda^*$), for the actual law and every bounded measurable $h$:
$$\mathbb E\big[h(g(f))\ \big|\ R\big]\ =\ \mathbb E\big[h(g(f))\ \big|\ \kappa(R)\big]\quad\text{a.s.},\qquad\text{i.e.}\qquad P\big(g(f)\in\cdot\mid O_*,Q\big)=P\big(g(f)\in\cdot\mid\kappa(O_*),Q\big).$$
*Proof.* $\kappa(R)$ is $\sigma(R)$-measurable, so $\sigma(\kappa(R))\subseteq\sigma(R)$; conditional independence $g(f)\perp R\mid\kappa(R)$ states precisely that for bounded $h$, $\mathbb E[h(g(f))\mid \sigma(R)]$ is $\sigma(\kappa(R))$-measurable and equals $\mathbb E[h(g(f))\mid\kappa(R)]$ a.s. (tower property applied to the definition of conditional independence between $\sigma$-algebras with one containing the conditioning algebra). No marginal object appears: the left side *is* the current-observation conditional, and the identity replaces it by a fiber conditional **only because the declared independence says so**. $\square$

**Theorem MC-12 (estimability of the right-hand side). [conditional]**
Under (C-IID-$\kappa$), the historical tasks in fiber $c_*=\kappa(O_*)$ are IID draws from the same conditional population; the forced/compatible interval polytope at multiset count $N_{c_*}$ with margins $\eta_{N_{c_*}}$ (the repaired simultaneous-coverage theorem) is a valid $1-\delta$ confidence class for $P(g(f)\in\cdot\mid\kappa(O_*))$ — hence, by MC-11, **for the current-observation-conditioned object**. Zero fiber: MC-5 fallback. $\square$

**Theorem MC-13 (support consistency — and a falsification signal, not a patch). [proved]**
The true member always lies in its own identified set, so $P(g(f)\in g(I(R))\mid R)=1$ a.s.; combined with MC-11, (SUFF-$\kappa$-$\lambda^*$) *implies* that the fiber conditional is a.s. supported inside the current identified image. Consequently: if the estimated fiber class places material mass outside $g(\widehat J)$ (beyond its confidence slack), this is **evidence against the sufficiency declaration** — the support-restriction step then functions as a declared projection and must be flagged, not silently applied. Support restriction itself remains likelihood-free and always valid (frozen DE-H2); what MC-13 adds is that *needing it badly* is a SUFF-$\kappa$ audit statistic. $\square$

**Theorem MC-14 (posterior classes; the ladder; the floor). [proved / impossible]**
(i) A declared likelihood maps the prior class to the **posterior class** $\{P_0(\cdot\mid O_*;\lambda^*):P_0\in\widehat{\mathcal Q}\}$ — a singleton iff the prior class is (the 8.2 "single posterior" remains retracted).
(ii) The rung ladder stands: nothing / zero fiber → rung 1 (all laws on the identified support; frozen minimax endpoint); (C-IID-$\kappa$) → rung 2, **marginal-typed**, consumable only by decisions declared insensitive to residual-record dependence; + (SUFF-$\kappa$-$\lambda^*$ or U) → rung 3, the proved conditional (MC-11/12); + declared $\lambda^*$ on a class → rung 4, posterior class.
(iii) **[impossible]** Below rung 3 the gap between the $\kappa$-conditional and the $O_*$-conditional is not estimable within the frozen noise model (DE-H4: one prior, many posteriors). Any interface consuming a rung-2 object at rung 3 makes the forbidden marginal-to-conditional leap — excluded here by the type system, not by good intentions. $\square$

## 3. Routes A and C — declared alternatives, same normal form

**Route A [conditional].** Learn the joint law $P(g(f),\,s(R),\,Q)$ for a declared observable statistic $s$ evaluable on the current record, then condition on the observed $s$-cell. Requirements: the strong per-task gate (each historical task's own data (interval-)identify *both* $s$ and $g$), a declared positive-mass bound $q_0$ for the conditioning cell, rates degraded by $1/q_0$. Formally this is Route B with $\kappa=$ the declared discretization of $s$ and (SUFF-$\kappa$) holding by construction when the cell *is* the conditioning event — A normalizes to B.
**Route C [declared].** An explicit conditional fiber model: declare the family $\{\pi(\cdot\mid c)\}_{c\in C_\kappa}$ up to declared parameters and estimate within fibers — again B's structure with a stronger parametric declaration replacing the nonparametric interval class.
**Chosen route: B** — the weakest declared machinery for which the conditioning identity is a theorem (MC-11); A and C are permitted strengthenings, and both are consumed through the same rung typing.
