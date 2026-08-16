# Parameterized Meta-Learning Theorem (§3)

> **Status:** Phase-13, 2026-08-03. The missing bridge from the **canonical** operator $A_\phi$ (a construction) to a **trainable parameterized family** $\{A_\theta\}_{\theta\in\Theta}$ — the audit's engineering-handoff failure. No architecture: $\Theta$ is an abstract set and every assumption is a declared capacity or regularity property. New results **PM-6–PM-9**, tagged **[proved] / [conditional] / [declared]**. Notation alignment (audit's request): $A_\phi$ denotes the canonical operator $H\mapsto\widehat M(H)$; the population target remains $M^\dagger$; the mandated claim $d_{\mathbb M}(A_\theta,A_\phi)\to0$ is **operator-to-operator, uniform over histories**, and composes with the statistical theorem $A_\phi(H_N)\to M^\dagger$ — the two convergences are different statements and are never merged.

---

## 1. The structural theorem that makes the bridge finite

**Theorem PM-6 (the canonical operator factors through a finite sufficient statistic). [proved]**
On the Route-A class, $A_\phi(H)$ depends on $H$ only through
$$s(H)\ =\ \Big(\ N_c,\ \ \textstyle\sum_i l_i(E;\iota),\ \ \sum_i u_i(E;\iota)\ \Big)_{c\in C_\kappa,\ \iota\in\mathcal I,\ E\in\mathcal E_\iota}\ \in\ S,$$
a vector of multiset counts of dimension $\le|C_\kappa|\,(1+2\,\bar e\,|\mathcal I_{/c}|)$ — finite by Route A and the finite declared index/atlas. Moreover the evaluation map $s\mapsto A_\phi(H)$ is **explicit** (empirical endpoints $\pm$ margins $\to$ polytope; schedule $\to$ confidence; fiber counts $\to$ rung/fallback) and, on each stratum of fixed integer counts, Lipschitz into $(\mathbb M,d_{\mathbb M})$: endpoints enter affinely, and endpoint-to-polytope is $\tfrac12\bar H$-Lipschitz (PM-2a); the rung/fallback coordinate is constant per stratum.
*Proof.* Inspection of the canonical construction: every quantity it forms is a function of the listed counts; the Lipschitz claim is the Hoffman transfer applied to endpoint perturbations at fixed pattern. $\square$

*Consequence.* Approximating $A_\phi$ is approximating **one explicitly known, stratum-wise Lipschitz map on a compact stratified domain** (normalized counts lie in $[0,1]^{\dim}$; strata indexed by the finite fiber-count pattern up to a declared horizon, plus a tail stratum where normalized statistics are the arguments). This — not any statistical property — is what a parameterized family must express.

## 2. The parameterized family and the approximation theorem

**Definition PM-7 (trainable family — abstract). [declared]** $\Theta\ne\emptyset$; for each $\theta$, $A_\theta:\bigcup_N\mathbb T^N\to\mathbb M$ measurable, factoring through $s(H)$ (the family reads the same statistic — a typing choice, not an architecture), with values satisfying the $\mathbb M$-constraints (coherence enforced by constructing values through the same pullback-closed constraint formation — i.e. $A_\theta$ outputs endpoint vectors, and the polytope/confidence/rung formation is shared canonical postprocessing; then $A_\theta\in\mathbb M$-valued for every $\theta$ by the Phase-11 coherence theorem applied verbatim). **[coherence proved by inheritance]**

**Assumption (P-CAP: declared expressiveness).** For every $\varepsilon>0$ there exists $\theta\in\Theta$ whose endpoint map is uniformly within $\varepsilon$ (sup over the compact statistic domain, all coordinates) of the canonical endpoint map. This is a capacity declaration about the family on a **compact finite-dimensional domain against an explicitly known target** — the weakest possible form of a universal-approximation property, and deliberately architecture-free.

**Theorem PM-8 (uniform operator approximation). [conditional on P-CAP]**
$$\inf_{\theta\in\Theta}\ \sup_{H}\ d_{\mathbb M}\big(A_\theta(H),\,A_\phi(H)\big)\ \le\ \big(\alpha\,\tfrac12\bar H+\beta L_\psi\big)\,\varepsilon\ =:\ \varepsilon',$$
where the sup is over all histories (they enter only through $s(H)$, PM-6), the $\alpha$-term is the endpoint-to-polytope transfer (PM-2a), the $\beta$-term the declared Lipschitz constant of the shared confidence schedule in its inputs, and the rung coordinate contributes $0$ (shared postprocessing ⇒ identical rungs at identical counts). Hence the mandated bound $\inf_\theta\sup_H d_{\mathbb M}\le\varepsilon'$ holds with an explicit constant, and $\varepsilon'\to0$ as the capacity $\varepsilon\to0$. $\square$

**Theorem PM-9 (total error of a trained operator; the four tiers). [conditional]**
Let $\theta_N$ satisfy $\sup_H d_{\mathbb M}(A_{\theta_N}(H),A_\phi(H))\le\varepsilon'_N+\gamma^{\mathrm{opt}}_N$ (achievable: the target is computable, so the fit is against a known oracle; $\gamma^{\mathrm{opt}}_N$ = declared optimization tolerance). Then by the triangle inequality and PM-5:
$$d_{\mathbb M}\big(A_{\theta_N}(H_N),\ M^\dagger\big)\ \le\ \underbrace{\varepsilon'_N+\gamma^{\mathrm{opt}}_N}_{\text{approximation}}\ +\ \underbrace{\alpha\tfrac12\bar H(\eta_N+\rho)+\beta\delta_N}_{\text{statistical}}\ \big[+\ (\alpha{+}\beta{+}\gamma)\ \text{on }\mathrm{Miss}_N\cup\text{bad events, with their probabilities}\big],$$
and $d_{\mathbb M}(A_{\theta_N}(H_N),M^\dagger)\to0$ a.s. when $\varepsilon'_N,\gamma^{\mathrm{opt}}_N,\delta_N\to0$, $\rho=0$. **The four tiers, separated:**
- **Existence [proved]:** $A_\phi$ exists as a construction; $M^\dagger$ exists from $\Pi_{\mathrm{obs}}$; $M^\star$ only from a declared marked law — all carried, none used below.
- **Identification [proved/scoped]:** $M^\dagger$ is the identified target; **$\theta$ itself is deliberately not identified** — only $A_\theta$'s value in $(\mathbb M,d_{\mathbb M})$ matters, and parameter gauge freedom is explicitly meaningless (the frozen gauge philosophy, one level up).
- **Statistical learning [conditional]:** PM-5 — canonical estimator to target, with the zero-fiber term.
- **Approximation learning [conditional]:** PM-8 — family to canonical operator, uniform over histories, against a computable oracle; its error is a capacity/optimization quantity carrying **no** statistical confidence and requiring **no** samples. $\square$

*Scope note (carried from the audit):* all of this is the Route-A finite-outcome class; continuous scalar outcomes require Route B's declared stability class before PM-6's compactness argument applies.
