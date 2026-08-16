# Complete Codomain Metric and Target Typing (§2–§3)

> **Status:** Phase-11, 2026-08-03. Repairs the audit's codomain findings: the value-space metric covered only the polytope coordinate (pseudometric risk), the confidence factor $(0,1]$ was incomplete, and the target $M^\dagger$ lacked assigned confidence/rung coordinates. New results **OM-3–OM-5**, tagged **[proved] / [declared]**.

---

## 1. The complete value-space metric (§2)

**Definition OM-3 (retyped value spaces and their metric). [declared; completeness proved]**
For each index $\iota=(c,Q,\gamma)$:
$$\mathbb V_\iota\ =\ \underbrace{\mathcal C\big(\Delta(\Omega_Q)\big)}_{\text{probability object}}\ \times\ \underbrace{[0,1]}_{\text{confidence }1-\delta}\ \times\ \underbrace{\{1,2,3,4\}}_{\text{rung}},$$
with the **confidence coordinate retyped to the closed $[0,1]$** (the audit's incompleteness witness $1/n\to0$ in $(0,1]$ is thereby closed; confidence $0$ is the legitimate value of vacuous claims, consistent with the zero-fiber fallback). Metric on $\mathbb V_\iota$:
$$d_{\mathbb V_\iota}\big((K,q,r),(K',q',r')\big)\ =\ \max\Big\{\ d_H^{TV}(K,K'),\ \ |q-q'|,\ \ \mathbf 1\{r\ne r'\}\ \Big\}.$$
*Properties [proved]:* each factor is a complete metric space — the Hausdorff hyperspace of the compact $\Delta(\Omega_Q)$ under TV (compact, hence complete); $[0,1]$ Euclidean (compact); the finite rung set discrete (complete) — so the max-product is a **complete metric**, and it is a genuine metric, not a pseudometric: two values differing in *any* coordinate (including tags alone) have positive distance. Boundedness: $d_{\mathbb V_\iota}\le\max(1,1,1)=1$ (TV-Hausdorff $\le1$).

**Definition OM-3′ (operator metric, retyped).**
$$d_{\mathbb M}(M,M')\ =\ \sup_{\iota\in\mathcal I}\ d_{\mathbb V_\iota}\big(M(\iota),M'(\iota)\big)\ \in[0,1].$$
*[proved]* $(\mathbb M,d_{\mathbb M})$ is a complete metric space: uniform limits of maps into uniformly bounded complete spaces exist pointwise and converge uniformly; the defining constraints (rung-consistency, projective coherence, zero-fiber convention, admissibility $K\subseteq\Delta(\Omega_Q)$) are all preserved under uniform limits (each is a closed condition per index or per index-pair: coherence inclusions $h_*K_Q\subseteq K_{Q'}$ are closed under $d_H$ limits since $h_*$ is TV-nonexpansive and Hausdorff limits preserve inclusions). Identity of indiscernibles holds coordinate-wise — **no pseudometric remains**. Evaluation maps stay $1$-Lipschitz; the evaluation σ-algebra and measurability statements of Phase 10 carry over verbatim (the two added coordinates are measurable into their Polish factors). $\square$

## 2. The target, fully typed (§3)

**Definition OM-4 ($M^\dagger$ as a genuine element of $\mathbb M$). [declared; membership proved]**
At every index $\iota=(c,Q,\gamma)$, assign **all three coordinates**:
$$M^\dagger(\iota)\ =\ \Big(\ K^\dagger_\iota,\ \ 1,\ \ r^{\mathrm{decl}}(\iota)\ \Big),$$
- $K^\dagger_\iota$ = the population constraint polytope with endpoints $\big[\Pi_{\mathrm{obs}}(\text{forces }E\mid c),\ \Pi_{\mathrm{obs}}(\text{compatible with }E\mid c)\big]$, $E\in\mathcal E_\iota$ — nonempty (contains the pushforward of every admissible lift's conditional; LC-12(i));
- confidence coordinate $=1$: $M^\dagger$ is the deterministic $N\to\infty$ functional of $\Pi_{\mathrm{obs}}$ — a population object, not a random estimate; the confidence coordinate quantifies *sampling* uncertainty, of which it has none (identification width lives inside $K^\dagger$, where it belongs);
- rung coordinate $r^{\mathrm{decl}}(\iota)$ = the rung fixed by the **declared assumption stack** at $\iota$ (rung 3 where query- and kernel-indexed sufficiency is declared, rung 2 where only C-IID-$\kappa$, etc.) — a deterministic declared function of the index, identical for target and estimator by construction, so the rung coordinate never contributes to the estimation distance while still separating differently-tagged operators in $d_{\mathbb M}$.
*Membership [proved]:* $K^\dagger$ satisfies admissibility (subset of the simplex by construction) and projective coherence (the population forcing/compatibility indicators for a coarse event coincide with those of its pullback — the same argument as OM-5, applied at the population level); rung-consistency holds by the declared assignment; indices whose context has $\Pi_{\mathrm{obs}}$-mass zero receive the vacuous rung-1 value with confidence $1$ (consistent with the estimator's zero-fiber fallback). Hence $M^\dagger\in\mathbb M$: **the target is a true element of the operator space, all coordinates assigned, not a constraint list.** $\square$

**Alignment lemma OM-5-pre. [proved]** With OM-4's assignments, for the canonical estimator $\widehat M=A_\phi(H_N)$ (whose confidence coordinate is its declared $1-\delta$ and whose rung coordinate is $r^{\mathrm{decl}}$):
$$d_{\mathbb M}(\widehat M,M^\dagger)\ =\ \max\Big\{\sup_\iota d_H^{TV}\big(\widehat K_\iota,K^\dagger_\iota\big),\ \ \delta\Big\},$$
since the rung coordinates cancel and the confidence coordinates differ by exactly $\delta$. Thus operator convergence requires $\delta_N\to0$ along the sample path — handled in the final theorem by letting the declared confidence schedule $\delta_N\downarrow0$ slowly enough that $\eta_N(\delta_N)\to0$ still holds (e.g. $\delta_N=1/N$ changes $\eta_N$ only logarithmically). The metric forces the bookkeeping honesty the audit demanded: an estimator claiming fixed confidence $1-\delta$ never converges to the deterministic target in $d_{\mathbb M}$, and now the formalism says so. $\square$
