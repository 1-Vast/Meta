# Estimator Coherence (§4)

> **Status:** Phase-11, 2026-08-03. Repairs the audit's finding that estimator membership in $\mathbb M$ was unproved: independently built confidence polytopes might violate projective coherence, and the fiber-wise confidence allocation across contexts was unaccounted. New results **OM-5–OM-6**, tagged **[proved] / [declared]**.

---

## 1. The coherent-atlas construction rule

**Declared construction rule (pullback closure).** The event atlas is built so that for every declared coarsening $h:\Omega_Q\to\Omega_{Q'}$ ($Q'\subseteq Q$: order restriction, coordinate projection) and every coarse event $E'\in\mathcal E_{(c,Q',\gamma')}$, the pullback $h^{-1}(E')$ belongs to $\mathcal E_{(c,Q,\gamma)}$; and margins at paired events are computed from the same fiber count and the same declared $\delta$-schedule.

**Theorem OM-5 (the canonical estimator is coherent — membership in $\mathbb M$ proved). [proved]**
Under the pullback-closure rule, the canonical estimator's values satisfy, for every coarsening $h$:
$$h_*\widehat K_{(c,Q,\gamma)}\ \subseteq\ \widehat K_{(c,Q',\gamma')},$$
together with admissibility and rung/zero-fiber consistency — so $A_\phi(H_N)\in\mathbb M$, genuinely.
*Proof.* The key identity is at the per-task indicator level. Let $\Sigma^{(i)}_Q$ be task $i$'s identified object at $Q$; the frozen projection property gives $\Sigma^{(i)}_{Q'}=h\big(\Sigma^{(i)}_Q\big)$ (restriction of the same admissible member set). Hence for a coarse event $E'$:
$$\text{forces: }\ \Sigma^{(i)}_{Q'}\subseteq E'\iff h(\Sigma^{(i)}_Q)\subseteq E'\iff\Sigma^{(i)}_Q\subseteq h^{-1}(E');\qquad \text{compatible: }\ \Sigma^{(i)}_{Q'}\cap E'\ne\emptyset\iff\Sigma^{(i)}_Q\cap h^{-1}(E')\ne\emptyset.$$
So the per-task indicators — and therefore the empirical endpoints and (by the shared schedule) the margins — **coincide exactly** between the coarse event $E'$ and its pullback $h^{-1}(E')$. Now take $p\in\widehat K_{(c,Q,\gamma)}$: it satisfies the pullback constraint $\ell\le p(h^{-1}E')\le u$; since $(h_*p)(E')=p(h^{-1}E')$, the pushforward satisfies the identical coarse constraint; this holds for every $E'$ in the coarse atlas (pullback closure), and $h_*$ maps the simplex to the simplex — so $h_*p\in\widehat K_{(c,Q',\gamma')}$. Admissibility is by construction (simplex rows in every constraint system); rung and zero-fiber coordinates follow the declared assignments (OM-4). Nonemptiness: the empirical fiber law's pushforwards satisfy all constraints simultaneously (its evaluations *are* the constraint centers), witnessing every value nonempty — and coherently so, since pushforward of the empirical law at $Q$ is the empirical law at $Q'$. $\square$

*Remark.* Coherence is thus a property of the **atlas discipline**, not a lucky accident of the data: identical indicators at paired events make the coarse constraints logically implied by the fine ones. An atlas violating pullback closure loses the proof — which is why the rule is part of the declared construction, echoed in the ledger.

## 2. Confidence allocation across contexts (the audit's secondary accounting gap)

**Theorem OM-6 (context-complete simultaneous coverage). [conditional on the declared stack]**
Let the declared index class have $|C_\kappa|$ contexts and per-index atlases within the (RA1) bounds. Allocate $\delta$ as: $\delta_c=\delta/|C_\kappa|$ per context (or include the context-conditioned indicators $\mathbf 1\{\kappa(r)=c\}\cdot\mathbf 1\{\ldots\}$ in the declared VC class, enlarging $d^*$ by the standard composition bound). Then with probability $\ge1-\delta$, **simultaneously over all contexts, queries, specifications, and events**:
$$\max_{c\,:\,N_c\ge1}\ \sup_{(Q,\gamma)}\ \sup_{E}\ \big|\widehat{\text{endpoint}}-\text{endpoint}\big|\ \le\ \eta_N\ :=\ C\sqrt{\frac{d^*\ln(N_{\min}+1)+\ln(|C_\kappa|/\delta)}{N_{\min}}},$$
the $\ln|C_\kappa|$ term making the previously implicit maximization over fibers explicit and legitimate. Zero fibers are exempt (vacuous values carry no coverage claim). $\square$
