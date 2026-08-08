# Query-Indexed Meta-Object (§2)

> **Status:** Phase-10, 2026-08-03. Repairs audit finding 1 (`../11_meta_learning_final_audit/`): the Phase-9 transferable object $M(\kappa,\gamma)$ dropped the separately defined query argument — untyped composition. New results **LC-5–LC-7**, tagged **[proved] / [declared]**.

---

## 1. The repaired object

**Definition LC-5 (query-indexed transferable object).**
$$M:\ C_\kappa\times\mathcal Q_0\times\Gamma\ \longrightarrow\ \bigsqcup_{Q\in\mathcal Q_0}\Big(\mathfrak Q(\Omega_Q)\times(0,1]\times\mathrm{Rung}\Big),\qquad (c,Q,\gamma)\ \mapsto\ \big(\widehat{\mathcal Q}_{c,Q,\gamma},\ 1-\delta,\ r\big),$$
with the **typing constraint** $M(c,Q,\gamma)\in\mathfrak Q(\Omega_Q)\times(0,1]\times\mathrm{Rung}$ — the value's outcome space is the one belonging to *its* query index ($\Omega_Q$: the order space $S_{m_Q}$ for ranking specifications, value/difference spaces for prediction). Equivalent representation (permitted by the mandate): absorb $Q$ into $\gamma$ as a distinguished component, **with the coherence condition that the value's outcome space is determined by that component** — the representations are mathematically identical; what is forbidden is any form in which two specifications differing only in $Q$ receive one value.

**Coherence across queries (part of the definition).** For nested or overlapping queries, values are **projectively consistent as outer classes**: if $Q'\subseteq Q$ with quotient map $h:\Omega_Q\to\Omega_{Q'}$ (restriction of orders, coordinate projection), then $h_*\widehat{\mathcal Q}_{c,Q,\gamma}\subseteq\widehat{\mathcal Q}_{c,Q',\gamma'}$ — a finer query may not contradict a coarser one; outer slack is permitted (validity is one-sided, as everywhere in this program).

## 2. Different queries provably require different outputs

**Theorem LC-6 (the query index is load-bearing — audit witness, adopted and proved). [proved]**
There is a population and a single (context, specification) pair for which two queries force different values; hence no $Q$-free object can represent the target.
*Witness.* Let the population be degenerate at one member with values $v(x_0)=2$, $v(x_1)=0$, $v(x_2)=1$ (constant context $c$; specification $\gamma$ = pairwise $0$–$1$ ranking of the first-listed query point against the second). Then the true population value at $Q=(x_0,x_1)$ is "first wins with probability $1$", and at $Q'=(x_1,x_2)$ it is "first wins with probability $0$": $\ P\big(v(x_0)>v(x_1)\big)=1\ne0=P\big(v(x_1)>v(x_2)\big)$. A map receiving only $(c,\gamma)$ receives identical inputs for these two required outputs and can represent at most one — contradiction. So $Q$ must appear as an explicit index (or as the distinguished $\gamma$-component of LC-5's equivalent form). $\square$

## 3. The query survives composition

**Typing rule LC-7 (no arrow erases $Q$). [declared, enforced; verified in §6]**
Adaptation-time evaluation is at the full index: $\Delta_{\mathrm{pop}}=[A_\phi(H_N)]\big(\kappa(O_*),\ Q_*,\ \gamma_*\big)$, and its value lives in $\mathfrak Q(\Omega_{Q_*})\times(0,1]\times\mathrm{Rung}$ — the same $\Omega_{Q_*}$ on which the identification arm's order projection $\widehat\Sigma_{Q_*}$ and the decision operator's loss are typed. Every subsequent arrow (support restriction, criterion evaluation, ledger emission) is $Q_*$-typed; a composition step whose domain or codomain forgets $Q_*$ is ill-formed. The Phase-9 conditioning theorem composes at the same index: the sufficiency declaration and the fiber estimate are **query-indexed** ($g_{Q_*}(f)\perp O\mid\kappa(O)$, per query specification), exactly as the audit's conditional-information section required. $\square$
