# Operator Learnability Theorem (§5)

> **Status:** Phase-10, 2026-08-03. The learning tier at operator level — repairing the audit's "pointwise concentration is not operator learning". New results **LC-15–LC-17**, tagged **[proved] / [conditional] / [declared]**. The theorem controls a supremum over the **entire** valid index class $\mathcal I=C_\kappa\times\mathcal Q_0\times\Gamma_0$ and its full event atlas — not a finite manually selected family.

---

## 1. Assumptions (each declared, each priced)

- **(A1) Task sampling.** Task-level **(IID)** — or **(C-IID-$\kappa$)** with multiset fiber counts $N_c$, the relevant fibers having $N_c\ge N_{\min}$. Task-exchangeability alone is admitted only at its true strength (mixture reformulation; **no rates** — the shared-Bernoulli witness is carried). Within-task noise remains the frozen adversarial model (no within-task IID anywhere).
- **(A2) Complexity control.** The **induced event-indicator class**
$$\mathcal F^*=\Big\{r\mapsto\mathbf 1\{\text{identified object of record }r\ \text{forces }E\},\ \ r\mapsto\mathbf 1\{\text{compatible with }E\}\ :\ E\in\mathcal E_{Q,\gamma},\ (Q,\gamma)\in\mathcal Q_0\times\Gamma_0\Big\}$$
has declared finite VC dimension $d^*$ (equivalently a declared polynomial covering bound). *Automatic case:* a finite atlas gives $d^*\le\log_2|\mathcal F^*|$; for infinite atlases the declaration is substantive regularity of the interplay between task designs and the query/specification class, and it is echoed.
- **(A3) Concentration.** Bounded indicators (automatic under (A1)); any declared dependence structure must bring its own declared inequality with constants.
- **(A4) Transport.** Declared radius $\rho$ if historical and current populations differ ($\rho=0$ under (IID) to the current population).

**The operator-level deviation metric (what is actually estimated and consumed).** For $M,M'\in\mathbb M$ with values described by eventwise constraint intervals,
$$d_{\mathrm{desc}}(M,M')\ =\ \sup_{\iota\in\mathcal I}\ \sup_{E\in\mathcal E_\iota}\ \max\big(|l\text{-endpoint gap}|,\ |u\text{-endpoint gap}|\big),$$
the uniform constraint-description metric — finer-grained than $d_{\mathbb M}$-via-$d_H$ for the estimator's purposes, since every downstream consumption (LP risk brackets, robustness checks, $\Gamma$-minimax) reads the constraints.

## 2. The theorem

**Theorem LC-15 (uniform finite-task operator learning). [conditional on (A1)–(A4)]**
With probability $\ge1-\delta$, **simultaneously over the entire index class and event atlas**,
$$d_{\mathrm{desc}}\big(A_\phi(H_N),\ M^\dagger_{\Pi_{\mathrm{obs}}}\big)\ \le\ \eta_N+\rho,\qquad \eta_N\ =\ C\,\sqrt{\frac{d^*\ln(N_{\min}+1)+\ln(1/\delta)}{N_{\min}}}$$
($C$ the absolute constant of the classical VC uniform-deviation inequality, applied per relevant fiber and maximized). Moreover, under (A1)–(A2) the strong uniform Glivenko–Cantelli property holds: $d_{\mathrm{desc}}\big(A_\phi(H_N),M^\dagger\big)\to0$ almost surely as $N_{\min}\to\infty$ — **consistency of $A_\phi$ in the operator metric**, toward the *identified* operator $M^\dagger$ (never toward any single lift's $M^\star$: the identification width of LC-12/13 is not closed by sampling and does not appear in this bound).
*Proof.* Each endpoint of each constraint at each index is an empirical mean of one member of $\mathcal F^*$ over the relevant fiber; the VC inequality bounds $\sup_{\mathcal F^*}|\text{empirical}-\text{population}|$ by $\eta_N$ with probability $\ge1-\delta$; the sup over $\mathcal I$ and $\mathcal E$ is a sup over $\mathcal F^*$, so the bound transfers verbatim; transport adds $\rho$ per DE-T2. Almost-sure convergence is the VC strong-GC theorem. Measurability of the suprema: $\mathcal F^*$ is countable (countable atlas over countable $\mathcal I$), so all suprema are over countable families. $\square$

**Fallback LC-15′ (no complexity declaration). [proved]** Without (A2), the weighted union bound $\delta_j=\delta 2^{-j}$ over an enumeration of the atlas still yields *simultaneous coverage at every index* with index-dependent radii $\eta_{N,j}=\sqrt{\ln(4\cdot2^{j+1}/\delta)/2N_c}$ — valid uniform-in-index **coverage**, but no uniform **rate** ($\sup_j$ need not vanish). The gap between LC-15 and LC-15′ is exactly what the complexity declaration purchases; the contract permits either, tagged.

## 3. Transfer to the decision layer

**Theorem LC-16 (from description deviation to decision quantities). [proved / conditional as marked]**
(i) *Risk brackets:* for bounded loss with $\|\ell\|_\infty\le\bar L$, the LP risk brackets computed from the estimated class differ from those of $M^\dagger$'s class by at most the constraint perturbation transferred through the LP — bounded by $\bar L$ times the total-variation displacement of the feasible set; **[conditional]** under a declared interiority (Hoffman/condition-number) constant $\rho_0>0$ for the target polytopes, the feasible-set displacement is $\le C(\rho_0)\,(\eta_N+\rho)$, giving explicit end-to-end risk-bracket perturbation $\le\bar L\,C(\rho_0)(\eta_N+\rho)$. Without the interiority declaration, one-sided validity still holds (outer constraints widened by $\eta_N+\rho$ remain outer — conservative reporting needs no condition number).
(ii) *Robustness decisions:* the per-law LP robustness criterion (ML-X3) evaluated on the $(\eta_N{+}\rho)$-widened class is conservative: robustness certified there implies robustness for $M^\dagger$'s class. **[proved]** $\square$

**LC-17 (what the theorem does and does not claim). [declared]** LC-15 is a statement about the canonical estimator in the operator metric — the mathematically closed learnability of the *interface object*. It does not assert anything about any parameterized approximation family beyond measurability (LC-10); consistency of a *constrained* approximator within a declared subfamily of $\mathbb M$ is an instantiation-time question whose *target and metric are now fixed by this file* — which is exactly the boundary between mathematics and implementation that the mandate draws.
