# Bayes Optimality (Item 3, rewritten to $g^\star_\mu$)

> **Status:** Phase-23.1, 2026-08-03. The optimality statement, rewritten so it refers only to $g^\star_\mu$ and to the *regularized* risk $R_\mu$ — never to an unregularized "same risk" claim (the audit's invalid transfer). New result **TI-5**, tagged **[proved]**.

---

**Theorem TI-5 (regularized Bayes optimality — one target). [proved]**
Define, for measurable $F:Z\to\Delta_m$, the regularized population risk
$$R_\mu(F)\ =\ \mathbb E_\zeta\big[\,L_0(\zeta,B\,F(\zeta))+\tfrac\mu2\|F(\zeta)\|^2\,\big]\ =\ \mathbb E_\zeta\big[J_\mu(\zeta,F(\zeta))\big]$$
(pointwise separability by the tower property over $\zeta$; $F(\zeta)$ is $\zeta$-measurable). Then
$$R_\mu(F)\ \ge\ \mathbb E_\zeta\big[\min_{p\in\Delta_m}J_\mu(\zeta,p)\big]\ =\ \mathbb E_\zeta\big[J_\mu(\zeta,g^\star_\mu(\zeta))\big]\ =\ R_\mu(g^\star_\mu),$$
with equality iff $F=g^\star_\mu$ $\mu_\zeta$-a.e. Hence **$g^\star_\mu$ is the unique (a.e.) minimizer of the regularized risk $R_\mu$ over all measurable maps** — the risk-optimal target *for the objective that defines it*. $\square$

**Scope of the claim (Item 4, enforced here). [declared]**
This optimality is with respect to $R_\mu$, the regularized risk. No statement is made that $R_\mu(g^\star_\mu)$ equals, or that $g^\star_\mu$ minimizes, the unregularized risk $R_0(F)=\mathbb E_\zeta[L_0(\zeta,BF(\zeta))]$. The Phase-22.1 phrase "Bayes optimality transfers with the same risk" is retracted (TI-3); it conflated $R_0$ and $R_\mu$ and implicitly identified two coordinate-dependent targets. Here the risk and the target are one matched pair $(R_\mu,g^\star_\mu)$, and the optimality theorem is stated only for that pair.

**Consequence for the assembled band. [proved]**
Because $g^\star_\mu$ is the unique argmin, the assembled target band is $\beta^\star_\mu(z)=B\,g^\star_\mu(z)\in\mathbb B$, a single well-defined valid band per $z$; the operator value is $K(\beta^\star_\mu(z))$ restricted to $\mathrm{supp}\,I(S)$, with the $\omega$-invariant certificate/confidence/rung channels unchanged. One target coefficient map ⇒ one target operator. $\square$
