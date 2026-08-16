# Meta-Learning Objective (Deliverable 4)

> **Status:** Phase-20, 2026-08-03. Defines the empirical and population objectives over $\Omega$ and proves minimizer existence and measurable selection. New results **TF-7–TF-8**, tagged **[proved] / [conditional]**.

---

## 1. The objects

- **Task distribution:** $T\sim P(T)$, the observable task law $\Pi_{\mathrm{obs}}$ (marked lifts stay declared-model objects; the objective touches observables only).
- **Support / query / target:** $S_T$ (finite set, mandatory $\varepsilon$, optional declared label), $Q_T$ (query + pushforward), $A_T$ the task's **identified** target information — value / forced-compatible interval / admissible-order set (observable; on the point channel a value, on the censored channel an interval, per the Phase-19.5 two-channel scope).
- **Loss:** $L:\mathbb M\text{-value}\times\mathcal A\to[0,\infty)$, the declared operator-level score — convex and Lipschitz in the band argument, bounded by $\bar L$; the elicitable interval score on the point channel, the forced/compatible estimation contribution on the censored channel (calibration claimed only on the point channel — DT-6).

## 2. The objective and its solvability

**Definition (population and empirical risk).**
$$R(\omega)\ =\ \mathbb E_{T\sim P(T)}\big[\,L\big(F_\omega(S_T,Q_T;z_H),\,A_T\big)\,\big],\qquad \widehat R_N(\omega)\ =\ \tfrac1N\sum_{t=1}^N L\big(F_\omega(S_{T_t},Q_{T_t};z_H),\,A_{T_t}\big),$$
$$\omega^\star\ \in\ \operatorname*{arg\,min}_{\omega\in\Omega}\ R(\omega),\qquad \hat\omega_N\ \in\ \operatorname*{arg\,min}_{\omega\in\Omega}\ \widehat R_N(\omega).$$

**Theorem TF-7 (existence and measurable selection — no convexity needed). [proved]**
(i) *Continuity of the risk:* $\omega\mapsto L(F_\omega(\cdot),A_T)$ is continuous (TF-5.2 + $L$ Lipschitz in the band argument, which moves continuously with the class by the stability constant); bounded by $\bar L$; so $R$ and $\widehat R_N$ are continuous on the compact $\Omega$ (dominated convergence for $R$).
(ii) *Existence:* both minima are attained (Weierstrass on compact $\Omega$) — the minimizing sets are nonempty and compact. **No convexity is assumed** (the realization $G$ is generic); existence rests only on compactness of $\Omega$ and continuity, which the admissibility conditions (G1)–(G3) supply. The Phase-18 convexity claim is *not* invoked here — a deliberate weakening, since a generic architecture is not convex in its parameters, and the audit's convexity correction (DT-0.3) is thereby moot for this objective.
(iii) *Measurable selection:* $\hat\omega_N$ can be chosen as a measurable function of the sample (the argmin correspondence of a Carathéodory function on a compact metric space admits a measurable selector — Kuratowski–Ryll-Nardzewski), so $\hat\omega_N$ is a well-defined random element and $R(\hat\omega_N)$ a random variable — required for the generalization statement to be typed. $\square$

**Remark TF-8 (what the minimizer targets — carried honesty). [declared]** $R$'s minimizer over $\mathcal H$ elicits, on the point channel, the calibrated conditional band (classical interval-score elicitation); it does **not** target the outer identified band, which is carried by the $\omega$-invariant certificate row (DT-6, MR-13). Training moves advice; certificates are fixed. The objective is thus honest by type regardless of $\omega^\star$: a poorly-fit $\omega$ yields bad advice priced by $R$, never an invalid or falsely-certified output (TF-5.4). $\square$
