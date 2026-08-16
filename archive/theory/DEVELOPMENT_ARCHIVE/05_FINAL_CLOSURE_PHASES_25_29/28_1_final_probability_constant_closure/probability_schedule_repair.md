# Probability Schedule Repair (Item 2)

> **Status:** Phase-28.1, 2026-08-03. Closes the audit's second blocker: Phase-27.1 (MR-5) wrongly claimed $\log(1/\delta_N)/N\to0$ implies $\delta_N\to0$ and dropped the latter. The audit's witness $\delta_N=1/2$ refutes the implication. Both conditions are retained explicitly. New results **CL-5–CL-7**, tagged **[proved] / [declared]**.

---

## 1. The two conditions are independent — the implication is retracted

**Retraction CL-5. [retracted]** The Phase-27.1 statement "$\log(1/\delta_N)/N\to0$ implies $\delta_N\to0$" is **withdrawn**; it is false. Witness (the audit's): $\delta_N=\tfrac12$ gives $\log(1/\delta_N)/N=\log2/N\to0$, yet $\delta_N\not\to0$ and the confidence level $1-\delta_N=\tfrac12\not\to1$. The two conditions are logically independent:
- $\delta_N\to0$ — makes the **confidence level** $1-\delta_N\to1$ (the probability the bound holds tends to one);
- $\log(1/\delta_N)/N\to0$ — makes the **confidence contribution to $\Gamma_N$** ($\sqrt{\log(1/\delta_N)/N}$) vanish (the bound's *value* tends to its floor).
Neither implies the other ($\delta_N=\tfrac12$: second holds, first fails; $\delta_N=e^{-N}$: first holds, second fails). **Both are required.**

## 2. The corrected high-probability schedule

**Declaration CL-6 (high-probability consistency schedule — both conditions). [declared]**
$$\boxed{\ \text{Require BOTH}\quad \delta_N\to0\quad\text{AND}\quad \frac{\log(1/\delta_N)}{N}\to0,\quad\text{together with}\quad \frac{D_N\log(\Lambda N)}{N}\to0,\ \ \mathrm{mesh}(r_N)\to0,\ \ \gamma^{\mathrm{opt}}_N\to0.\ }$$
It is stated explicitly that **both** $\delta_N\to0$ and $\log(1/\delta_N)/N\to0$ are needed for high-probability consistency: the first for the confidence level, the second for $\Gamma_N\to0$. A concrete admissible choice satisfying both: $\delta_N=1/N$ ($\delta_N\to0$ and $\log N/N\to0$).

**Theorem CL-7 (high-probability consistency, correctly conditioned). [proved]**
Under (S-IID), CL-6's schedule, and $L_p^\star$ (CL-1): $\Gamma_N=C_0(\bar L+\tfrac\mu2)\sqrt{(D_N\log(\Lambda N)+\log(1/\delta_N))/N}\to0$ (both under-root terms vanish, using $\log(1/\delta_N)/N\to0$), and with probability $\ge1-\delta_N$,
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(2\Gamma_N+\gamma^{\mathrm{opt}}_N+L_p^\star\,\varepsilon_{\mathrm{approx}}(N)\big)+2h.$$
Because $\delta_N\to0$, the confidence level $1-\delta_N\to1$, so this is a genuine high-probability-tending-to-one statement; because $\log(1/\delta_N)/N\to0$ (and the other schedule terms vanish), the $\Phi$-argument $\to0$ and the bound tends to its floor $2h$. Both roles are covered by the two retained conditions. $\square$

**Declaration CL-8 (almost-sure clause — unchanged, correct). [declared]**
For an almost-sure eventual $\limsup_N\|d_{\mathbb M}(F_{\hat\omega_N},g^\star_\mu)\|_{L^2(\mu_\zeta)}\le2h$, require **$\sum_N\delta_N<\infty$** (Borel–Cantelli; this supplies $\delta_N\to0$ automatically) **together with $\log(1/\delta_N)/N\to0$** (for $\Gamma_N\to0$). Summability alone is insufficient ($\delta_N=e^{-N}$ is summable but violates the rate); the logarithmic-rate condition is independently required. A choice meeting both: $\delta_N=N^{-2}$. This clause is retained from Phase-27.1 unchanged; the audit found it correct. $\square$
