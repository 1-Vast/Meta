# Zero-Fiber Finite-Sample Repair (§2)

> **Status:** Phase-13, 2026-08-03. Repairs the audit's finite-sample defect: the finite-$N$ bound took a supremum over all contexts while exempting zero fibers, and the rung alignment silently cancelled a coordinate that differs on unseen positive-mass contexts. New results **PM-4–PM-5**, tagged **[proved] / [conditional]**.

---

## 1. The missing-fiber event, priced

**Definition.** Relevant contexts: $C^+=\{c\in C_\kappa:\Pi_{\mathrm{obs}}(\kappa=c)>0$ and $r^{\mathrm{decl}}$ assigns rung $\ge2$ somewhere at $c\}$. Missing-fiber event: $\mathrm{Miss}_N=\{\exists c\in C^+:N_c=0\}$.

**Lemma PM-4 (probability of the missing-fiber event). [proved]**
Under task-(IID) with context marginal $\pi(c)=\Pi_{\mathrm{obs}}(\kappa=c)$:
$$\Pr(\mathrm{Miss}_N)\ \le\ \sum_{c\in C^+}(1-\pi(c))^N\ \le\ |C^+|\,e^{-N\pi_{\min}},\qquad \pi_{\min}=\min_{c\in C^+}\pi(c)>0,$$
(union bound; $C^+$ finite). $\pi_{\min}$ is a population constant; if it is to appear numerically in a certificate, it must be declared or lower-bounded from data with its own tagged confidence — otherwise only the exact-sum form is emitted. On $\mathrm{Miss}_N$, the estimator's rung coordinate at the missing context is $1$ while the target's is $\ge2$: $d_R=1$ there — this discrepancy is **real and is not cancelled** (the Phase-11 cancellation is retracted, PM-2c(ii)). $\square$

## 2. The repaired finite-$N$ consistency theorem

**Theorem PM-5 (finite-$N$ bound with the zero-fiber term; a.s. consistency). [conditional on the declared stack (S1)–(S6) of Phase 11, with PM-1's metric]**
(i) *(finite $N$)* With probability at least $1-\delta_N-\Pr(\mathrm{Miss}_N)$:
$$d_{\mathbb M}\big(\widehat M_N,M^\dagger\big)\ \le\ \alpha\,\tfrac12\bar H\,\big(\eta_N+\rho\big)\ +\ \beta\,\delta_N,$$
and **unconditionally** (all sample paths):
$$\mathbb E\ d_{\mathbb M}\big(\widehat M_N,M^\dagger\big)\ \le\ \alpha\,\tfrac12\bar H\,(\eta_N+\rho)+\beta\,\delta_N\ +\ (\alpha+\beta+\gamma)\big(\delta_N+\Pr(\mathrm{Miss}_N)\big),$$
using the trivial bound $d_{\mathbb M}\le\alpha+\beta+\gamma$ on the bad events. **No uniform finite-$N$ claim is made without the $\Pr(\mathrm{Miss}_N)$ term** — the mandated repair, and the honest reading: before every relevant fiber has been seen, the operator is provably at rung-distance $\gamma$ from its target at the unseen contexts, and the bound says so.
(ii) *(almost-sure consistency, $\rho=0$)* $C^+$ finite with positive masses ⇒ by the strong law of large numbers $N_c/N\to\pi(c)>0$ a.s. for each $c\in C^+$, so every relevant fiber is observed from some finite time on, a.s., and $\mathbf 1\{\mathrm{Miss}_N\}\to0$ a.s.; combined with the strong uniform GC ($\eta$-part) and the schedule ($\delta_N\to0$): $d_{\mathbb M}(\widehat M_N,M^\dagger)\to0$ a.s. — the asymptotic statement the audit judged supportable, now derived *through* the explicit event rather than around it. $\square$

*Bookkeeping note.* The estimator's own ledger surfaces $\mathrm{Miss}$: at any unseen relevant context it emits the vacuous rung-1 value **flagged "fiber unobserved"** — so the finite-sample discrepancy the bound charges is also visible per-emission, not only in the aggregate constant.
