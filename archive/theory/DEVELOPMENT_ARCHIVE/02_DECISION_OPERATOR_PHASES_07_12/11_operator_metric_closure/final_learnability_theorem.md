# Final Learnability Theorem (§5)

> **Status:** Phase-11, 2026-08-03. The operator-level learnability theorem, now with a **proved** metric transfer — closing the audit's central gap. New results **OM-7–OM-9**, tagged **[proved] / [conditional]**. The three tiers remain separated; none is used to establish another.

---

## 1. The assumption stack (all declared, all echoed)

- **(S1)** Task-level (IID) from the declared task distribution $\Pi_{\mathrm{obs}}$ — or (C-IID-$\kappa$) with fiber counts $N_c$, minimum relevant count $N_{\min}\to\infty$; exchangeability admitted only as mixture reformulation (no rates).
- **(S2)** Route-A operator class (metric transfer): bounded outcome complexity $\bar n,\bar e$; hence uniform Hoffman constant $\bar H$ (OM-1/OM-2). *(Route B substitutes the declared stability inequality with its verification burden.)*
- **(S3)** Finite effective dimension: declared VC bound $d^*$ for the context-conditioned forcing/compatibility indicator class over the full index-and-event atlas (OM-6's enlarged class).
- **(S4)** Pullback-closed atlas discipline (OM-5) with the shared $\delta$-schedule.
- **(S5)** Confidence schedule $\delta_N\downarrow0$ with $\delta_N\ge e^{-o(N_{\min})}$ (e.g. $\delta_N=1/N$), so $\eta_N(\delta_N)\to0$ while the confidence coordinate converges (OM-5-pre).
- **(S6)** Transport radius $\rho$ declared if the current population differs ($\rho=0$ under (S1) to the current population).

## 2. The theorem

**Theorem OM-7 (finite-history operator learnability, in the true operator metric). [conditional on (S1)–(S6)]**
Let $\widehat M_N=A_\phi(H_N)$ (canonical estimator, coherent by OM-5, an element of the complete metric space $(\mathbb M,d_{\mathbb M})$ of OM-3′) and let $M^\dagger$ be the fully-typed identified target (OM-4). Then:
(i) *(finite-$N$ bound)* with probability $\ge1-\delta_N$:
$$d_{\mathbb M}\big(\widehat M_N,\ M^\dagger\big)\ \le\ \max\Big\{\ \tfrac12\bar H\,\big(\eta_N+\rho\big),\ \ \delta_N\ \Big\},\qquad \eta_N=C\sqrt{\frac{d^*\ln(N_{\min}+1)+\ln(|C_\kappa|/\delta_N)}{N_{\min}}};$$
(ii) *(consistency)* if $\rho=0$: $\ d_{\mathbb M}\big(\widehat M_N,M^\dagger\big)\longrightarrow0$ almost surely as $N_{\min}\to\infty$.
*Proof.* (i) OM-6 gives the uniform endpoint event ($d_{\mathrm{desc}}\le\eta_N+\rho$, all contexts/queries/specifications/events simultaneously); on that event OM-2 transfers it to the polytope coordinate: $\sup_\iota d_H^{TV}\le\tfrac12\bar H(\eta_N+\rho)$ (both polytopes nonempty: empirical and population witnesses); OM-5-pre assembles the three coordinates ($\text{rung gap}=0$; confidence gap $=\delta_N$). (ii) The strong uniform Glivenko–Cantelli property of the (S3) class gives $d_{\mathrm{desc}}\to0$ a.s.; the transfer constant $\bar H$ is deterministic and uniform; $\delta_N\to0$ by (S5). No step appeals to the retracted Phase-10 sentence: the passage from descriptions to the operator metric is OM-2, a proved theorem with the audited counterexample excluded by (S2) — and shown irreducible without it. $\square$

**Corollary OM-8 (what converges to what — the mandate's $d_{\mathbb M}(\widehat A_\phi(H_N),A_\phi)\to0$, typed honestly). [proved]**
The limit object is the **identified operator** $M^\dagger_{\Pi_{\mathrm{obs}}}$ — "the true operator" in the only sense the observable task distribution defines. Convergence to any single marked lift's ideal $M^\star_{\Pi^\bullet}$ is impossible in general (the lift class is non-degenerate whenever censoring persists, LC-3(iii)/LC-12), and the distance $d_{\mathbb M}(M^\dagger,M^\star_{\Pi^\bullet})$ — the identification width, now measured in the same metric — is a population constant that sampling provably does not reduce. $\square$

## 3. The three tiers, final form

| Tier | Statement | Status |
|---|---|---|
| **Existence** | $M^\star_{\Pi^\bullet}\in\mathbb M$ for every declared marked law (regular conditionals; all coordinates assigned as in OM-4 with confidence $1$); $\widehat M_N\in\mathbb M$ unconditionally (OM-5, totality) | proved; no data |
| **Identification** | $\Pi_{\mathrm{obs}}$ determines $M^\dagger$ exactly (a genuine $\mathbb M$-element, OM-4) and the lift class only up to it; point-latent-identification iff a.s. zero censoring; eventwise sharpness proved, joint sharpness not claimed (outer semantics is the contract) | proved / scoped |
| **Learning** | OM-7: finite-$N$ uniform bound and a.s. consistency in the **complete operator metric**, with proved metric transfer, coherent estimator, context-complete confidence accounting, and a convergent confidence schedule | conditional on (S1)–(S6), proved under them |

Existence is never cited for learning; learning converges to identification's object; identification's residual width is a floor, not an error term. **[proved discipline, verified per statement]**
