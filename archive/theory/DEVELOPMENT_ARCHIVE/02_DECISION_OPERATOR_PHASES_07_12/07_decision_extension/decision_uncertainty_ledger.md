# Decision Uncertainty Ledger (Parts X–XI)

> **Status:** Phase-7, 2026-08-03. Frozen corpus cited, not modified. New results carry **DE-U** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Two mandates: (X) a decision-error decomposition with a which-term-shrinks audit; (XI) proof that the decision extension's uncertainty quantities are **distinct, mutually non-determining objects** that must be reported side by side.

---

## 1. The four ledger quantities

For a decision problem at $(O,Q,\mathcal A,L)$ with declared information $\Delta=(\mathcal Q,\lambda\text{-declaration},\text{criterion})$:

| Symbol | Quantity | Depends on | Epistemic grade |
|---|---|---|---|
| $\rho_{\mathrm{id}}$ | identification radius $\tfrac12\operatorname{diam}$ of the relevant section (frozen; per query or per functional, e.g. $\tfrac12\operatorname{diam}\Delta_{ab}$) | $O$, family, $\varepsilon$ | unconditional (Theorem 1 certificate) |
| $\sigma_{\mathrm{prob}}$ | probabilistic spread of the conditioned weighting (posterior quantile width / variance) | $O$ **and** $(\pi,\lambda)$ | conditional on (EXCH/IID)+(LIK) |
| $r_{\mathrm{dec}}(a)$ | decision regret of the selected action over the declared class: $\sup_{\mu\in\mathcal Q}\big[\int L(a)\,d\mu-\inf_b\int L(b)\,d\mu\big]$ | $O$, $\Delta$, $L$ | conditional on the class declaration |
| $u_{\mathrm{rank}}$ | ranking tier (DE-R6) + the $p_{ab}$ interval | $O$, $\Delta$, history | mixed: tier 1 unconditional; tiers 2–3 conditional |

---

## 2. Mutual non-determination (mandate XI): proofs and counterexamples

**DE-U1 ($\rho_{\mathrm{id}}$ does not determine $\sigma_{\mathrm{prob}}$). [proved]** Same identified interval $[0,1]$; weighting uniform vs weighting concentrated near $0$: posterior spreads differ arbitrarily. The radius is weighting-free by construction. $\square$

**DE-U2 ($\sigma_{\mathrm{prob}}$ does not bound $\rho_{\mathrm{id}}$ — the central honesty exhibit). [proved]** Query off the covered set: the frozen section is all of $\mathbb R$, $\rho_{\mathrm{id}}=+\infty$ (F18/C1). A declared prior (e.g. Gaussian pushforward on the query value) nevertheless yields a proper posterior with **finite, arbitrarily small** $\sigma_{\mathrm{prob}}$. Every bit of that concentration is manufactured by the declaration, none by observation. A report showing $\sigma_{\mathrm{prob}}$ without $\rho_{\mathrm{id}}$ is indistinguishable from evidence — which is precisely what the ledger exists to prevent. $\square$

**DE-U3 ($\rho_{\mathrm{id}}=0\Rightarrow$ everything collapses; the converse fails). [proved]** If the relevant section is a singleton, all weightings on it coincide, every criterion selects the same action, $r_{\mathrm{dec}}=0$, $\sigma_{\mathrm{prob}}=0$: identification determines decision. Conversely $r_{\mathrm{dec}}=0$ does **not** imply $\rho_{\mathrm{id}}=0$: on a symmetric section with a reflection-closed declared class $\mathcal Q$ and symmetric loss, every $\mu\in\mathcal Q$-Bayes action is the center (DE-S4 symmetry argument applied inside $\mathcal Q$), so regret over the class vanishes while the radius is positive. Decision consensus is not identification. $\square$

**DE-U4 (regret does not determine the radius, quantitatively). [proved]** Scale the symmetric example: on the section $[-M,M]$ with a reflection-closed class and symmetric loss, the selected center has $r_{\mathrm{dec}}=0$ while $\rho_{\mathrm{id}}=M$ is unbounded. Oppositely, tiny sections can carry maximal regret under tie-tier weightings: two-point section $\{\pm\delta\}$ with $0$–$1$ ranking loss gives regret $\tfrac12$ for any strict choice under the symmetric class, radius only $\delta$. Neither direction transfers. $\square$

**DE-U5 (ranking uncertainty is not a function of marginal radii). [proved]** DE-J6/DE-R2 witnesses: identical marginal radii $\tfrac12$ with $u_{\mathrm{rank}}$ ranging from Tier-1-identified-tie (diagonal) to fully ambiguous (square/anti-diagonal); and sign-identified problems with huge value radii (diagonal scaled by $M$), versus sign-ambiguous problems with radii $\delta$ ($\{(\delta,-\delta),(-\delta,\delta)\}$). Value uncertainty and ordering uncertainty are independent coordinates of the ledger. $\square$

**DE-U6 (the certificate is criterion-free; the criterion adds no certificate). [proved]** $\rho_{\mathrm{id}}$ is invariant under every element of $\Delta$ (it is computed before $\Delta$ enters); and no choice of criterion produces an unconditional guarantee below it (DE-L5(i)). Hence the four quantities are pairwise non-redundant, and the ledger is minimal: deleting any row loses information the others provably cannot reconstruct (DE-U1–U5). $\square$

---

## 3. Decision-error decomposition (mandate X)

Setting: implemented action $\hat a$ obtained by $\gamma$-approximately minimizing the estimated robust risk $\hat\rho$ built from: history of $n$ members (possibly censored, DE-L4), declared transport class of radius $\rho_{TV}$ (DE-T2(ii)), declared likelihood, criterion $\Gamma$-minimax over the resulting class. Benchmark: the same criterion under the true conditioned population law. Loss bounded, $\bar V$-BV in the scalar decision functional (DE-L3(iii) regime).

**Theorem DE-U7 (five-term decomposition). [conditional on the declared axioms; each term proved under them]**
$$\underbrace{\rho^{\mathrm{true}}(\hat a)-\inf_a\rho^{\mathrm{true}}(a)}_{\text{criterion-relative excess risk}}\ \le\ \underbrace{\gamma}_{\text{(v) action/optimization}}\;+\;2\Big[\underbrace{\bar V\,\eta_n}_{\text{(iii) finite-sample}}\;+\;\underbrace{\bar V\,\tfrac12(\hat p^+-\hat p^-)\ \text{[width]}}_{\text{(ii) decision-object estimation: censoring}}\;+\;\underbrace{2\rho_{TV}\|L\|_\infty}_{\text{(iv) distribution shift}}\Big],$$
by the chain $\rho^{\mathrm{true}}\to\rho^{\pi_h\text{-exact}}\to\rho^{\text{empirical}}\to\hat\rho$ (triangle inequality through the three intermediate risks; each gap bounded by DE-T2(ii), DE-L4, DE-L3(iii) respectively) plus the standard $2\sup|\hat\rho-\rho|+\gamma$ argument. Standing **outside** this chain, not summable with it:
$$\text{(i) identification ambiguity: the unconditional worst-case error of *any* action remains }\ \ge\rho_{\mathrm{id}}\ \text{(Theorem 1)},$$
which no term above touches — the excess-risk chain is *relative to the criterion*, while (i) is *absolute*. $\square$

**Corollary DE-U8 (which terms more historical members reduce). [proved]**
- (iii) finite-sample $\eta_n=O(\sqrt{\ln(1/\delta)/n})$: **reduced by $n$** — the only one.
- (ii) censoring width: reduced only by **per-member coverage** (richer historical designs), never by $n$ (DE-L4).
- (iv) shift term: reduced only by **tighter declared transport** (or actually re-observing the new population), never by $n$ from the old population (DE-T3 in the limit).
- (v) optimization $\gamma$: computational, orthogonal.
- (i) identification ambiguity: reduced by **nothing in the decision layer** — only by new observations of the current member or stronger declared family structure, i.e. by the frozen theory's own levers. $\square$

The ledger row order is deliberately the audit order: first the unconditional floor, then what assumptions bought, then what the selected action risks, then the ordering tier. A Phase-7 output is complete iff all four rows are present with their tags (`final_decision_operator.md`, axiom H6).
