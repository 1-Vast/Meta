# Meta-Learning Theory: Existence, Identification, Estimation (Part IV)

> **Status:** Phase-9, 2026-08-03. Phases 0–7 frozen and cited. New results carry **ML-L** numbers, tagged **[proved] / [conditional] / [declared]**. This file supplies the separation the audit found incomplete: existence of the meta-operator, identification of the meta-operator, and finite-task estimation are three different statements with three different assumption stacks.

---

## 0. Typing warning (mandated): task exchangeability $\ne$ within-task IID

All probabilistic structure in this program lives **across tasks**. Within a task, the frozen noise model is adversarial with bounded support — no within-task stochastic assumption exists anywhere, and none is introduced here (introducing one would be a likelihood declaration, priced separately: DE-H4, ML-K2). "The tasks are exchangeable" is a statement about the sequence $(T_1,\dots,T_N,T_*)$ in $\mathcal T$; it implies nothing about the noise inside any record. Conflating the two levels is a type error, flagged here once and checked in the echo row.

---

## 1. Existence of the meta-operator

**Theorem ML-L1 (existence — population level, assumption-minimal). [proved]**
Let $\Pi$ be any law on $\mathcal T$ (tasks-with-latent-members; $\mathcal T$ standard Borel, ML-T1). Then the **ideal meta-operator**
$$M^*_\Pi(c,\gamma)\ =\ \text{the rung-appropriate conditioned decision-information object of }\Pi\ \text{at }(\kappa=c,\ \gamma)$$
exists: regular conditional laws exist on standard Borel spaces, and each value is a well-defined element of $\mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung}$ (singleton class at the top rung). Under task-**exchangeability** of an infinite task sequence, de Finetti gives a directing random measure $\Pi$; $M^*_\Pi$ exists conditionally on it, so the existence tier needs only exchangeability (or less: any single $\Pi$). Existence consumes **no data and no rates**. $\square$

## 2. Identification of the meta-operator

**Theorem ML-L2 (identification — partial, by inheritance from the frozen theory). [proved]**
What infinite history reveals is the law $\Pi_{\mathrm{rec}}$ of observable **records**, not of members. Since each member is seen only through its finite noisy record, $\Pi$ is not determined by $\Pi_{\mathrm{rec}}$; the object identified in the $N\to\infty$ limit is the **class**
$$\mathcal M^*(\Pi_{\mathrm{rec}})\ =\ \{M^*_\Pi:\ \Pi\ \text{compatible with}\ \Pi_{\mathrm{rec}}\},$$
whose $(c,\gamma)$-values are exactly the zero-sampling-width forced/compatible constraint polytopes: for each event $E$, $P(E)\in[\Pi_{\mathrm{rec}}(\text{record forces }E\mid c),\ \Pi_{\mathrm{rec}}(\text{record compatible with }E\mid c)]$. **The meta-operator is point-identified iff the per-task censoring width vanishes $\Pi_{\mathrm{rec}}$-a.s.** (each task's own record decides the decision-relevant functional — the frozen per-task gate holding almost surely). Partial identification of the meta-operator is thus a *theorem-level identification statement*, not an estimation artifact: it is the frozen per-task partial identification, integrated over the population. This is the exact content of "second-order partial identification" (DE-R5), now placed at its correct tier. $\square$

## 3. Finite-task estimation

**Theorem ML-L3 (estimation — rates under the declared cross-task stack). [conditional]**
Under task-**(IID)** — or **(C-IID-$\kappa$)** with the fiber count $N_c$ from the **multiset** (ML-T3/T4) — and the declared finite event family $\mathcal E$: the empirical forced/compatible intervals with margins $\eta_{N_c}=\sqrt{\ln(4|\mathcal E|/\delta)/2N_c}$ cover, simultaneously over $\mathcal E$, the identified-class values of ML-L2, with probability $\ge1-\delta$; a declared transport class adds its radius $\rho$. Total interval width $=$ **identification width** (ML-L2; $N$-irreducible, reduced only by richer per-task designs) $+\ 2\eta_{N_c}$ (reduced by tasks in the fiber) $+\ 2\rho$ (reduced by tighter transport). Bare task-exchangeability yields **no** rate (8.1 shared-Bernoulli witness, carried); it yields the mixture-conditional reformulation only. Zero fibers: the vacuous fallback (ML-T4). $\square$

## 4. The three-tier separation, stated once

| Tier | Object | Assumptions | Data role |
|---|---|---|---|
| **Existence** | $M^*_\Pi$ | a task law (or exchangeability + de Finetti); standard Borel $\mathcal T$ | none |
| **Identification** | the class $\mathcal M^*(\Pi_{\mathrm{rec}})$; point iff a.s. zero censoring | the frozen per-task information model; nothing statistical | infinite-history limit |
| **Estimation** | $\widehat{\mathcal Q}_{c,\gamma}$ with $\eta_{N_c}$ margins | (IID) or (C-IID-$\kappa$) [+ transport]; finite $\mathcal E$ or declared uniform bound | finite $N$, multiset-typed |

**Corollary ML-L4 (stop-condition 2). [proved]** Finite-history learning (tier 3) is separated from existence (tier 1) by an intermediate, non-trivial identification tier (tier 2) with its own iff-condition; no statement of one tier is used at another tier's type. In particular: the meta-operator can exist and be unidentified (censored populations); be identified and not yet estimated (small $N_c$); and estimation error vanishing does **not** close the identification width — the three failure modes are distinct, and the ledger reports them distinctly. $\square$
