# Joint Ranking Object Repair (Closure Target 3)

> **Status:** Phase-8.2, 2026-08-03. Repairs the audit's finding 3: the 8.1 confidence object was mis-typed (intervals over a queried subset $S\subsetneq S_m$ read as "a law on $S$"; overlapping pairwise events treated as outcomes). New results **DC-R1–R5**, tagged **[proved] / [classical, scoped] / [conditional]**. The identification-level order object $\Sigma(J)$ (DR-J3/J4) is unchanged; this file fixes the **population/probability layer** above it.

---

## 1. The correctly typed object

**Definition DC-R1 (outcome space and constraint class).**
For $m$ queries, the outcome space is the **full** permutation set
$$\Omega_m=S_m\qquad(\text{ties handled by a declared convention: either a declared tie-breaking refinement, or }\Omega_m^\pm=\text{the set of weak orders — declared once}).$$
Population ranking information is a **law $P\in\Delta(\Omega_m)$** — or, honestly under finite data, a **constraint class**
$$\widehat{\mathcal Q}\;=\;\Big\{P\in\Delta(\Omega_m):\ \ell(E)\le P(E)\le u(E)\ \ \forall E\in\mathcal E\Big\},$$
where $\mathcal E$ is the declared finite family of **events** $E\subseteq\Omega_m$ — single orders $\{\pi\}$, pairwise events $E_{ab}=\{\pi:\pi\text{ ranks }a\text{ before }b\}$, top-$k$ events — entering as *constraints on one law over the full space*, never as outcomes of a law on a subset. Un-queried mass is unconstrained (beyond $\sum_\pi P(\pi)=1$), which is exactly the honest reading of partial information: the audit's pathology ($S=\{123\}$ forcing $P(123)=1$) is gone by construction — the interval $[\ell,u]$ for the single queried event constrains, and mass elsewhere remains free. $\widehat{\mathcal Q}$ is a polytope: convex, compact, closed under the LP operations below. **Pairwise probabilities are marginals:** $p_{ab}=P(E_{ab})$ — projections of $P$, carrying strictly less information than $P$ (DC-R4).

**Support ceiling (identification layer, carried).** When the current-task identified order object $\widehat\Sigma\supseteq\Sigma(J_Q(O))$ is imposed, the constraint $P(\widehat\Sigma)=1$ joins the class. Outerness of $\widehat\Sigma$ makes this a valid relaxation of the true ceiling $P(\Sigma)=1$ (DE-H2 lifted to orders); positive-mass failure ($\widehat{\mathcal Q}$ empties) is the DC-A4(3) failure flag.

---

## 2. Pairwise consistency

**Theorem DC-R2 (realizability of a pairwise matrix as marginals). [classical, scoped]**
A matrix $(p_{ab})$ arises as $\big(P(E_{ab})\big)$ for some $P\in\Delta(S_m)$ iff it lies in the **linear ordering polytope** $P^m_{LO}$ (convex hull of the $0/1$ order indicators). Always necessary: $p_{ab}+p_{ba}=1$ and the $3$-dicycle inequalities $p_{ab}+p_{bc}-p_{ac}\le 1$. These are **sufficient for $m\le5$** (Grötschel–Jünger–Reinelt); for $m\ge6$ the polytope has further facets (fence/Möbius-ladder classes) and the dicycle conditions no longer suffice. Consequences: (i) checking pairwise consistency is finite LP feasibility in all cases; (ii) at the program's budget-relevant sizes ($m\le5$) the explicit inequalities characterize it; (iii) an *interval-valued* pairwise object is jointly realizable iff the constraint polytope $\widehat{\mathcal Q}$ (with those interval constraints) is nonempty — again LP feasibility, and emptiness is a declared-information inconsistency → failure flag, not renormalization. $\square$

---

## 3. Listwise decisions from $P$ — what marginals do and do not decide

**Theorem DC-R3 (pairwise-decomposable losses: marginals are Bayes-sufficient). [proved]**
For Kendall discordance $K(\sigma,\pi)=\sum_{a<b}\mathbf 1\{\sigma,\pi\text{ disagree on }(a,b)\}$ (and any loss that is an affine combination of pairwise disagreement indicators):
$$\mathbb E_{P}\,K(\sigma,\pi)\;=\;\sum_{(a,b):\,\sigma\text{ ranks }a\text{ before }b}\big(1-p_{ab}\big),$$
linear in the pairwise marginal matrix. Hence Bayes (and $\Gamma$-minimax over marginal-defined classes) listwise decisions under such losses are functions of $(p_{ab})$ alone. (The *argmin over $\sigma$* is a linear-ordering optimization; its mathematical well-posedness is all that is claimed here.)

**Theorem DC-R4 (non-decomposable losses: marginals are insufficient — equal-marginal witness). [proved]**
Let $P_1=\tfrac12\delta_{123}+\tfrac12\delta_{321}$ and $P_2=$ uniform on $S_3$. Both have every pairwise marginal $p_{ab}=\tfrac12$. Under exact-match $0$–$1$ loss ($\ell(\sigma,\pi)=\mathbf 1\{\sigma\ne\pi\}$): Bayes risk under $P_1$ is $\tfrac12$ with optimal set $\{123,321\}$; under $P_2$ it is $\tfrac56$ with all six orders tied. Identical marginals, different optimal actions and different achievable risks: exact-match, top-$k$-type, and other non-decomposable listwise losses require the law on $\Omega_m$ (or its relevant non-pairwise events in $\mathcal E$), not the pairwise matrix. $\square$

---

## 4. Confidence bounds, correctly composed

**Theorem DC-R5 (from event intervals to decision bounds). [conditional on the 8.1 rungs]**
Under (IID)/(C-IID) with the per-task forced/compatible bounds and the union-bound allocation over the declared finite $\mathcal E$ (DR-L3-R, unchanged as an *event-interval* theorem), with probability $\ge1-\delta$ the true law lies in $\widehat{\mathcal Q}$ of DC-R1 — now a correctly-typed object. Then for any listwise action $\sigma$ and loss $\ell$:
$$\Big[\min_{P\in\widehat{\mathcal Q}}\mathbb E_P\,\ell(\sigma,\cdot),\ \max_{P\in\widehat{\mathcal Q}}\mathbb E_P\,\ell(\sigma,\cdot)\Big]$$
are attained linear programs over the polytope (objective linear in $P$), giving valid simultaneous risk brackets for all $\sigma\in S_m$; the robust ($\Gamma$-minimax) listwise action minimizes the upper end, and **decision robustness** (Tier 2, DE-R6 lifted listwise) holds iff one $\sigma$'s upper end is below all others' lower ends — an LP-checkable condition. Every quantity lives on $\Delta(\Omega_m)$; no subset-law or event-as-outcome typing appears anywhere in the chain. $\square$

---

$$\boxed{\begin{array}{c}\Omega_m=S_m;\ \text{population ranking information}=\text{a constraint polytope of laws on the full space, events as constraints;}\\ \text{pairwise consistency}=\text{linear-ordering-polytope membership (dicycle-complete for }m\le5\text{); Kendall-type decisions need only marginals (DC-R3);}\\ \text{exact-match/top-}k\text{ need the law (DC-R4 witness); all confidence and robustness statements are LPs over the polytope (DC-R5).}\end{array}}$$
