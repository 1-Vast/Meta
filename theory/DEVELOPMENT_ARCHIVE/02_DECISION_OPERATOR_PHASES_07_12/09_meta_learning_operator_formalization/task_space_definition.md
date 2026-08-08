# Task Space Definition (Part I)

> **Status:** Phase-9 (meta-learning operator formalization), 2026-08-03. Phases 0–7 frozen and cited. Audit of record: `../11_final_closure_audit/FINAL_VERDICT.md` (verdict `DECISION_OPERATOR_INVALID`; the archive-type defect is repaired *here*, the remaining defects in the companion files). New results carry **ML-T** numbers, tagged **[proved] / [conditional] / [declared]**.

---

## 1. Tasks

**Definition ML-T1 (task space).** Fix the frozen ambient structure ($\mathcal X$, family assumptions as declared). The task space $\mathcal T$ consists of tuples
$$T=(\,O,\ S,\ Q,\ \gamma\,)$$
- **observation object** $O\in\mathcal O=(\mathcal X\times\mathbb R)^{\le k}\times[0,\infty)\times C^{\le1}$: the complete observable record — design, values, noise level $\varepsilon$, optional auxiliary label (CI-A);
- **support object** $S\subseteq O$: the sub-record released for adaptation (in the frozen setting $S$ is the $(x_i,\tilde y_i)$ list; the distinction matters when parts of $O$ are metadata);
- **query object** $Q\in\mathcal X^{\le m}\times\mathcal G_{\mathrm{push}}$: the query points and the demanded pushforward map $g$ (values, differences, order type — DR-J lattice); optional for historical tasks ($\bot$ allowed);
- **decision specification** $\gamma\in\mathrm{Ctx}$: action set (incl. declared abstention cost), loss, criterion, tie-break, tolerances, confidence (8.2 $\mathrm{Ctx}$); optional for historical tasks ($\bot$ allowed).
$\mathcal T$ is declared standard Borel (needed once, for regular conditional laws in `meta_learning_theory.md`). Behind each task stands an unobserved member $f_T\in\mathcal F$; only $O$ is ever seen.

## 2. The meta-training sample is a sequence, never a set

**Definition ML-T2.** The historical meta-training sample is the **ordered sequence**
$$H_N=(T_1,\dots,T_N)\ \in\ \mathcal T^N,$$
with its quotient lattice
$$\mathcal T^N\ \twoheadrightarrow\ \mathrm{mult}(H_N)\ (\text{multiset})\ \twoheadrightarrow\ \mathrm{set}(H_N)\ (\text{set of distinct records}).$$

**Theorem ML-T3 (why a set is forbidden — the audit's defect, closed at the type level). [proved]**
(i) *Frequencies are the learning target.* The population learning theorems (DR-L3-R and successors) consume empirical frequencies $\tfrac1N\#\{i:\cdot\}$ and confidence radii $\eta_N=\sqrt{\ln(4|\mathcal E|/\delta)/2N}$: both are functionals of the **multiset** — deduplication changes the frequencies, the effective sample size, and the radii. Since two independent tasks produce identical records with positive probability under any atomic record law, an exact duplicate is *evidence of population mass, not redundancy*: mapping $H_N\mapsto\mathrm{set}(H_N)$ before the frequency channel destroys precisely the information the channel estimates, and no theorem downstream survives (this is the Phase-8.2 interface error: $\mathcal H$ was typed as a finite set for both operators).
(ii) *Order is the carrier of declared non-exchangeable structure.* Under declared task exchangeability the sample law is order-invariant and the multiset is a **sufficient quotient** (order adds nothing — proved by invariance); under a declared drift/transport structure indexed by arrival order, the sequence itself is required. Hence the primitive is the sequence; the multiset is its canonical quotient *under EXCH*, and the set is a valid quotient for **exactly one** consumer:
(iii) *The feasibility channel, and only it, uses the set.* Identification-level archive information (windows/traces, frozen F17 and Phase-6) is a property of *which* traces occurred, not how often: the frozen channel-separation theorem (DE-H1) makes the feasibility channel duplication-invariant **by theorem**. So the 8.2 axiom $V_I$(iii) (set semantics for $I_\theta$) is retained — but as the *channel projection* $\mathrm{set}(\cdot)$ applied inside a sequence-typed interface, not as the type of $H_N$ itself. $\square$

**Typing rule ML-T4 (no frequency destruction — stop condition 3). [declared, enforced]**
$$I_\theta\ \text{consumes}\ \mathrm{set}(H_N);\qquad A_\phi/M_\phi\ \text{consume}\ \mathrm{mult}(H_N)\ \big(\text{the full sequence under declared drift}\big).$$
Applying $\mathrm{set}(\cdot)$ anywhere on the path into the frequency channel is a type violation. Conversely, letting multiplicities influence $I_\theta$ would let frequencies shrink the identified set — forbidden by the frozen DE-H2/H3. The two projections are the *entire* legal interaction between the sample's algebraic structure and the two channels.

## 3. Fiber bookkeeping (used throughout Phase 9)

For a declared context map $\kappa$ (`conditional_population_repair.md`): $N_c=\#\{i\le N:\kappa(O_i)=c\}$ — a multiset count. The **zero-fiber fallback** is declared here once: $N_{\kappa(O_*)}=0$ entitles the population layer to *no* claim at the current context — its output degrades to the vacuous rung-1 object (all laws on the identified support; frozen minimax endpoint), never to an undefined or silently-borrowed quantity.
