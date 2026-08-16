# Task Space and Historical Information (§1)

> **Status:** Phase-9 closure, 2026-08-03. Phases 0–7 frozen; repaired decision-theory statements (Phases 8.1/8.2 as corrected, `../09_meta_learning_operator_formalization/`) cited, not modified. Audit of record: `../11_final_closure_audit/` (`DECISION_OPERATOR_INVALID`); its archive-type defect is closed here at the type level. New results carry **MC-** numbers, tagged **[proved] / [conditional] / [declared]**.

---

## 1. The task space

**Definition MC-1 (task space $\mathbb T$).** A task is a tuple
$$T=(\,O,\ S,\ Q,\ \gamma\,)\in\mathbb T,$$
- **observation object** $O$: the complete observable record — design $D$, values $\tilde y$, noise level $\varepsilon$, optional auxiliary label (CI-A);
- **support information** $S\subseteq O$: the sub-record released for adaptation;
- **query object** $Q$: finite query points plus the demanded pushforward map $g$ (values / differences / order type — the DR-J quotient lattice); $\bot$ permitted for historical tasks;
- **decision specification** $\gamma$: action set with declared abstention cost, loss, criterion, tie-break, tolerance, confidence; $\bot$ permitted for historical tasks.

$\mathbb T$ is declared standard Borel (regular conditional laws exist — consumed once, in §5). Behind each task stands an unobserved member $f_T\in\mathcal F$; only $O$ is observable.

## 2. Historical meta-training information

**Definition MC-2.** $H_N=(T_1,\dots,T_N)\in\mathbb T^N$ — an **ordered sequence**, with the quotient lattice
$$\mathbb T^N\twoheadrightarrow\mathrm{mult}(H_N)\ \text{(multiset)}\twoheadrightarrow\mathrm{set}(H_N)\ \text{(distinct records)}.$$

**Theorem MC-3 (why $H_N$ cannot be a simple set). [proved]**
(i) The meta-learning layer's estimation theorems consume empirical frequencies $\tfrac1N\#\{i:\cdot\}$ and radii $\eta_N=\sqrt{\ln(4|\mathcal E|/\delta)/2N}$ — functionals of the **multiset**: deduplication alters frequencies, effective sample size, and every confidence constant. Under any record law with atoms, independent tasks collide with positive probability, so a duplicate record is *evidence of population mass*, not redundancy. A set-typed $H_N$ makes the learning theorems inapplicable to their own input — the audited defect, now excluded by the type of $H_N$ itself.
(ii) Order: under declared task-exchangeability the sample law is order-invariant and $\mathrm{mult}(H_N)$ is a sufficient quotient [proved by invariance]; the sequence is retained as primitive because declared drift/transport structure is order-indexed. $\square$

**Theorem MC-4 (identification invariance vs meta-learning statistical information — the separation). [proved]**
The two channels consume different quotients, and must:
- **Identification invariance.** Archive feasibility information (windows/traces; frozen F17, Phase 6) is a property of *which* traces occur — the frozen channel-separation theorem (DE-H1) makes it invariant under permutation **and duplication**. Hence the feasibility channel legally consumes $\mathrm{set}(H_N)$, and *must* (multiplicity influence on the identified set would let frequencies shrink it — forbidden by DE-H2/H3).
- **Meta-learning statistical information.** Population frequencies are exactly what multiplicities carry; the frequency channel consumes $\mathrm{mult}(H_N)$ (the sequence, under declared drift).
**Typing rule:** $\mathrm{set}(\cdot)$ is applied only on the feasibility path; multiset counts only on the frequency path. Either crossing is a type violation — one fabricates identification from frequencies, the other destroys the statistical sample. This single rule is the entire legal interaction between $H_N$'s algebraic structure and the two channels, and it discharges stop-condition 3. $\square$

**Convention MC-5 (fiber counts and the zero-fiber fallback). [declared]** For the declared context map $\kappa$ (§4): $N_c=\#\{i:\kappa(O_i)=c\}$, a multiset count. $N_{\kappa(O_*)}=0$ entitles the population layer to no claim at the current context: its output is the vacuous rung-1 object (all laws on the identified support — the frozen minimax endpoint), never an undefined value and never a silently borrowed neighbor fiber.
