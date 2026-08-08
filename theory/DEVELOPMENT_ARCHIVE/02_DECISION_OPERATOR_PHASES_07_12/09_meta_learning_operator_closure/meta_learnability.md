# Meta-Learning Learnability (§5)

> **Status:** Phase-9 closure, 2026-08-03. New results **MC-16–MC-19**, tagged **[proved] / [conditional] / [declared]**. The three questions — existence, identification (uniqueness), finite-history learning — are answered at three different types with three different assumption stacks; no statement of one tier is used at another tier's type.

---

## 0. Typing warning (mandated, enforced once)

**Within-task IID does not exist in this program.** The within-task noise model is the frozen adversarial bounded-support model; introducing a within-task stochastic assumption is a likelihood declaration priced at rung 4 (MC-14). All IID/exchangeability statements below are **task-level**: about the sequence $(T_1,\dots,T_N,T_*)$ in $\mathbb T$. Conflating the levels is a type error checked in the echo row.

---

## 1. Existence: does $A_\phi$ exist?

**Theorem MC-16 (existence of the target, and of the operator). [proved]**
(i) *Ideal target.* For any task law $\Pi$ on the standard Borel $\mathbb T$ (or an exchangeable law on $\mathbb T^\infty$ via its directing measure), the ideal object $M^*_\Pi\in\mathbb M$ — $(c,\gamma)\mapsto$ the rung-appropriate conditional decision-information object of $\Pi$ — exists by existence of regular conditional laws. No data, no rates, no uniqueness claim.
(ii) *Operator.* $A_\phi$ exists constructively for every finite $H_N$: the canonical forced/compatible interval-polytope construction (MC-8) is a total, assumption-tagged map $\bigcup_N\mathbb T^N\to\mathbb M$, defined even at empty history and empty fibers (vacuous rung-1 values). Existence is thus unconditional at both levels. $\square$

## 2. Identification: is $A_\phi$ uniquely determined?

**Theorem MC-17 (identification is partial, and exactly quantifiable). [proved]**
Uniqueness must be asked of the **target**, and the answer has three layers:
(i) *Given $\Pi$:* $M^*_\Pi$ is unique (conditional laws are a.s.-unique) — the ideal object is well-posed.
(ii) *Given the observable record law $\Pi_{\mathrm{rec}}$ (the $N\to\infty$ limit of what history reveals):* $\Pi$ is **not** determined — members are seen only through finite noisy records — and the identified object is the class $\{M^*_\Pi:\Pi\ \text{compatible with}\ \Pi_{\mathrm{rec}}\}$, whose $(c,\gamma)$-values are exactly the zero-sampling-width forced/compatible polytopes: $P(E)\in[\Pi_{\mathrm{rec}}(\text{record forces }E\mid c),\ \Pi_{\mathrm{rec}}(\text{record compatible with }E\mid c)]$. **Point identification holds iff the per-task censoring width vanishes $\Pi_{\mathrm{rec}}$-a.s.** — the frozen per-task partial identification, integrated over the population; this is the theorem-level home of "second-order partial identification".
(iii) *The estimator:* within the identified class, the canonical interval-polytope choice is a **declared canonical selection** (it is the class itself, reported honestly), so no spurious uniqueness is asserted where none exists. $\square$

## 3. Learning: can finite $H_N$ estimate $A_\phi$?

**Theorem MC-18 (finite-history estimation). [conditional on the declared stack]**
Assumptions, each declared: **(task-IID)** — or task-exchangeability *demoted to its true strength* (mixture reformulation only; **no rates**: the shared-Bernoulli witness is carried — all-task-identical copies of one Bernoulli are exchangeable and never concentrate), or **(C-IID-$\kappa$)** with fiber counts; **(finite task complexity)** — a declared finite event family $\mathcal E$ (union-bound constant $\ln(4|\mathcal E|/\delta)$) or a declared uniform-convergence condition with its constants for infinite families; **(concentration)** — Hoeffding for bounded event indicators (or a separately declared inequality with constants for declared dependence). Then, simultaneously over $\mathcal E$, with probability $\ge1-\delta$, the estimated class covers the tier-(ii) identified values within $\eta_{N_c}=\sqrt{\ln(4|\mathcal E|/\delta)/2N_c}$ (+ declared transport radius $\rho$). Total width $=$ **identification width** (tier (ii); $N$-irreducible, reduced only by richer per-task designs) $+\,2\eta_{N_c}$ (reduced by tasks in the fiber) $+\,2\rho$ (reduced by tighter transport). $\square$

## 4. The separation, certified

**Corollary MC-19 (stop-condition 2). [proved]**
The three tiers are logically independent in both directions along the ladder: the target can exist and be unidentified (censored populations); be identified and unestimated (small $N_c$); and estimation error $\to0$ never closes the identification width. Each tier's failure has its own observable signature (vacuous class vs wide-but-stable interval vs shrinking interval) and its own remedy (declare structure vs improve per-task designs vs collect fiber tasks) — reported distinctly in the ledger. Finite-history learning is thereby separated from existence by a non-trivial identification tier with an iff-condition, and no theorem crosses tiers. $\square$
