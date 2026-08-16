# Conditional Population Learnability (Closure Target 4)

> **Status:** Phase-8.2, 2026-08-03. Repairs the audit's finding 4: DR-L3-R estimates **population** event frequencies, but $V_M$/DR-M1-R consumed them as a law **conditioned on the current observations $O$** — a step no theorem established (a marginal law of $g(f)$ does not determine $P(g(f)\mid O)$). New results **DC-C1–C4**, tagged **[proved] / [conditional] / [impossible]**. The mandate offers routes A (joint-distribution learning), B (conditional fiber construction), C (declared sufficient statistic). **Chosen route: B, with the sufficiency declaration made explicit** (which is where C's content lives); A is stated as the declared alternative. All assumptions explicit.

---

## 1. The gap, stated precisely

**DC-C1 (what was missing). [proved]** Let $g$ be the decision-relevant map. History under (C-IID) identifies (interval-)information about the population marginal $\pi_g=\mathrm{law}(g(f))$ (or fiberwise $\pi_g(\cdot\mid c)$). The decision layer needs $P(g(f_\beta)\in\cdot\mid O)$. These differ whenever $g(f)$ is dependent on the observable record: witness — tasks with $f(x_1)\in\{0,1\}$ equiprobable and $g(f)=f(x_q)=f(x_1)$; the marginal is $\mathrm{Ber}(\tfrac12)$, but $O$ with $\tilde y_1\approx1$ makes the conditional degenerate at $1$. Using the marginal where the conditional belongs is a real typing error, not a formality. Two legitimate mechanisms exist for $O$ to act, and they are different in kind: (i) **support restriction** — $P(f_\beta\in I(O)\mid O)=1$ under the frozen bounded-noise semantics, likelihood-free (DE-H2); (ii) **density reweighting** *within* the support — which requires statistical machinery beyond the frozen model (DE-H4). The repair must supply (ii) honestly. $\square$

---

## 2. Route B: declared conditional fiber construction

**Declared objects and assumptions.**
- **($\kappa$-DESIGN)** A declared measurable **context map** $\kappa$ from task records to a finite set $C_\kappa$, computable from any task's *observable* record alone — historical: $\kappa(\text{record}_i)$; current: $c_\beta=\kappa(O)$. ($\kappa$ may read the design, the observed values, declared auxiliary labels $c_i$ — never unidentified quantities.)
- **(C-IID-$\kappa$)** Conditional on $\kappa=c$, the task records-and-members $(\mathrm{rec}_i,f_i)$ are IID, and the current task is a further independent draw from the same conditional population (transport class of radius $\rho$ declared if not).
- **(SUFF-$\kappa$) [a declaration, not a derivation]** $\kappa$ is **decision-sufficient**: in the population model, $g(f)\perp \mathrm{record}\,\mid\,\kappa(\mathrm{record})$ — given the fiber, the residual record carries no further population information about $g(f)$.

**Theorem DC-C2 (the conditioning theorem). [conditional on the declared assumptions]**
Under ($\kappa$-DESIGN) + (C-IID-$\kappa$) + (SUFF-$\kappa$), for every admissible noise kernel $\lambda\in\Lambda_{\mathrm{adm}}$ (frozen support semantics):
$$P_\lambda\big(g(f_\beta)\in\cdot\mid O\big)\;=\;\underbrace{\pi_g\big(\cdot\mid \kappa(O)\big)}_{\text{fiber conditional}}\ \big|\ \text{restricted to the identified support}\ g\big(I(O)\big),$$
in the honest set-valued sense: every $\lambda$-conditional is supported in $g(I(O))$ and, on events measurable w.r.t. $\kappa$-sufficient structure, agrees with the fiber conditional; consequently the DR-L3-R/DC-R5 fiber estimate — computed from the $n_{c_\beta}$ historical tasks in the current fiber, with simultaneous event intervals at $\eta_{n_{c_\beta}}$ and the support constraint $P(\widehat\Sigma)=1$ adjoined — is a valid $1-\delta$ confidence class **for the current-observation-conditioned object**, because the conditioning statistic $\kappa(O)$ is a function of $O$ and the only further $O$-dependence permitted by (SUFF-$\kappa$) is the support restriction, which is likelihood-free.
*Proof sketch (assumption-tracking is the content).* (C-IID-$\kappa$) makes the fiber sample a draw from $\pi(\cdot\mid c_\beta)$; (SUFF-$\kappa$) collapses $P(g\in\cdot\mid \mathrm{record})$ to $P(g\in\cdot\mid\kappa(\mathrm{record}))$ at the population level; the frozen support argument (DE-H2, verbatim) forces every admissible posterior inside $I(O)$; DR-L3-R supplies the finite-sample class within the fiber. Each step consumes exactly one declared assumption; none is hidden. $\square$

**Theorem DC-C3 (the conditioning ladder — what each declaration buys, and the impossibility floor). [proved / impossible]**
Without (SUFF-$\kappa$), the gap between the $\kappa$-conditional and the full-$O$-conditional is **not estimable** from the declared data: by DE-H4 (one prior, many posteriors), even exact knowledge of the joint population law leaves the $O$-conditional undetermined within the frozen noise model — the set of admissible conditionals shares only the support restriction. Hence the exhaustive ladder, each rung a strictly stronger declaration:
| Declared | Valid conditioned object | Decision consequence |
|---|---|---|
| nothing (population layer absent/failed) | all laws on $g(I(O))$ | $\Gamma$-minimax collapses to the frozen minimax endpoint (DE-T4) — graceful, honest |
| (C-IID-$\kappa$) alone | fiber class **as a marginal-typed object**; may inform only decisions declared insensitive to residual-record dependence | must be tagged marginal; consuming it as conditional is forbidden |
| + (SUFF-$\kappa$) | fiber class **as the conditional** (DC-C2) | full Phase-7/8 conditional machinery legitimate |
| + declared likelihood (LIK) | single posterior | Bayes tier |
No rung's output may be consumed at a higher rung's type — this typing rule *is* the repair of DR-M1-R's over-claim. $\square$

**Route A (declared alternative, stated for completeness). [conditional]** Learn the **joint** pushforward $\mathrm{law}(s(f),g(f))$ for a declared observable-statistic map $s$ evaluable on the current design (requires the strong per-task gate: each historical task's own data (interval-)identify $s(f_i)$ *and* $g(f_i)$ — an archive-coverage condition at the current $(D,Q)$), then condition on the observed $s$-cell with a declared positive-mass bound $q_0$ (rates $\eta_{n}/q_0$-degraded). Route A is Route B with $\kappa=$ a declared discretization of $s$ — the fiber construction is the common normal form, which is why B is the chosen route.

---

## 3. What this closes

**DC-C4 (repaired $V_M$ semantics). [proved given DC-C2/C3]** $M_\phi$'s output is re-typed as a pair (class, **rung tag**). $V_M$ now reads: the class covers, with probability $\ge1-\delta$, the true object *of its tagged rung* — marginal fiber law at rung 2, conditional law at rung 3, posterior at rung 4 — and $D_\psi$ may consume it only at that rung (type check, DC-C3). The composition theorem can now legitimately say "conditioned": at rung 3 the conditioning is proved (DC-C2); at lower rungs the operator provably degrades to the conservative endpoint instead of over-claiming. The audit's finding — a marginal used as a conditional — is structurally impossible in the re-typed interface. $\square$
