# Failure Modes (Route B)

> **Status:** Phase-21, 2026-08-03. Adversarial pass over the single-target construction: what breaks each assumption, and confirmation that no failure re-introduces target-switching or a hidden second target. Results **PT-12–PT-16**, tagged **[attack refuted] / [found — priced] / [found — scoped]**.

---

## PT-12 — the regularizer bias (strong convexity secured by ridge)

**Attack.** (A-SC) via an adjoined ridge makes $g^\star=g^\star_\mu$ a *different* object from the unregularized risk-optimal $g^\star_0$ — a covert target switch.
**Outcome: [found — priced, single target preserved].** $g^\star_\mu$ is the unique, explicitly declared target of the deployment (PT-1 ownership); it is *one* object, not two. The gap to the unregularized minimizer is bounded: if $g^\star_0$ is unique and $\ell_0(z,\cdot)$ is $L_c$-Lipschitz, $\|g^\star_\mu(z)-g^\star_0(z)\|\le\sqrt{2L_c\,\mathrm{diam}(C)/\mu}$-type control, and $g^\star_\mu\to g^\star_0$ pointwise as $\mu\downarrow0$ (standard Tikhonov/Moreau argument). The $\mu$-trade is explicit: small $\mu$ → small bias but large calibration constant $C_{\mathrm{stab}}\sqrt{2/\mu}$ (PT-10); large $\mu$ → the reverse. The deployment declares $\mu$; the target is unambiguous given it. If the base score is *natively* strongly convex, $\mu$ is its modulus and there is no bias at all. $\square$

## PT-13 — continuity failure (A-CONT violated)

**Attack.** The conditional loss field $z\mapsto\ell_0(z,\cdot)$ is discontinuous (a task-law discontinuity across a statistic boundary); then PT-6 fails and $g^\star$ may jump, defeating PT-9.
**Outcome: [found — scoped honestly].** PT-6 requires (A-CONT), declared and (in principle) checkable as weak-continuity of the task-conditional law plus bounded-Lipschitz loss. Where it holds, $g^\star$ is Lipschitz by proof. Where a deployment cannot assert it, the phase makes **no** approximation claim — rather than assuming continuity of $g^\star$ directly (the Phase-20 error). This is the correct behavior: the failure is declared as an unmet hypothesis, not hidden. Partial remedy within scope: on each closed cube stratum of $Z$ (A-SKEL) continuity need only hold stratum-wise; cross-stratum jumps are absorbed by the finite stratification (the witness family reads the stratum label), so only *within-stratum* discontinuity is fatal — a genuinely weaker requirement. $\square$

## PT-14 — the two-target relapse (the audit's core issue)

**Attack.** Does any theorem secretly reintroduce the canonical operator $A^\star$ as a target, or use its observability where $g^\star$'s is needed?
**Outcome: [attack refuted].** Audit of every use: the objective (PT-7/PT-11) is task risk against observable identified targets $A_T$ — $A^\star$ absent; the target (PT-1) is $\arg\min_c\ell$ — defined by that same risk; approximation (PT-9) targets $g^\star$ via its *proved* modulus (PT-6), not a computed $A^\star$; calibration (PT-10) controls $d_{\mathbb M}(F,g^\star)$ by excess *risk*, needing no evaluation of any oracle. The fixed $b^{\mathrm{pop}}$ is an *input constant to the assembly*, not a target. **No $A^\star$ evaluation, no imitation label, no second target appears anywhere.** The Phase-20 relapse channel (metric imitation of a canonical operator) is structurally absent because the objective is risk, not imitation. $\square$

## PT-15 — interpolation-as-learnability relapse

**Attack.** Is PT-9 (approximation witness) doing the work that only calibration+generalization should?
**Outcome: [attack refuted].** PT-9 supplies exactly one term ($\varepsilon_{\mathrm{approx}}$) of PT-11; the passage from empirical minimization to operator error runs through excess *population risk* and the calibration inequality $\Phi$ (PT-10), with generalization $\Gamma_N$ carrying the finite-sample step. Interpolation existence is never equated with learnability — the prohibition is honored, and removing PT-9 would only remove the approximation term, not collapse the chain. $\square$

## PT-16 — design-dependence of the sup upgrade

**Attack.** The sup-metric calibration hides a design assumption.
**Outcome: [found — scoped].** Exactly (A-DESIGN): the sup upgrade needs a declared mesh mass floor $q_0$; without it, the honest claim is $L^2(\mu_\zeta)$ (average-case) operator convergence. Stated, not upgraded silently. Average-case convergence to the risk-optimal operator is itself a complete Route-B result; sup is an optional strengthening at declared cost. $\square$

---

## Summary

| Failure | Result |
|---|---|
| ridge bias | priced; single declared target ($g^\star_\mu$), gap bounded, $\mu$-trade explicit |
| continuity failure | scoped; (A-CONT) declared, only within-stratum discontinuity fatal, no silent assumption |
| two-target relapse | refuted; $A^\star$ absent from target, objective, approximation, calibration |
| interpolation-as-learnability | refuted; PT-9 is one term, calibration+generalization carry learnability |
| sup design-dependence | scoped; (A-DESIGN) or honest $L^2(\mu_\zeta)$ claim |

No failure reintroduces a second target or an undeclared assumption. The construction has exactly one target under exactly one declared assumption set.
