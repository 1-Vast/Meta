# Stopping Criterion (Target Alignment Repair)

> **Status:** Phase-21 terminal decision, 2026-08-03. Sources: the five deliverable files (PT-1–PT-16) and the audit `../21_trainable_operator_audit/` (`TRAINABLE_OPERATOR_FOUNDATION_INVALID`; sole obstruction: two targets, never proved identical). Phases 0–19.5 unmodified; Phase-20's two-target alternation abandoned. **Route B chosen; Route A not used.** No architecture, application, or implementation appears. Stop condition: *an auditor cannot claim "the learner has no single mathematical target."*

---

## The single target, and the four proofs that make it a foundation

$$\boxed{\ g^\star(z)=\operatorname*{arg\,min}_{c\in C}\ \ell(z,c)\quad\text{— the risk-optimal coefficient map, defined by the task risk alone (PT-1).}\ }$$

| Mandate item | Delivered | Status |
|---|---|---|
| **One target, no switching** | $g^\star$ defined solely by the declared task risk (PT-1); $A^\star$ absent from target, objective, approximation, and calibration (PT-14 audit) — the fixed $b^{\mathrm{pop}}$ is an assembly input constant, not a target | one target, proved |
| **1. Define $g^\star:Z\to C$** | PT-1; well-defined single-valued function by strong convexity (PT-2); Bayes-optimal over all measurable maps (PT-3) | proved |
| **2. Minimal explicit assumptions** | A-SKEL, A-STAT, A-LOSS, A-SC, A-CONT (+ optional A-DESIGN), each declared with what it buys and what fails without it (`assumptions.md`); continuity of $g^\star$ is **not** among them | declared, complete |
| **3. Approximation $\inf_F d(F,g^\star)\to0$** | Target continuity **derived** (PT-6: Lipschitz with modulus $\varpi_\ell/\mu$), then uniform approximation of that same $g^\star$ by the witness family, transferred to $d_{\mathbb M}$ (PT-9) | proved |
| **4. Calibration $\text{op error}\le\Phi(\text{excess risk})$, $\Phi\to0$** | $\|d_{\mathbb M}(F,g^\star)\|_{L^2(\mu_\zeta)}\le C_{\mathrm{stab}}\sqrt{2\,\mathcal E(F)/\mu}$ (PT-10), via pointwise separability (PT-7) and the strong-convexity lower bound (PT-8); sup form under A-DESIGN | proved |
| **5. Empirical loss → population risk → operator error valid** | PT-11 end-to-end: excess risk = generalization ($2\Gamma_N$, covering numbers over compact $\Omega$) + approximation ($L_c\varepsilon_{\mathrm{approx}}$, PT-9) + optimization ($\gamma^{\mathrm{opt}}$), fed through $\Phi$; consistency PT-11.1; all to the single $g^\star$ | proved / tagged conditional |

## The two forbidden moves, both avoided

- **Interpolation-existence as learnability:** PT-9 supplies only the approximation *term* of PT-11; learnability runs through calibration (PT-10) and generalization — excess *population risk*, not interpolation, controls operator error (PT-15).
- **Assuming continuity of the target:** not assumed — PT-6 *derives* $g^\star$'s Lipschitz continuity from declared regularity of the risk field (A-CONT) plus strong convexity (A-SC), which is precisely the gap the audit identified in Phase 20's multilinear argument.

## Residual, named (not definitional)

The $\mu$ of (A-SC) — native modulus or declared ridge — fixes *which* risk-optimal target (PT-12, bias bounded, $\mu$-trade explicit); (A-CONT) is a checkable hypothesis whose failure is declared, not hidden (PT-13, only within-stratum discontinuity fatal); the sup-metric upgrade needs (A-DESIGN), else the honest claim is $L^2(\mu_\zeta)$ (PT-16); optimization *speed* remains the implementer's. None of these is a second target or a tacit assumption.

## Verdict

$$\boxed{\textbf{TARGET\_ALIGNMENT\_COMPLETE}}$$

The learner has exactly one mathematical target, $g^\star$, defined by the declared task risk. Its continuity is proved (not assumed), its approximability is proved for that same object, its operator error is bounded by an explicitly vanishing function of excess population risk, and the empirical→population→operator chain converges to it. An auditor inspecting Phase 21 cannot claim "the learner has no single mathematical target": the target is one named object, the guarantees all refer to it, and the previously-conflated canonical operator plays no role in the objective. The one degree of freedom that remains — the strong-convexity modulus $\mu$ — selects *which* declared risk-optimal target, not *whether* there is a single one.

**Verdict: `TARGET_ALIGNMENT_COMPLETE`.**
