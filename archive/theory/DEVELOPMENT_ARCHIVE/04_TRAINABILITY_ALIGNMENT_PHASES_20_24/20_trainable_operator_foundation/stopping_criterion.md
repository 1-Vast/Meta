# Stopping Criterion

> **Status:** Phase-20 terminal decision, 2026-08-03. Sources: the six deliverable files (TF-1–TF-20) and the audit `../20_final_theory_audit/BLOCKING_ISSUES.md` (minimal obstruction: $F_\omega$ introduced with $\Omega$, $\mathcal H$, the empirical objective undefined, and C3 imposed rather than derived). Phases 0–19.5 unmodified; no architecture, network name, biological example, application, or implementation choice appears. Stop condition: *an independent auditor can no longer claim "$F_\omega$ is undefined."*

---

## The obstruction, item by item, now removed

| Audit item | Resolution |
|---|---|
| **Parameter space $\Omega$ undefined** | $\Omega=\Xi\times\mathbb B^m$, $\Xi\subset\mathbb R^D$ nonempty compact — finite-dimensional, compact, explicit (TF-1); the realization map $G$ and its regularity (G1)–(G3) are the sole architecture interface, and $F_\xi=\pi_C\circ G$ is $C$-valued by the free $1$-Lipschitz projection (TF-2/3) |
| **Hypothesis class $\mathcal H$ undefined** | $\mathcal H=\{F_\omega\mid\omega\in\Omega\}$ with each $F_\omega$ a fully-specified input$\to$(valid operator $+$ certificate) map (TF-4); admissibility — measurability, continuity, support permutation invariance, set-valued-output compatibility, operator-metric consistency — all proved (TF-5), with an explicit class Lipschitz constant $L_{\mathcal H}$ (TF-6) |
| **Empirical objective undefined** | Population $R(\omega)$ and empirical $\widehat R_N(\omega)$ defined over observable tasks with identified targets; minimizers exist (Weierstrass on compact $\Omega$, no convexity needed) and admit a measurable selection (TF-7); the honesty firewall (training moves advice, not certificates) proved (TF-8) |
| **No theorem connects optimization to operator approximation; C3 imposed** | C3 **derived**, not imposed: a specified admissible witness family (piecewise-multilinear, TF-9) achieves $\sup_z\|F_{\xi^\star}-g^\star\|\le\varepsilon$ with explicit $D(\varepsilon)$ (TF-10); the imitation objective closes the empirical-to-operator arrow in $L^1(\mu_Z)$ unconditionally and in $\sup$ under a declared design floor (TF-14); the full chain $d_{\mathbb M}(F_{\hat\omega_N},A^\star)\le C_{\mathrm{stab}}\varepsilon+2\Gamma_N+\gamma^{\mathrm{opt}}$ is assembled from proved terms (TF-15); the fixed-capacity floor (TF-11) and all-resolution impossibility under uniform regularity (TF-12, the audit-corrected hypothesis) delimit it honestly |

## What remains outside the theorems (named, not hidden)

- **Optimization efficiency:** attainability of $\gamma^{\mathrm{opt}}$ is in-principle by compactness; *speed* is the implementer's — a scalar tolerance now, no longer an undefined bridge.
- **Capacity choice:** $D$ vs accuracy is a declared trade (TF-11/12); any $D$ is valid, higher $D$ buys accuracy, and the floor is displayed.
- **Design declaration:** the $\sup$-metric upgrade (TF-14(ii)) needs a declared mass floor $q_0$ on $Z$; without it, $L^1(\mu_Z)$ (average-case) approximation is what is claimed.
- **Carried scope limits (Phase-19.5):** skeleton-relative deployment; calibration from point-identified supervision only; ranking separately supervised. Unchanged.

None of these is the definition of $F_\omega$, $\Omega$, $\mathcal H$, or the objective — each of which is now a written mathematical object with proved properties. They are declared trades and tagged conditions, echoed in the contract.

## The interface an engineer instantiates

Choose any $D$ and any construction $G:\Xi\times Z\to\mathbb R^{\dim C}$ satisfying (G1)–(G3) — **any architecture whatsoever**, since (G1)–(G3) are its entire mathematical contract. Compose with the fixed $\pi_C$, the fixed convex assembly $\mathsf{asm}$, the fixed band-to-class $K$, support restriction, and the canonical side/certificate channels. Minimize $\widehat R_N$ (task-supervised) or $\widehat R^{\mathrm{im}}_N$ (imitation, oracle-free) over the compact $\Omega$. The result is a valid operator at every parameter (by type), generalizing at $\dim\Omega$-classical rates (TF-13), approximating the target to the specified-witness floor plus generalization plus optimization tolerance (TF-15) — with the capacity floor and scope limits stated in advance.

## Verdict

$$\boxed{\textbf{TRAINABLE\_OPERATOR\_FOUNDATION\_COMPLETE}}$$

The minimal mathematical completion is delivered: "$F_\omega$" now denotes a specific member of a defined hypothesis class over a defined compact parameter space, trained by a defined empirical objective, with a proved chain from empirical optimization to operator-metric approximation (C3 derived for a specified witness family, not assumed), and its limits proved and displayed rather than hidden. An independent auditor inspecting this phase cannot claim "$F_\omega$ is undefined": every symbol in that expression is a written object with a proof of its stated properties, and the one expression the previous audit called an unproved obligation (C3) is now a corollary of a construction. The only statements not reduced to theorems are explicitly named engineering trades — dimension, optimization speed, design resolution — none of which is a definitional gap.
