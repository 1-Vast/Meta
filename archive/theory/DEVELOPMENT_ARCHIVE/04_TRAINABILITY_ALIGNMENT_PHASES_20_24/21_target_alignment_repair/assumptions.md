# Assumptions (Route B)

> **Status:** Phase-21, 2026-08-03. Every assumption used to define, approximate, and calibrate the single target $g^\star$ — declared explicitly, none tacit. The audit's prohibition ("assuming continuity without declaring it") is honored: **continuity of $g^\star$ is not assumed here — it is derived in PT-6 from the assumptions below.** New labels **A-***, tagged **[declared]**. Each carries what it buys and what fails without it (the latter detailed in `failure_modes.md`).

---

## The declared assumption set

- **(A-SKEL) [declared].** A fixed finite deployment skeleton: compact metric statistic domain $Z$ (finite union of compact cubes); finite anchor set with $\operatorname{conv}(\{b^{\mathrm{pop}}\}\cup\{b_j\})=\mathbb B$; declared value grid/mesh (Route-B continuous scalar scope, Phase-19.5). *Buys:* the learnable object is a map $Z\to C$ with no loss of expressiveness (PT-1 §1).

- **(A-STAT) [declared].** Task sampling is (IID) — or (C-IID-$\kappa$) with fiber counts — from the observable task law $\Pi_{\mathrm{obs}}$; a declared transport radius $\rho$ if the deployment population differs ($\rho=0$ otherwise). Within-task noise remains the frozen adversarial-bounded model (no within-task stochastic assumption). *Buys:* the local risk $\ell_0(z,\cdot)$ is a well-defined conditional expectation and the empirical risk concentrates (generalization bridge).

- **(A-LOSS) [declared].** The band-score loss $L(\cdot,a)$ is bounded by $\bar L$, convex and $L_{\mathrm{Lip}}$-Lipschitz in the band argument, uniformly in the observable target $a$; on the point-supervised channel it is the elicitable interval score (Phase-19.5 DT-6). *Buys:* convexity and continuity of $\ell_0$ in $c$; calibration compatibility.

- **(A-SC) [declared].** $\ell(z,\cdot)$ is $\mu$-strongly convex in $c$ on $C$, $\mu>0$: either $L\circ\mathsf{asm}$ is already strongly convex in $c$ (then $\mu$ is its modulus, declared) or a ridge $\tfrac\mu2\|c\|^2$ is adjoined (then $g^\star=g^\star_\mu$ is the regularized target, PT-3 ownership). *Buys:* uniqueness and — crucially — Lipschitz continuity of $z\mapsto g^\star(z)$ (PT-6) and the quadratic risk-to-coefficient inequality underlying calibration (PT-8).

- **(A-CONT) [declared].** The conditional loss field is continuous in the statistic: $z\mapsto\ell_0(z,c)$ is continuous uniformly in $c$, with a declared modulus $\varpi_\ell$ (equivalently: $z\mapsto\mathrm{law}(A_T\mid\zeta=z)$ is weakly continuous and $L$ is bounded-Lipschitz, which delivers $\varpi_\ell$). *Buys:* continuity of the risk field, from which continuity of the argmin $g^\star$ follows by strong convexity (PT-6). **This is the honest replacement for Phase-20's assumed target continuity: we declare regularity of the *risk field* (a statement about the data-generating law and the loss, checkable in principle) and *prove* the target map's continuity, rather than assuming the latter directly.**

- **(A-DESIGN) [declared, optional].** For sup-metric (rather than $L^1(\mu_\zeta)$) approximation claims: a design measure on $Z$ with a mass floor $q_0>0$ on the skeleton mesh. *Buys:* upgrade of average-case to worst-case coefficient control (generalization bridge); without it, only $L^1(\mu_\zeta)$ is claimed.

## Dependency map (which theorem needs which)

| Result | Assumptions |
|---|---|
| PT-2 target well-defined (unique) | A-STAT, A-LOSS, A-SC |
| PT-3 Bayes optimality | A-LOSS, A-SC |
| PT-6 target continuity (**derived**) | A-SKEL, A-LOSS, A-SC, A-CONT |
| PT-9 approximation $\inf_F d(F,g^\star)\to0$ | + A-SKEL (witness family) |
| PT-8/PT-10 calibration inequality | A-LOSS, A-SC (strong convexity is the load-bearing one) |
| PT-11 empirical→population→$d_{\mathbb M}$ chain | all of the above (+ A-DESIGN for sup form) |

No assumption beyond this list is used anywhere in Phase 21. In particular, **no continuity of $g^\star$ is assumed** (it is PT-6), and **no interpolation-existence argument is used as a learnability proof** (the witness family establishes approximability; learnability is the calibration + generalization chain, PT-10/PT-11) — the two prohibitions the audit named.
