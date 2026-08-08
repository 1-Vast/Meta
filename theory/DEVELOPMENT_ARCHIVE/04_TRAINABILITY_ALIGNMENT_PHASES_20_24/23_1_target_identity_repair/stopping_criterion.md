# Stopping Criterion (Target Identity Repair)

> **Status:** Phase-23.1 terminal decision, 2026-08-03. Sources: the four repair files (TI-1–TI-11) and the audit `../23_final_theory_freeze_audit/FINAL_AUDIT.md` (`THEORY_STILL_INVALID`; minimal obstruction: the ridge in barycentric coordinates changed the regularized minimizer, so the package held two targets). Phases 0–21 unmodified; operator not redesigned; scope unchanged. Only target-identity repair performed.

---

## The obstruction, removed

The audit's singular obstruction: Phase 22.1 simultaneously (a) preserved the Phase-21 regularized target $g^{\mathrm{old}}_\mu$ (ridge on $c=(\lambda,w)$) and (b) asserted a new-coordinate target (ridge on $p$) with the "same risk" — and these differ (one-anchor, zero-loss witness: $\lambda=0$ vs $\lambda=\tfrac12$). The repair:

1. **One target defined (Item 1):** $g^\star_\mu(z)=\arg\min_{p\in\Delta_m}[\ell(z,Bp)+\tfrac\mu2\|p\|^2]$ — TI-1, well-defined, unique, everywhere-defined (TI-2).
2. **Equivalence claim removed (Item 2):** every statement implying $g^{\mathrm{old}}_\mu=g^\star_\mu$ is explicitly retracted (TI-3); the old target is not carried, not referenced, not claimed equal; the ridge non-invariance under the nonlinear reparameterization is stated as the reason.
3. **All theorems rewritten to $g^\star_\mu$ (Item 3):** target definition (TI-1/2), Bayes optimality against the *regularized* risk $R_\mu$ (TI-5, "same risk" purged), continuity (TI-6), approximation (TI-7), calibration with explicit floor (TI-8), consistency (TI-9) plus the typings the audit found missing — interpolation resolution tied to $\dim\Omega_N$ by node-count, $\delta_N\to0$ declared, and a defined common limiting space with a coherent target family for the mesh-refining case (TI-10/11).
4. **Regularized-target honesty (Item 4):** $g^\star_\mu$ is declared the regularized risk-optimal target and explicitly **not** claimed equal to the unregularized Bayes target $g^\star_0$ (TI-4, TI-5 scope).
5. **Scope kept (Item 5):** fixed deployment $z_H=z_H^0$, continuous point-valued affinity regression only, no ranking guarantee — all carried unchanged from Phase 22.1 (RP-10), restated in TI-4/TI-11.

## One-target check, theorem by theorem

| Theorem | Sole target referenced |
|---|---|
| TI-1 definition / TI-2 well-posed | $g^\star_\mu$ |
| TI-5 Bayes optimality | $g^\star_\mu$ (vs regularized risk $R_\mu$) |
| TI-6 continuity | $g^\star_\mu$ |
| TI-7 approximation | $g^\star_\mu$ |
| TI-8 calibration | $g^\star_\mu$ |
| TI-9 consistency (fixed mesh) | $g^\star_{\mu,h}$ |
| TI-10 consistency (refining mesh) | $g^\star_{\mu,0}$, with declared rate $g^\star_{\mu,h_N}\to g^\star_{\mu,0}$ |

Every theorem names exactly one target. TI-9 and TI-10 name different targets *for different declared deployments* (fixed vs refining mesh), each single within its deployment, related by an explicit proved convergence in a named common space — not two targets in one statement, and not a silent identification. The mesh index $h$ joins $\mu$ and $z_H^0$ as part of the deployment declaration (TI-11).

## Residual, named (not a target ambiguity)

Regularized (not unregularized) target under the declared $\mu$; square-root continuity/approximation rate absent (A-GRAD); a design floor $2h$ at fixed mesh, vanishing only along a declared refinement (TI-10); consistency conditional on the declared joint schedule; single-deployment, regression-only scope. None introduces a second target within any deployment.

## Verdict

$$\boxed{\textbf{TARGET\_IDENTITY\_FIXED}}$$

Every theorem in Phase 23.1 refers to one and only one target: $g^\star_\mu$ (or, for the optional mesh-refining deployment, the single continuum target $g^\star_{\mu,0}$ reached at a declared rate). The old barycentric/original-coordinate regularized target is retracted, not preserved; Bayes optimality is stated against the matched regularized risk with no "same-risk" transfer; and the previously unstated consistency typings are supplied. An auditor cannot claim the theory holds two targets: within any declared deployment $(z_H^0,B,\mu,h)$ there is exactly one.

**Verdict: `TARGET_IDENTITY_FIXED`.**
