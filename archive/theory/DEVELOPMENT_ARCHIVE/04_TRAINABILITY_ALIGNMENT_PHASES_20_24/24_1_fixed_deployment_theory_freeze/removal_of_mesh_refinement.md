# Removal of Mesh-Refinement Claims (Item 1)

> **Status:** Phase-24.1 (fixed deployment theory freeze repair), 2026-08-03. Phases 0–21 unmodified; the Phase-23.1 fixed-mesh results (TI-1–TI-9, minus the TI-9 limsup wording) are retained; the mesh-refinement content is deleted, not fixed. Audit of record: `../24_final_target_identity_audit/FINAL_AUDIT.md` (`FINAL_THEORY_INVALID`; sole obstruction: TI-10's unproved mesh-target convergence). No mesh theory is extended, no argmin-stability is added, no new asymptotic claim is made. New results carry **FD-** numbers, tagged **[retracted] / [declared] / [proved]**.

---

## 1. What is removed (Item 1)

**Retraction FD-1. [retracted]** The following Phase-23.1 statements are withdrawn in full and are not part of the final theory:
- **TI-10** in its entirety — the mesh-indexed target family $g^\star_{\mu,h_N}$, the continuum target $g^\star_{\mu,0}$, the asserted convergence $g^\star_{\mu,h_N}\to g^\star_{\mu,0}$, the common limiting space construction, and the zero-error continuum consistency conclusion $\big\|d_{\mathbb M}(F_{\hat\omega_N},g^\star_{\mu,0})\big\|_{L^2(\mu_\zeta)}\to0$.
- **TI-11**'s "mesh-refining deployment" clause and the phrase "$h$ joins $\mu$ and $z_H^0$ as part of the deployment declaration" is **kept** (mesh is declared and fixed), but its use to license a *refining* sequence is withdrawn.
- Any sentence, in any Phase-23.1 file, asserting that the design floor $2h$ can be driven to $0$, or that a continuum/zero-error target is reached.

**Reason (the audit's point, adopted). [proved]** The estimate $d_{\mathbb M}(K_h(\beta),K_0(\beta))\le 2h$ compares the discretized and continuum operators **for the same band $\beta$**; it does **not** bound $d_{\mathbb M}(K_h(\beta^\star_{\mu,h}),K_0(\beta^\star_{\mu,0}))$, because the two minimizers differ. Closing that gap would require a uniformly convergent objective family $J_{\mu,h}\to J_{\mu,0}$ plus an argmin-stability theorem (yielding coefficient displacement $O(\sqrt{h/\mu})$, not $2h$) — precisely the mesh theory / stability argument the mandate forbids adding. The only mathematically honest action, given both the gap and the prohibition, is to **remove the claim**. The final theory makes no continuum statement.

## 2. The TI-9 wording correction (the smaller overstatement)

**Correction FD-2. [proved]** Phase-23.1 TI-9 wrote the fixed-mesh error as "$\to 2h$". The upper bound proves only
$$\limsup_{N\to\infty}\ \big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ 2h,$$
not convergence to exactly $2h$. The final statement (`fixed_deployment_theory.md`, FD-6) uses this $\limsup\le 2h$ form. No other change to TI-9.

## 3. What remains (pointer)

Everything the audit passed is retained and re-collected, unchanged in content, in `fixed_deployment_theory.md`: the single regularized target $g^\star_\mu$ (TI-1/2), regularized Bayes optimality (TI-5), strong convexity (TI-2/RP-2), continuity at the square-root modulus (TI-6), approximation (TI-7), calibration with the fixed design floor (TI-8), fixed-mesh consistency with the corrected $\limsup$ (TI-9 → FD-6), and the regression-only fixed-deployment scope. Nothing in the retained set depends on TI-10/TI-11's removed content — the fixed-mesh chain never used the continuum limit, so removal leaves it intact.
