# Stopping Criterion (Final Probability and Constant Closure)

> **Status:** Phase-28.1 terminal decision, 2026-08-03. Sources: the two content files (CL-1–CL-8) and the audit `../28_final_theory_audit/FINAL_AUDIT.md` (`FINAL_THEORY_NOT_READY`; two theorem-alignment issues, theory content passed). Phases 0–24 unmodified; positive-ridge target, operator $\mathsf A(F,z)=K(B(z)F(z))$, and all scope unchanged. Only the two mandated closures performed.

---

## The two issues, closed

| Item | Audit issue | Repair | Where |
|---|---|---|---|
| **1. Merge $L_p$ constants** | $L_p$ defined base-only, but the retained theorem needs the ridge term too; choice between $L_p$ and $L_p'$ left unmade | **One** constant $L_p^\star=L_{\mathrm{base}}+\mu\operatorname{diam}(\Delta_m)$ defined (CL-1), proved to satisfy $|R_\mu(F)-R_\mu(G)|\le L_p^\star\|F-G\|$ (CL-2); the approximation-excess-risk step $\le L_p^\star\varepsilon_{\mathrm{approx}}(N)$ now follows formally (CL-3); $L_p$ and $L_p'$ retired, $L_p^\star$ used in approximation theorem, consistency theorem, and symbol index (CL-4) | `Lp_star_merge.md` |
| **2. Probability schedule** | $\log(1/\delta_N)/N\to0$ wrongly claimed to imply $\delta_N\to0$; latter dropped ($\delta_N=1/2$ witness) | False implication retracted (CL-5); **both** $\delta_N\to0$ **and** $\log(1/\delta_N)/N\to0$ retained and stated explicitly as jointly required for high-probability consistency — first for confidence level $\to1$, second for $\Gamma_N\to0$ (CL-6/CL-7); a.s. clause retains $\sum_N\delta_N<\infty$ plus the log-rate (CL-8) | `probability_schedule_repair.md` |

## The final theorems

- **Target (unchanged):** $g^\star_\mu(z)=\arg\min_{p\in\Delta_m}[L_0(z,B(z)p)+\tfrac\mu2\|p\|^2]$, positive ridge, $\ge\mu$-strongly convex; existence, uniqueness, continuity intact. *(The Phase-28 mandate's displayed $J_\mu$ carried a stray minus sign; per its own checklist "positive ridge everywhere" and this repair's Item-1 constraint "positive ridge target", the positive ridge is authoritative — as in Phase-27.1 MR-1/MR-2, which excluded the negative form by proof.)*
- **Calibration (AL-5, unchanged):** $\|d_{\mathbb M}(F,g^\star_\mu)\|_{L^2(\mu_\zeta)}\le\Phi(\mathcal E_\mu(F))+2h$ for $\mathsf A(F,z)=K(B(z)F(z))$.
- **Consistency (final):** under (S-IID) and the schedule of CL-6 (both $\delta_N\to0$ and $\log(1/\delta_N)/N\to0$, plus $D_N\log(\Lambda N)/N\to0$, mesh$\to0$, $\gamma^{\mathrm{opt}}_N\to0$), with probability $\ge1-\delta_N\to1$,
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(2\Gamma_N+\gamma^{\mathrm{opt}}_N+L_p^\star\,\varepsilon_{\mathrm{approx}}(N)\big)+2h,$$
every term inside $\Phi$ vanishing on the schedule; a.s. eventual $\limsup\le2h$ under CL-8.

## Scope (unchanged, as mandated)

Positive-ridge target; operator $\mathsf A(F,z)=K(B(z)F(z))$; one fixed deployment $\mathcal D$; continuous point-valued affinity regression; fixed output mesh $h$; fixed $z_H^0$; no ranking; no continuum refinement; no varying $z_H$.

## Verdict

$$\boxed{\textbf{FINAL\_THEORY\_CLOSURE\_COMPLETE}}$$

The two theorem-alignment issues are closed with no change to operator, target, or scope: a single coefficient-to-regularized-risk Lipschitz constant $L_p^\star=L_{\mathrm{base}}+\mu\operatorname{diam}(\Delta_m)$ replaces the ambiguous $L_p/L_p'$ pair and makes the approximation-to-excess-risk step follow formally; and the high-probability consistency schedule retains **both** $\delta_N\to0$ (confidence level $\to1$) and $\log(1/\delta_N)/N\to0$ ($\Gamma_N\to0$), with the false implication between them retracted and the a.s. clause's summability condition retained. Every retained theorem is now symbol-closed and correctly conditioned.

**Verdict: `FINAL_THEORY_CLOSURE_COMPLETE`.**
