# Stopping Criterion (Minimal Final Theory Repair)

> **Status:** Phase-27.1 terminal decision, 2026-08-03. Sources: the two content files (MR-1–MR-6) and the audit `../27_final_theory_audit/FINAL_AUDIT.md` (`FINAL_THEORY_NOT_READY`; three blockers). Phases 0–24 unmodified; Phase-26.1 aligned theorems retained; operator not redesigned; scope not enlarged. Only the three mandated minimal repairs performed.

---

## The three blockers, resolved

| Item | Blocker | Repair | Where |
|---|---|---|---|
| **1. Ridge sign** | mandate wording had $-\tfrac\mu2\|p\|^2$; package used $+\tfrac\mu2\|p\|^2$ | Positive ridge affirmed as the intended and retained target (MR-1); the negative-ridge form **explicitly rejected** with proof that it destroys strong convexity and hence existence/uniqueness/continuity/calibration (MR-2). One unambiguous target. | `ridge_sign_alignment.md` |
| **2. Undefined $L_p$** | $L_p$ used in AL-12, never defined | $L_p=\sup_{z,p_1\ne p_2}\dfrac{|L_0(z,B(z)p_1)-L_0(z,B(z)p_2)|}{\|p_1-p_2\|}$ defined, proved finite ($\le L_{\mathrm{Lip}}\kappa_B$), and its role in the approximation-to-excess-risk conversion made explicit (MR-3/MR-4). | `Lp_definition_and_probability_schedule.md` |
| **3. Probability schedule** | $\delta_N\to0$ does not give $\Gamma_N\to0$ ($e^{-N}$ witness) | Schedule condition corrected to $\log(1/\delta_N)/N\to0$ (implies $\delta_N\to0$, excludes the witness); almost-sure clause requires **both** $\log(1/\delta_N)/N\to0$ **and** $\sum_N\delta_N<\infty$ (MR-5/MR-6). | same |

## The retained theorems, now unblocked

- **Target (MR-1/MR-2):** one positive-ridge $g^\star_\mu(z)=\arg\min_p[L_0(z,B(z)p)+\tfrac\mu2\|p\|^2]$; negative-ridge form excluded by proof. Existence, uniqueness, continuity (AL-8), all resting on the positive-$\mu$ strong convexity, stand.
- **Calibration (AL-5, unchanged):** $\|d_{\mathbb M}(F,g^\star_\mu)\|_{L^2(\mu_\zeta)}\le\Phi(\mathcal E_\mu(F))+2h$ for the declared operator $\mathsf A(F,z)=K(B(z)F(z))$.
- **Consistency (AL-12, now grounded):** with $L_p$ defined (MR-3) and the schedule $\log(1/\delta_N)/N\to0$ (MR-5), with probability $\ge1-\delta_N\to1$,
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(2\Gamma_N+\gamma^{\mathrm{opt}}_N+L_p\,\varepsilon_{\mathrm{approx}}(N)\big)+2h,$$
with every term inside $\Phi$ vanishing on the schedule and $\Gamma_N\to0$ now genuinely implied; the a.s. eventual $\limsup\le2h$ available under the two declared conditions (MR-6).

## Nonblocking audit remarks, absorbed

The audit's forward-only correction to AL-7 is adopted: only "$\mathrm{mesh}(r_N)\to0\Rightarrow\varepsilon_{\mathrm{approx}}(N)\to0$" is claimed and used (a constant target can be represented without refinement — the reverse implication is neither claimed nor needed). Operator alignment, hypothesis-class definition, meta-learning typing, and scope were passed and are unchanged.

## Scope (unchanged, as mandated)

One fixed deployment $\mathcal D$; continuous point-valued affinity regression; fixed output mesh $h$; fixed $z_H^0$; no ranking; no continuum refinement; no varying $z_H$.

## Verdict

$$\boxed{\textbf{FINAL\_THEORY\_REPAIR\_COMPLETE}}$$

The three blockers are repaired minimally: the retained target is the positive-ridge $g^\star_\mu$ with the negative-ridge form explicitly excluded by a strong-convexity argument; $L_p$ is a defined, finite coefficient-loss Lipschitz constant used in the excess-risk conversion; and the probability schedule is corrected to $\log(1/\delta_N)/N\to0$, with the almost-sure strengthening gated on both that condition and $\sum_N\delta_N<\infty$. No operator was redesigned, no scope enlarged, and every other passed result is retained unchanged.

**Verdict: `FINAL_THEORY_REPAIR_COMPLETE`.**
