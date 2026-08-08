# Calibration and Scope Repair (Items 5–6)

> **Status:** Phase-22.1, 2026-08-03. Repairs findings 4 (calibration: omitted mesh floor; unconstrained sieve growth) and item 6 (scope). New results **RP-7–RP-10**, tagged **[proved] / [declared]**.

---

## 1. Calibration with the design floor included (Item 5)

**Theorem RP-7 (coefficient calibration — unchanged, valid under genuine strong convexity). [proved]**
By RP-2 (genuine $\mu$-strong convexity) and pointwise separability (PT-7), for any measurable $F:Z\to\Delta_m$,
$$\tfrac\mu2\,\|F-g^\star\|_{L^2(\mu_\zeta)}^2\ \le\ R(F)-R(g^\star)=\mathcal E(F)\quad\Rightarrow\quad \|F-g^\star\|_{L^2(\mu_\zeta)}\le\sqrt{2\mathcal E(F)/\mu}.$$
This step the audit accepted as valid; it is retained verbatim (now genuinely justified because strong convexity is genuine).

**Theorem RP-8 (operator calibration with the explicit design floor). [proved]**
The transfer to the operator metric carries the **additive Route-B mesh floor** $\varepsilon_{\mathrm{design}}=2h$ (mesh $h$ in value units, from the CDF-band stability $d_H^{W_1}\le\varepsilon D_V+2h$ — Phase-15 DM-11, Phase-19.5 DT-8):
$$\boxed{\ \big\|d_{\mathbb M}(F,g^\star)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(\mathcal E(F)\big)\ +\ \varepsilon_{\mathrm{design}},\qquad \Phi(t)=D_V\sqrt{2t/\mu},\ \ \varepsilon_{\mathrm{design}}=2h.\ }$$
**$\Phi(t)\to0$ as $t\to0$, but the total operator error does not vanish unless $\varepsilon_{\mathrm{design}}=2h\to0$ as well** — i.e. only along a mesh-refinement schedule $h\to0$. The audit's objection ("with a fixed positive additive floor, $\Phi(t)$ does not tend to zero") is met by writing the floor into the bound and declaring its vanishing as a *separate* requirement, never folded into $\Phi$. On Route A (finite outcomes) $\varepsilon_{\mathrm{design}}=0$ and the floor is absent; the floor is the honest price of the continuous scalar output. $\square$

**Theorem RP-9 (end-to-end chain with an explicit sieve schedule). [conditional on the declared schedule]**
The audit noted PT-11.1 let $\varepsilon_{\mathrm{approx}},\Gamma_N,h$ vary freely without a joint schedule. Impose the **declared sieve schedule** linking mesh $r_N$, dimension $\dim\Omega_N$, and $N$:
$$\dim\Omega_N\ \to\infty,\qquad \frac{\dim\Omega_N\,\ln N}{N}\ \to\ 0,\qquad \mathrm{mesh}(r_N),\,h_N\ \to\ 0,$$
(e.g. $\dim\Omega_N=\lfloor N^{1/2}\rfloor$, $h_N=N^{-1/4}$ — any schedule meeting the display). Then, with probability $\to1$,
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},g^\star)\big\|_{L^2(\mu_\zeta)}\ \le\ \underbrace{\Phi\big(L_p\,\varepsilon_{\mathrm{approx}}(r_N)+2\Gamma_N(\dim\Omega_N)+\gamma^{\mathrm{opt}}_N\big)}_{\to0}\ +\ \underbrace{2h_N}_{\to0}\ \longrightarrow\ 0,$$
because each interior term vanishes on the schedule ($\varepsilon_{\mathrm{approx}}(r_N)\to0$ by RP-5; $\Gamma_N=C\bar L\sqrt{(\dim\Omega_N\ln N+\ln1/\delta)/N}\to0$ by the schedule; $\gamma^{\mathrm{opt}}_N\to0$ declared) and $\Phi$ is continuous at $0$, and the floor $2h_N\to0$ separately. **Consistency now holds on an explicit joint schedule, not by simultaneous free limits** — the audit's sieve-growth requirement, supplied. $\square$

## 2. Scope restriction (Item 6)

**Declaration RP-10 (exact scope). [declared]**
$$\boxed{\ \text{This theory covers continuous point-valued affinity regression only. It provides no ranking guarantees.}\ }$$
Precisely:
- **In scope:** a single fixed deployment ($z_H$ constant, Item 1); observable point-valued identified affinity targets in a declared bounded interval $V$; support/query-conditioned statistic $z(S,Q,\gamma)$; Route-B $W_1$-closed CDF-band outputs at declared mesh; genuine $\mu$-strong convexity (RP-2); risk-field continuity (A-CONT); the calibration + scheduled-sieve chain (RP-8/RP-9). When a ridge supplies $\mu$, the target is the **regularized** affinity regressor $g^\star_\mu$ (declared bias, Phase-21 PT-12), *not* the unregularized regression target.
- **Out of scope (no theorem claimed):** any ranking target, joint ordering object, ranking loss, or ranking calibration. Pairwise/listwise/metric-ranking objectives generally violate the uniqueness and strong-convexity conditions (RP-2 does not apply to them), so **none of RP-1–RP-9 is asserted for ranking**. Ranking remains, as in Phase 19.5, a separately supervised Route-A objective outside this repair — not derived from affinity regression (the frozen impossibility DT-7 stands). Varying-$z_H$ operators are likewise out of scope (Item 1). $\square$
