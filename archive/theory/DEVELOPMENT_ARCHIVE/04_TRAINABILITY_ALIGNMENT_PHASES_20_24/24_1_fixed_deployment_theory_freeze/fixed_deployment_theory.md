# Fixed Deployment Theory (Items 2–4)

> **Status:** Phase-24.1, 2026-08-03. The complete final theory on **one fixed finite deployment**, re-collecting the Phase-23.1 fixed-mesh results the audit passed, with the continuum content removed (FD-1) and the TI-9 wording corrected (FD-2). New results **FD-3–FD-6**, each a re-statement of a passed Phase-23.1 theorem restricted to the fixed deployment; tagged **[proved] / [declared]**.

---

## Item 2 — the fixed deployment skeleton

**Declaration FD-3 (the deployment). [declared]**
The theory is stated for a single fixed tuple
$$\boxed{\ \mathcal D\ =\ \big(\,z_H^0,\ B,\ \Delta_m,\ \mu,\ h\,\big)\ }$$
— frozen deployment state $z_H^0$; fixed band matrix $B=[\beta_0\cdots\beta_m]$ (population band $\beta_0=b^{\mathrm{pop}}_{\kappa(z)}$, fixed anchors $\beta_1,\dots,\beta_m$); fixed coefficient simplex $\Delta_m$; fixed ridge modulus $\mu>0$; fixed Route-B output value mesh $h>0$. $\mathcal D$ is declared once; **nothing in the theory varies any coordinate of $\mathcal D$.** The statistic is $z=z(S_T,Q_T,\gamma)\in Z$ (support/query conditioned, as before).

## Item 3 — the resolution statement

**Declaration FD-4 (fixed-resolution guarantee). [declared]**
$$\textbf{The theory provides guarantees for a fixed finite deployment resolution }\mathcal D.$$
All errors below are relative to $\mathcal D$; in particular the design floor $2h$ is a fixed positive constant of $\mathcal D$, not a quantity driven to zero. No statement quantifies over meshes, over deployments, or over a continuum limit.

## Item 4 — the retained theory (one target, all guarantees)

**FD-5 (the retained theorem set — content identical to the passed Phase-23.1 results). [proved]**
1. **Unique target.** $g^\star_\mu(z)=\arg\min_{p\in\Delta_m}J_\mu(z,p)$, $J_\mu(z,p)=L_0(z,Bp)+\tfrac\mu2\|p\|^2$ — single-valued, everywhere-defined (TI-1/TI-2).
2. **Strong convexity.** $J_\mu(z,\cdot)$ is $\ge\mu$-strongly convex on $\Delta_m$ (linear assembly $Bp$ + convex base loss + Euclidean ridge; RP-2/TI-2).
3. **Regularized Bayes optimality.** $g^\star_\mu$ uniquely (a.e.) minimizes the regularized risk $R_\mu$ over measurable maps $Z\to\Delta_m$; no claim about the unregularized risk $R_0$ or target $g^\star_0$ (TI-5, TI-4).
4. **Continuity.** $\|g^\star_\mu(z)-g^\star_\mu(z')\|\le\sqrt{2\varpi_\ell(d_Z(z,z'))/\mu}$ (square-root modulus; linear only under the optional declared A-GRAD) (TI-6).
5. **Approximation.** $\inf_{F\in\mathcal H_r}\sup_z\|F(z)-g^\star_\mu(z)\|\le\sqrt{2\varpi_\ell(\mathrm{mesh}(r))/\mu}\to0$ as the coefficient-map resolution $r\to\infty$ (TI-7). *(This is the coefficient interpolation resolution $r$ — internal to the hypothesis class at fixed output mesh $h$; it is not the output mesh and does not refine $h$.)*
6. **Calibration.** $\big\|d_{\mathbb M}(F,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\le\Phi(\mathcal E_\mu(F))+2h$, $\Phi(t)=D_V\sqrt{2t/\mu}\to0$ as $t\to0$, with the design floor $2h$ additive and separate (TI-8).

**Theorem FD-6 (fixed-deployment consistency, corrected wording). [conditional on the declared schedule]**
On the declared sieve schedule ($\dim\Omega_N\to\infty$ with coefficient resolution $r_N$ tied to $\dim\Omega_N$ by node-count; $\dim\Omega_N\ln N/N\to0$; $\delta_N\to0$ with $\ln(1/\delta_N)/N\to0$ — all at the **fixed** output mesh $h$), with probability $\ge1-\delta_N\to1$:
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(L_p\,\varepsilon_{\mathrm{approx}}(r_N)+2\Gamma_N+\gamma^{\mathrm{opt}}_N\big)\ +\ 2h,$$
and therefore
$$\limsup_{N\to\infty}\ \big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ 2h.$$
*(Corrected per FD-2: $\limsup\le 2h$, not "$\to2h$".)* The estimator converges to the single target $g^\star_\mu$ up to the fixed, declared design floor $2h$ of $\mathcal D$. Every symbol is $g^\star_\mu$; there is no mesh-indexed or continuum target anywhere in the statement. $\square$

**Scope retained (Item 4 / carried).** Continuous point-valued affinity regression only, on the fixed deployment $\mathcal D$; ridge target $g^\star_\mu$ declared regularized, not the unregularized Bayes target; no ranking guarantee. Unchanged from Phase 22.1/23.1.
