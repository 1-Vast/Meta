# $L_p^\star$ Constant Merge (Item 1)

> **Status:** Phase-28.1 (final probability and constant closure repair), 2026-08-03. Phases 0–24 provenance; the retained theorems (AL-5 calibration, AL-12 consistency) and the positive-ridge target are unchanged in content — only the two audited alignment points are closed. Audit of record: `../28_final_theory_audit/FINAL_AUDIT.md` (`FINAL_THEORY_NOT_READY`; two theorem-alignment issues). No operator redesign, no target change, no scope enlargement. New results carry **CL-** numbers, tagged **[definition] / [proved] / [declared]**.

---

## 1. The single final constant $L_p^\star$

The audit's point: MR-3's boxed $L_p=\mathrm{Lip}_p[L_0(z,B(z)p)]$ is the **base-risk-only** constant, but the approximation contribution to the **regularized** excess risk also carries the ridge difference $\tfrac\mu2(\|F(z)\|^2-\|g^\star_\mu(z)\|^2)$; MR-4 noted this but left the choice between $L_p$ and $L_p'$ unmade. The choice is now made, once, formally.

**Definition CL-1 (the one theorem constant). [definition]**
$$\boxed{\ L_p^\star\ :=\ L_{\mathrm{base}}\ +\ \mu\,\operatorname{diam}(\Delta_m),\qquad L_{\mathrm{base}}:=\sup_{z\in Z}\ \mathrm{Lip}_p\big[\,L_0(z,B(z)p)\,\big]\ \big(\le L_{\mathrm{Lip}}\kappa_B\big).\ }$$
$L_{\mathrm{base}}$ is the former $L_p$ (base conditional risk, finite by MR-4's bound); $\mu\,\operatorname{diam}(\Delta_m)$ is the ridge Lipschitz constant of $p\mapsto\tfrac\mu2\|p\|^2$ on the compact $\Delta_m$. **The symbols $L_p$ and $L_p'$ are retired; $L_p^\star$ is the only coefficient-to-risk constant used henceforth** (approximation theorem, consistency theorem, symbol index).

**Theorem CL-2 ($L_p^\star$ is the regularized-risk-to-coefficient Lipschitz constant). [proved]**
For all measurable $F,G:Z\to\Delta_m$,
$$\big|R_\mu(F)-R_\mu(G)\big|\ \le\ L_p^\star\,\|F-G\|_{L^\infty}\qquad\Big(\text{and pointwise: } |J_\mu(z,F(z))-J_\mu(z,G(z))|\le L_p^\star\|F(z)-G(z)\|\Big),$$
where $\|F-G\|_{L^\infty}=\sup_z\|F(z)-G(z)\|$.
*Proof.* Pointwise, $J_\mu(z,p)=L_0(z,B(z)p)+\tfrac\mu2\|p\|^2$. The first term is $L_{\mathrm{base}}$-Lipschitz in $p$ (CL-1); the second is $\mu\operatorname{diam}(\Delta_m)$-Lipschitz on $\Delta_m$ (since $\nabla(\tfrac\mu2\|p\|^2)=\mu p$ and $\|p\|\le\operatorname{diam}(\Delta_m)$ on the simplex). Sum of Lipschitz constants: $|J_\mu(z,F(z))-J_\mu(z,G(z))|\le L_p^\star\|F(z)-G(z)\|$. Take $\mathbb E_\zeta$ (or sup): $|R_\mu(F)-R_\mu(G)|\le\mathbb E_\zeta|J_\mu(\zeta,F(\zeta))-J_\mu(\zeta,G(\zeta))|\le L_p^\star\|F-G\|_{L^\infty}$. $\square$

## 2. The approximation-to-excess-risk step, now exact

**Corollary CL-3 (the AL-12 approximation term, closed). [proved]**
With the witness $F_{\omega^\star}\in\mathcal H_N$ (AL-7, $\sup_z\|F_{\omega^\star}-g^\star_\mu\|=\varepsilon_{\mathrm{approx}}(N)$) and CL-2:
$$\inf_{\Omega_N}R_\mu-R_\mu(g^\star_\mu)\ \le\ R_\mu(F_{\omega^\star})-R_\mu(g^\star_\mu)\ \le\ L_p^\star\,\varepsilon_{\mathrm{approx}}(N).$$
This now **follows from the formal definition** (CL-1/CL-2), including the ridge term — the audit's gap ("the bound does not follow from the base-only definition") is closed because $L_p^\star$ contains $\mu\operatorname{diam}(\Delta_m)$ by construction. $\square$

**Declaration CL-4 (retained theorems, re-pointed to $L_p^\star$). [declared]**
- **Approximation theorem:** the excess-risk contribution of approximation is $\le L_p^\star\,\varepsilon_{\mathrm{approx}}(N)$ (CL-3).
- **Consistency theorem (AL-12), final form:** $\mathcal E_\mu(F_{\hat\omega_N})\le 2\Gamma_N+\gamma^{\mathrm{opt}}_N+L_p^\star\,\varepsilon_{\mathrm{approx}}(N)$, fed into calibration $\Phi$.
- **Symbol index:** $L_p$ and $L_p'$ entries replaced by the single $L_p^\star=L_{\mathrm{base}}+\mu\operatorname{diam}(\Delta_m)$ with definition location CL-1.
No other symbol changes; the ambiguity between $L_p$ and $L_p'$ is removed by using $L_p^\star$ exclusively.
