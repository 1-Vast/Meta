# Calibration Theorem (Route B)

> **Status:** Phase-21, 2026-08-03. The deliverable the audit demanded: **operator error $\le\Phi(\text{excess population risk})$ with $\Phi(t)\to0$**, for the single target $g^\star$, plus the validated empirical→population→operator chain. New results **PT-7, PT-8, PT-10, PT-11**, tagged **[proved] / [conditional]**. Load-bearing assumption: (A-SC) strong convexity — the same $\mu$ that made $g^\star$ well-defined and Lipschitz now makes excess risk control coefficient error.

---

## 1. Pointwise separability of the risk

**Theorem PT-7 (risk factorizes through the local risk). [proved]**
For any measurable coefficient map $F:Z\to C$, under (A-STAT, A-LOSS, A-SC):
$$R(F)\ :=\ \mathbb E_T\big[L(\mathsf{asm}(F(\zeta);\zeta),A_T)+\tfrac\mu2\|F(\zeta)\|^2\big]\ =\ \mathbb E_\zeta\big[\ell(\zeta,F(\zeta))\big].$$
*Proof.* Tower property over $\zeta=z(S_T,Q_T,\gamma)$: condition on $\zeta$, use $\ell_0(z,c)=\mathbb E[L(\mathsf{asm}(c;z),A_T)\mid\zeta=z]$ and the $\zeta$-measurability of $F(\zeta)$. $\square$

## 2. Two-sided local bounds (both from declared regularity, no smoothness assumed)

**Lemma PT-8 (strong-convexity lower bound and Lipschitz upper bound). [proved]**
For every $z\in Z$ and $c\in C$, with $g^\star(z)=\arg\min_c\ell(z,c)$:
$$\text{(lower)}\quad \ell(z,c)-\ell(z,g^\star(z))\ \ge\ \tfrac\mu2\|c-g^\star(z)\|^2\qquad\text{[}\mu\text{-strong convexity + optimality, A-SC]};$$
$$\text{(upper)}\quad \ell(z,c)-\ell(z,g^\star(z))\ \le\ L_c\,\|c-g^\star(z)\|\qquad\text{[}\ell(z,\cdot)\text{ is }L_c\text{-Lipschitz on compact }C\text{: convex \& finite ⇒ Lipschitz, A-LOSS + affine }\mathsf{asm}\text{]}.$$
No $\beta$-smoothness is invoked — the upper bound uses only Lipschitzness, so the piecewise-linear interval score is covered. $\square$

## 3. The calibration inequality

**Theorem PT-10 (operator error $\le\Phi(\text{excess risk})$, $\Phi(t)\to0$). [proved]**
For any measurable $F:Z\to C$, integrating PT-8(lower) against $\mu_\zeta$ and using PT-7:
$$\tfrac\mu2\,\mathbb E_\zeta\|F(\zeta)-g^\star(\zeta)\|^2\ \le\ R(F)-R(g^\star)\ =:\ \mathcal E(F),$$
so the coefficient error obeys $\ \|F-g^\star\|_{L^2(\mu_\zeta)}\le\sqrt{2\mathcal E(F)/\mu}$. Transferring through the fixed affine assembly and the operator stability constant $C_{\mathrm{stab}}$ (Hoffman on Route A; $D_V$, plus the declared mesh floor, on Route B):
$$\boxed{\ \big\|\,d_{\mathbb M}\!\big(F,\ g^\star\big)\,\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(\mathcal E(F)\big),\qquad \Phi(t)\ =\ C_{\mathrm{stab}}\sqrt{\tfrac{2t}{\mu}}\ \xrightarrow[t\to0]{}\ 0.\ }$$
Under (A-DESIGN) with mesh mass floor $q_0$, the $L^2(\mu_\zeta)$ operator error upgrades to a sup bound $\le\Phi(\mathcal E(F))/\sqrt{q_0}+$ mesh-modulus; without it, the calibrated statement is $L^2(\mu_\zeta)$ (average-case), stated as such. **This is the exact inequality the audit required: excess population risk controls operator-metric error through an explicit $\Phi$ vanishing at $0$ — for the risk-optimal target, not a canonical surrogate.** $\square$

## 4. The full valid chain: empirical task loss → population risk → operator error

**Theorem PT-11 (end-to-end, one target). [conditional on the declared stack]**
Let $\hat\omega_N$ minimize the empirical task risk $\widehat R_N$ over the compact finite-dimensional parameter space $\Omega$ (anchors fixed; only the coefficient-map parameters vary), to optimization tolerance $\gamma^{\mathrm{opt}}$. Then, with probability $\ge1-\delta$:
$$\underbrace{\big\|d_{\mathbb M}(F_{\hat\omega_N},g^\star)\big\|_{L^2(\mu_\zeta)}}_{\text{operator error to the single target}}\ \le\ \Phi\Big(\underbrace{L_c\,\varepsilon_{\mathrm{approx}}}_{\text{PT-8 upper, PT-9}}\ +\ \underbrace{2\Gamma_N}_{\text{generalization}}\ +\ \underbrace{\gamma^{\mathrm{opt}}}_{\text{optimization}}\Big),$$
where $\varepsilon_{\mathrm{approx}}=\inf_{F\in\mathcal H}\sup_z\|F-g^\star\|\to0$ by PT-9, and $\Gamma_N=C\bar L\sqrt{(\dim\Omega\,\ln(\mathrm{Lip}\,N)+\ln(1/\delta))/N}$ is the covering-number generalization gap over compact $\Omega$ (task-(IID)/(C-IID-$\kappa$), bounded Lipschitz loss; fiber-relative with the missing-fiber term under C-IID; DE-T3 reversal under undeclared shift, tagged).
*Proof.* Excess risk decomposes: $\mathcal E(F_{\hat\omega_N})=[R(F_{\hat\omega_N})-\inf_{\mathcal H}R]+[\inf_{\mathcal H}R-R(g^\star)]$. The first bracket $\le 2\Gamma_N+\gamma^{\mathrm{opt}}$ (uniform deviation + empirical optimality + tolerance). The second $\le L_c\,\varepsilon_{\mathrm{approx}}$ (PT-8 upper bound applied to the PT-9 witness, integrated). Feed the sum into PT-10's $\Phi$. Every term is a proved quantity or a tagged declaration; the target throughout is the single $g^\star$. $\square$

**Corollary PT-11.1 (consistency).** As $N\to\infty$, $\varepsilon_{\mathrm{approx}}\to0$ (family refined per PT-9), $\Gamma_N\to0$, $\gamma^{\mathrm{opt}}\to0$: the operator error to $g^\star$ tends to $0$ in $L^2(\mu_\zeta)$ (in sup under A-DESIGN). Empirical task learning provably converges to the risk-optimal operator in the operator metric — the claim Phase 20 asserted but could not prove, now proved for one explicitly named, continuity-derived target. $\square$
