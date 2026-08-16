# Calibration, Fully Defined (Item 1)

> **Status:** Phase-25.1, 2026-08-03. Every symbol in the calibration inequality is defined **here or in `preliminaries.md`** — the audit's finding 5 (undefined $\mathcal E_\mu$, undefined empirical risk, undefined $d_{\mathbb M}$/$D_V$, and the open question whether $\|B(z)\|$ is inside $D_V$) is closed. Definitions **FC-7–FC-10**, theorem restatement **FC-11**. Tagged **[definition] / [proved]**.

---

## 1. The population and excess risks

**FC-7 [definition] Population regularized risk.** For a measurable coefficient map $F:Z\to\Delta_m$,
$$R_\mu(F)\ =\ \mathbb E_{\zeta}\big[\,J_\mu(\zeta,F(\zeta))\,\big]\ =\ \mathbb E_{\zeta}\big[\,L_0(\zeta,B(\zeta)F(\zeta))+\tfrac\mu2\|F(\zeta)\|^2\,\big],$$
$\zeta=z(S_T,Q_T,\gamma)$ the random statistic induced by $T\sim P_T$ (FC-12); $\mu_\zeta$ its law on $Z$. Well-defined and finite (bounded loss + bounded ridge on compact $\Delta_m$).

**FC-8 [definition] Regularized excess risk.**
$$\mathcal E_\mu(F)\ =\ R_\mu(F)-R_\mu(g^\star_\mu)\ \ge0,$$
with $R_\mu(g^\star_\mu)=\mathbb E_\zeta[\min_{p\in\Delta_m}J_\mu(\zeta,p)]$ the regularized Bayes risk (attained by $g^\star_\mu$, FC-6). This is the explicit excess-risk object the audit found missing.

## 2. The operator metric and the transfer constant (with $\|B(z)\|$ resolved)

**FC-9 [definition] Operator metric and transfer constant.**
- **Operator metric $d_{\mathbb M}$** on operator values: since the compared maps differ only in the probability object (FC-5), $d_{\mathbb M}(F,g^\star_\mu)$ at a statistic $z$ is the **Hausdorff-$W_1$ distance** between the class objects, $d_H^{W_1}\big(K(B(z)F(z)),\,K(B(z)g^\star_\mu(z))\big)$; the $L^2(\mu_\zeta)$ norm is $\|d_{\mathbb M}(F,g^\star_\mu)\|_{L^2(\mu_\zeta)}=\big(\mathbb E_\zeta\,d_{\mathbb M}(\cdot)^2\big)^{1/2}$.
- **Two stability constants, composed (this is where $\|B(z)\|$ lives — the audit's precise question):**
 1. *band $\to$ class* (CDF-band stability, closed convention): $d_H^{W_1}(K(\beta),K(\beta'))\le C_{\mathrm{cls}}\,\|\beta-\beta'\|_{\mathbb B}+2h$, with $C_{\mathrm{cls}}=D_V^{\mathrm{val}}$ and the additive mesh floor $2h$;
 2. *coefficient $\to$ band* (linear assembly): $\|B(z)p-B(z)p'\|_{\mathbb B}\le\kappa_B\,\|p-p'\|$, $\kappa_B=\sup_z\|B(z)\|_{\mathrm{op}}$ (FC-4).
 **The composed transfer constant is $\displaystyle D_V\ :=\ C_{\mathrm{cls}}\cdot\kappa_B\ =\ D_V^{\mathrm{val}}\,\kappa_B$** — so $\|B(z)\|$ **is** inside $D_V$, explicitly, via the factor $\kappa_B$. This resolves the audit's open point: a coefficient-to-band bound needs the assembly-norm factor, and $D_V$ carries it.

**FC-10 [definition] Calibration function.**
$$\Phi(t)\ =\ D_V\,\sqrt{\tfrac{2t}{\mu}}\ =\ D_V^{\mathrm{val}}\,\kappa_B\,\sqrt{\tfrac{2t}{\mu}},\qquad \Phi:[0,\infty)\to[0,\infty),\ \ \Phi(0)=0,\ \Phi\text{ continuous, increasing},\ \Phi(t)\to0\ \text{as}\ t\to0.$$

## 3. The calibration theorem, every symbol defined

**Theorem FC-11 (calibration, self-contained). [proved]**
For any measurable $F:Z\to\Delta_m$:
$$\boxed{\ \big\|d_{\mathbb M}(F,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(\mathcal E_\mu(F)\big)\ +\ 2h,\qquad \Phi(t)=D_V\sqrt{2t/\mu},\ \ D_V=D_V^{\mathrm{val}}\kappa_B.\ }$$
*Proof (all steps from FC-definitions).* (i) *Coefficient control:* $J_\mu(z,\cdot)$ is $\ge\mu$-strongly convex (FC-6) and $g^\star_\mu(z)$ its minimizer, so $J_\mu(z,F(z))-J_\mu(z,g^\star_\mu(z))\ge\tfrac\mu2\|F(z)-g^\star_\mu(z)\|^2$; take $\mathbb E_\zeta$ and use FC-7/FC-8: $\tfrac\mu2\|F-g^\star_\mu\|_{L^2(\mu_\zeta)}^2\le\mathcal E_\mu(F)$, i.e. $\|F-g^\star_\mu\|_{L^2(\mu_\zeta)}\le\sqrt{2\mathcal E_\mu(F)/\mu}$. (ii) *Transfer, pointwise:* by FC-9, $d_{\mathbb M}(F,g^\star_\mu)(z)\le D_V^{\mathrm{val}}\kappa_B\|F(z)-g^\star_\mu(z)\|+2h=D_V\|F(z)-g^\star_\mu(z)\|+2h$. (iii) *$L^2$ norm:* by Minkowski, $\|d_{\mathbb M}(F,g^\star_\mu)\|_{L^2(\mu_\zeta)}\le D_V\|F-g^\star_\mu\|_{L^2(\mu_\zeta)}+2h\le D_V\sqrt{2\mathcal E_\mu(F)/\mu}+2h=\Phi(\mathcal E_\mu(F))+2h$. $\square$
Every symbol — $R_\mu$ (FC-7), $\mathcal E_\mu$ (FC-8), $d_{\mathbb M}$ and $D_V$ (FC-9, with $\|B(z)\|$ inside via $\kappa_B$), $\Phi$ (FC-10), $h$ (FC-1), $\mu$ (FC-1), $g^\star_\mu$ (FC-6) — is defined in this package. The design floor $2h$ is the fixed positive constant of $\mathcal D$; it is not driven to zero.
