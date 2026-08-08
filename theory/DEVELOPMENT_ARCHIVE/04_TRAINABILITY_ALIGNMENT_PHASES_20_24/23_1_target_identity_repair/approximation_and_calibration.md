# Approximation and Calibration (Item 3, rewritten to $g^\star_\mu$)

> **Status:** Phase-23.1, 2026-08-03. The approximation and calibration theorems, rewritten so every symbol refers to the single target $g^\star_\mu$. The mathematics is the Phase-22.1 content (which the audit passed on approximation, and passed on the coefficient calibration step); only the target symbol is unified and the "same risk" language purged. New results **TI-6–TI-8**, tagged **[proved]**.

---

## 1. Continuity of $g^\star_\mu$ (square-root modulus, carried and re-pointed)

**Theorem TI-6 (modulus of the one target). [proved]**
Under (A-LOSS, A-SC, A-CONT with uniform value-modulus $\varpi_\ell$), the single target $g^\star_\mu$ satisfies
$$\|g^\star_\mu(z)-g^\star_\mu(z')\|\ \le\ \sqrt{\tfrac{2\,\varpi_\ell(d_Z(z,z'))}{\mu}}.$$
*Proof.* The two-sided strong-convexity argument (RP-4), applied to $J_\mu(z,\cdot)$ and $J_\mu(z',\cdot)$ — both $\ge\mu$-strongly convex (TI-2) — with the two uniform value-modulus transfers: $\mu\|g^\star_\mu(z)-g^\star_\mu(z')\|^2\le2\varpi_\ell(d_Z(z,z'))$, giving the boxed bound. The audit accepted this boxed conclusion from the stated assumptions (§4). Linear modulus only under the optional declared (A-GRAD); not default. $\square$

## 2. Approximation of $g^\star_\mu$

**Theorem TI-7 (uniform approximability of the one target). [proved]**
Piecewise-multilinear coefficient maps on a mesh of resolution $r$, composed with $\pi_{\Delta_m}$, satisfy
$$\inf_{F\in\mathcal H_r}\ \sup_{z\in Z}\ \|F(z)-g^\star_\mu(z)\|\ \le\ \omega_{g^\star_\mu}(\mathrm{mesh}(r))\ \le\ \sqrt{2\varpi_\ell(\mathrm{mesh}(r))/\mu}\ =:\ \varepsilon_{\mathrm{approx}}(r)\ \xrightarrow[r\to\infty]{}\ 0,$$
by setting node values to $g^\star_\mu$ at the nodes (in $\Delta_m$, so $\pi$ is identity there) and bounding interpolation error by $g^\star_\mu$'s modulus (TI-6). Every occurrence of the target here is $g^\star_\mu$; no canonical operator, no old target, no "same-risk" transfer. $\square$

## 3. Calibration to $g^\star_\mu$, with the design floor explicit

**Theorem TI-8 (calibration inequality for the one target). [proved]**
For any measurable $F:Z\to\Delta_m$, using $\ge\mu$-strong convexity of $J_\mu(z,\cdot)$ and the pointwise separability of $R_\mu$ (TI-5):
$$\tfrac\mu2\|F-g^\star_\mu\|_{L^2(\mu_\zeta)}^2\ \le\ R_\mu(F)-R_\mu(g^\star_\mu)\ =:\ \mathcal E_\mu(F)\quad\Rightarrow\quad \|F-g^\star_\mu\|_{L^2(\mu_\zeta)}\le\sqrt{2\mathcal E_\mu(F)/\mu},$$
and, transferring through the fixed linear assembly $B$ and the Route-B stability constant with the **additive mesh floor** written separately (carried from RP-8, unchanged mathematics):
$$\boxed{\ \big\|d_{\mathbb M}(F,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(\mathcal E_\mu(F)\big)+\varepsilon_{\mathrm{design}},\qquad \Phi(t)=D_V\sqrt{2t/\mu}\xrightarrow[t\to0]{}0,\ \ \varepsilon_{\mathrm{design}}=2h.\ }$$
Excess risk $\mathcal E_\mu$ is measured against $R_\mu(g^\star_\mu)$ — the regularized optimum, matched to the target. Total operator error vanishes only when $\mathcal E_\mu\to0$ **and** $h\to0$ (Item-5 scope: $h$ is the declared output mesh; its vanishing is a separate schedule requirement, `consistency.md`). Every quantity in this inequality refers to the single $g^\star_\mu$; there is no second regularized minimizer anywhere in the statement. $\square$
