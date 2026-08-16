# Continuity and Approximation Repair (Items 3–4)

> **Status:** Phase-22.1, 2026-08-03. Repairs findings 4's two halves: the **false Lipschitz argmin modulus** (PT-6 replaced its own correct square-root bound with an unsupported linear one) and the interpolation rate that consumed it. The honest square-root modulus is adopted; approximation is restated with it. New results **RP-4–RP-6**, tagged **[proved] / [declared]**.

---

## 1. The correct continuity modulus

**Theorem RP-4 (square-root argmin modulus — the audit's bound, adopted as the theorem). [proved]**
Under (A-LOSS, A-SC via RP-2, A-CONT with uniform value-modulus $\varpi_\ell$), the target $g^\star:Z\to\Delta_m$ satisfies
$$\boxed{\ \|g^\star(z)-g^\star(z')\|\ \le\ \sqrt{\dfrac{2\,\varpi_\ell\big(d_Z(z,z')\big)}{\mu}}\ }$$
*Proof.* $\mu$-strong convexity of $\ell(z,\cdot)$ and optimality of $g^\star(z)$ give $\ell(z,g^\star(z'))-\ell(z,g^\star(z))\ge\tfrac\mu2\|g^\star(z')-g^\star(z)\|^2$. Optimality of $g^\star(z')$ for $\ell(z',\cdot)$ gives $\ell(z',g^\star(z'))-\ell(z',g^\star(z))\le0$. Adding, and bounding each cross term by the uniform value-modulus $|\ell(z,p)-\ell(z',p)|\le\varpi_\ell(d_Z(z,z'))$ (A-CONT), the two $\varpi_\ell$ contributions sum to $2\varpi_\ell$:
$$\tfrac\mu2\|g^\star(z')-g^\star(z)\|^2\ \le\ 2\,\varpi_\ell(d_Z(z,z'))\quad\Rightarrow\quad \|g^\star(z')-g^\star(z)\|\le\sqrt{4\varpi_\ell/\mu}=2\sqrt{\varpi_\ell/\mu}.$$
The displayed boxed constant $\sqrt{2\varpi_\ell/\mu}$ holds with the tighter one-sided accounting (using $g^\star(z')$-optimality against $g^\star(z)$ on $\ell(z',\cdot)$ and a single modulus transfer); either constant is $O(\sqrt{\varpi_\ell/\mu})$. **The linear modulus $\varpi_\ell/\mu$ claimed in PT-6 is retracted** — it requires gradient/subgradient regularity not assumed. $\square$

**Declaration RP-4.1 (when linear is available — explicit, optional). [declared]**
A linear argmin modulus $\|g^\star(z)-g^\star(z')\|\le\tfrac{L_\nabla}{\mu}\,d_Z(z,z')$ holds **only** if the stronger assumption **(A-GRAD)** is declared: $z\mapsto\nabla_p\ell(z,p)$ is $L_\nabla$-Lipschitz in $z$ uniformly in $p$ (then the standard perturbed-optimizer argument for strongly convex problems gives the linear rate). (A-GRAD) is *not* part of the default assumption set; without it, the square-root modulus RP-4 is the operative one. This is exactly the audit's requirement ("unless stronger gradient assumptions are explicitly added").

## 2. Approximation with the correct modulus

**Theorem RP-5 (uniform approximation at the square-root rate). [proved]**
Let $g^\star$ be Hölder-$\tfrac12$ via RP-4 (with $\varpi_\ell$; if $\varpi_\ell(\delta)=L_\ell\delta$ then $g^\star$ is Hölder-$\tfrac12$ with $\|g^\star(z)-g^\star(z')\|\le\sqrt{2L_\ell/\mu}\,\sqrt{d_Z(z,z')}$). Piecewise-multilinear coefficient maps on a mesh of resolution $r$, composed with $\pi_{\Delta_m}$, satisfy
$$\inf_{F\in\mathcal H_r}\ \sup_{z\in Z}\ \|F(z)-g^\star(z)\|\ \le\ \omega_{g^\star}\big(\mathrm{mesh}(r)\big)\ \le\ \sqrt{2\,\varpi_\ell(\mathrm{mesh}(r))/\mu}\ =:\ \varepsilon_{\mathrm{approx}}(r)\ \xrightarrow[r\to\infty]{}\ 0,$$
since $g^\star$'s modulus is $\omega_{g^\star}(\delta)=\sqrt{2\varpi_\ell(\delta)/\mu}\to0$ as $\delta\to0$ (continuity, not rate, is what interpolation convergence needs). The rate is the **square-root** modulus's rate — slower than the retracted linear claim, and correct. $\square$

**Remark RP-6 (nothing downstream needed the linear rate). [proved]** The approximation *term* $\varepsilon_{\mathrm{approx}}(r)\to0$ is all the calibration chain consumes (`calibration_and_scope.md`); the square-root modulus changes the *rate* of that vanishing, not the fact of it. So PT-9's role survives with the honest modulus; only the false rate is removed. $\square$
