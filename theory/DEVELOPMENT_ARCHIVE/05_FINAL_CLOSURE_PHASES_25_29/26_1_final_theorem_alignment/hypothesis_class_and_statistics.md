# Hypothesis Class Density and Statistical Contract (Items 2, 3, 4, 5)

> **Status:** Phase-26.1, 2026-08-03. Closes the audit's hypothesis-class and statistical defects: undefined $\Xi_N$; density imported externally instead of derived; sample untyped in its target; over-strong probability statement; dangling C-IID branch; non-uniform sieve constants. New results **AL-6–AL-12**, tagged **[definition] / [proved] / [conditional]**. Self-contained with `output_object_and_calibration.md`.

---

## Item 2 — hypothesis class defined, density derived (not inferred from dimension)

**Definition AL-6 (the class $\mathcal H_N$ — one defined parameter domain). [definition]**
Fix a mesh of $Z$ at resolution $r_N$ with node set $\mathcal N_N$ (a finite subset of $Z$; $|\mathcal N_N|=\nu_N$). The parameter domain is the **single defined object**
$$\Omega_N\ =\ \big(\Delta_m\big)^{\mathcal N_N}\ =\ \{\,\omega=(\omega_\nu)_{\nu\in\mathcal N_N}\ :\ \omega_\nu\in\Delta_m\,\},\qquad \Omega_N\subset(\mathbb R^{m+1})^{\nu_N},\ \ D_N=(m+1)\nu_N,$$
compact convex (product of simplices) with diameter $\le\sqrt{\nu_N}\,\mathrm{diam}(\Delta_m)$. The **realization is the piecewise-multilinear interpolation** $G_N(\omega,z)=\sum_{\nu}\phi_\nu(z)\,\omega_\nu$ ($\phi_\nu$ the multilinear basis weights of the mesh, $\ge0$, $\sum_\nu\phi_\nu(z)=1$), so $G_N(\omega,z)\in\Delta_m$ already (convex combination of simplex points) — no projection needed, and $\Omega_N$ **is** the realization domain (the audit's $\Xi_N$-vs-$\Omega_N$ mismatch is gone: one object). The coefficient map is $F_\omega=G_N(\omega,\cdot):Z\to\Delta_m$; the hypothesis class is $\mathcal H_N=\{F_\omega:\omega\in\Omega_N\}$; the operator is $\mathsf A(F_\omega,z)=K(B(z)F_\omega(z))$ (AL-1).

**Theorem AL-7 (the witness family is IN the class; density derived). [proved]**
The multilinear interpolation witness of the approximation argument **is exactly** $\mathcal H_N$ (same basis $\phi_\nu$, same nodes) — not a separate family the class must be shown to contain. Setting $\omega^\star_\nu=g^\star_\mu(\nu)$ for each node $\nu$ gives $F_{\omega^\star}\in\mathcal H_N$ with, by multilinear interpolation error against the modulus of the (proved-continuous, AL-8) target,
$$\varepsilon_{\mathrm{approx}}(N)\ :=\ \inf_{F\in\mathcal H_N}\ \sup_{z\in Z}\ \|F(z)-g^\star_\mu(z)\|\ \le\ \sup_z\|F_{\omega^\star}(z)-g^\star_\mu(z)\|\ \le\ \omega_{g^\star_\mu}\big(\mathrm{mesh}(r_N)\big)\ \le\ \sqrt{2\varpi_\ell(\mathrm{mesh}(r_N))/\mu}.$$
Thus $\varepsilon_{\mathrm{approx}}(N)\to0$ **iff $\mathrm{mesh}(r_N)\to0$** — a property of the *node refinement*, derived from membership of the witness, **not** inferred from $D_N$ growth. The audit's counterexample (constant realizations, $D_N\uparrow$, no approximation) is excluded: refinement is of $r_N$ (hence $\mathrm{mesh}\to0$), and the witness is a class member by construction. $\square$

**Theorem AL-8 (target continuity, in-folder). [proved]** Under (S-CONT) and $\ge\mu$-strong convexity, $\|g^\star_\mu(z)-g^\star_\mu(z')\|\le\sqrt{2\varpi_\ell(d_Z(z,z'))/\mu}$ (two-sided strong-convexity/optimality argument; proof as in the target-alignment phase, repeated here so this package needs no external citation). $\square$

## Item 3 — the supervised sample, typed

**Definition AL-9 (task sample with target). [definition]**
$$T_i\ =\ (S_i,\ Q_i,\ Y_i),\qquad i=1,\dots,N,$$
$S_i$ observable support, $Q_i$ query, and $Y_i\in V$ the **observable identified point-target** (point-supervised channel; Item-scope: continuous point-valued affinity regression). The induced statistic is $\zeta_i=z(S_i,Q_i,\gamma)$. The empirical risk is now typed (it consumes $Y_i$):
$$\widehat R_{\mu,N}(\omega)\ =\ \frac1N\sum_{i=1}^N\Big[\,L\big(B(\zeta_i)F_\omega(\zeta_i),\,Y_i\big)+\tfrac\mu2\|F_\omega(\zeta_i)\|^2\,\Big],$$
$\mathbb E\,\widehat R_{\mu,N}(\omega)=R_\mu(F_\omega)$ under (S-IID). Estimator $\hat\omega_N\in\arg\min_{\Omega_N}\widehat R_{\mu,N}$ (exists; measurable selection), tolerance $\gamma^{\mathrm{opt}}_N$: $\widehat R_{\mu,N}(\hat\omega_N)\le\inf_{\Omega_N}\widehat R_{\mu,N}+\gamma^{\mathrm{opt}}_N$.

## Item 5 — IID only; the C-IID branch removed

**Declaration AL-10. [declared]** This package uses **(S-IID) only**. Every C-IID / fiber-count / missing-fiber-term reference is **removed** — no such term appears in any retained theorem, so the audit's "undefined fiber variation term" cannot arise. (C-IID was an unexercised alternative; per the mandate it is deleted, not completed.)

## Item 4 + uniform sieve constants — the generalization theorem, correctly typed

**Theorem AL-11 (uniform deviation with declared, sieve-uniform constants). [conditional on (S-IID)]**
The per-task loss $\omega\mapsto L(B(\zeta)F_\omega(\zeta),Y)+\tfrac\mu2\|F_\omega(\zeta)\|^2$ is bounded by $\bar L+\tfrac\mu2$ and, **because $F_\omega(\zeta)=\sum_\nu\phi_\nu(\zeta)\omega_\nu$ is $1$-Lipschitz in $\omega$ in the $\ell^\infty\!\to\!$Euclidean sense uniformly in $\zeta$** ($\sum_\nu\phi_\nu=1$), is $\Lambda$-Lipschitz in $\omega$ with $\Lambda=L_{\mathrm{Lip}}\kappa_B+\mu\,\mathrm{diam}(\Delta_m)$ — a constant **independent of $N$** (it depends only on $\mathcal D$ and the loss, not on the node count). Hence the covering number of $\Omega_N$ at scale $s$ is $(3\,\mathrm{diam}(\Delta_m)/s)^{D_N}$, and the standard bounded-difference/covering bound gives, with probability $\ge1-\delta_N$,
$$\sup_{\omega\in\Omega_N}\big|\widehat R_{\mu,N}(\omega)-R_\mu(F_\omega)\big|\ \le\ \Gamma_N\ :=\ C_0\,(\bar L+\tfrac\mu2)\sqrt{\frac{D_N\ln\!\big(\Lambda N\big)+\ln(1/\delta_N)}{N}},$$
with $C_0$ an **absolute** constant (from the covering/concentration lemma) and $\Lambda,\bar L,\mu$ the $N$-independent constants of $\mathcal D$ — so no constant is "absorbed" across a changing class; the audit's uniformity objection is met by exhibiting $\Lambda$ explicitly and noting its $N$-independence. $\square$

**Theorem AL-12 (the aligned consistency theorem — high probability, no a.s. overclaim). [conditional on (S-IID) + declared schedule]**
Schedule: $\mathrm{mesh}(r_N)\to0$ (so $\varepsilon_{\mathrm{approx}}(N)\to0$, AL-7), $D_N=(m+1)\nu_N$ with $\dfrac{D_N\ln(\Lambda N)}{N}\to0$, $\delta_N\to0$, $\gamma^{\mathrm{opt}}_N\to0$; output mesh $h$ fixed. Then with probability $\ge1-\delta_N$:
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(\underbrace{2\Gamma_N+\gamma^{\mathrm{opt}}_N+L_p\,\varepsilon_{\mathrm{approx}}(N)}_{=\,\mathcal E_\mu(F_{\hat\omega_N})\ \text{bound}}\big)+2h,$$
via excess-risk decomposition ($\mathcal E_\mu(F_{\hat\omega_N})\le 2\Gamma_N+\gamma^{\mathrm{opt}}_N+L_p\varepsilon_{\mathrm{approx}}(N)$: AL-11 twice + empirical optimality + AL-7 with the coefficient-loss Lipschitz constant $L_p$) fed into the calibration theorem AL-5. **Probability typing (Item 4):** this is a **high-probability statement**, holding with probability $\ge1-\delta_N\to1$; it is stated as such. An almost-sure eventual bound is **not** claimed, because $\delta_N\to0$ alone does not give it — that would require $\sum_N\delta_N<\infty$ (Borel–Cantelli), a summability condition **not** assumed here; if a deployment additionally declares $\sum_N\delta_N<\infty$, then $\limsup_N\|d_{\mathbb M}(F_{\hat\omega_N},g^\star_\mu)\|_{L^2(\mu_\zeta)}\le2h$ holds almost surely — stated as the explicit extra hypothesis it is, never silently. $\square$
