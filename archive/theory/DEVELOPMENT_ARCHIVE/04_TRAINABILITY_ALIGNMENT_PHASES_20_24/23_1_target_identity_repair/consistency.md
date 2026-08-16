# Consistency (Item 3, rewritten to $g^\star_\mu$; typings completed)

> **Status:** Phase-23.1, 2026-08-03. The end-to-end consistency theorem, rewritten to the single target $g^\star_\mu$ and with the two typings the audit flagged as unstated (§5) now supplied: the changing-mesh limiting space, and the tie of interpolation resolution and confidence schedule to $\dim\Omega_N$. New results **TI-9–TI-11**, tagged **[proved] / [declared]**.

---

## 1. The fixed-target, fixed-mesh consistency (primary claim)

To keep **one** target across the sequence, the primary theorem holds the deployment — including the output mesh $h$ — **fixed**, matching Item 5's fixed-deployment scope. Then the target $g^\star_\mu$ is one function throughout, and only the estimator varies.

**Theorem TI-9 (consistency to $g^\star_\mu$ at fixed deployment). [conditional on the declared schedule]**
Fix $(z_H^0,B,\mu,h)$. Impose the declared sieve schedule tying resolution to dimension to sample size to confidence:
$$\dim\Omega_N\to\infty,\quad r_N\ \text{with node-count}(r_N)=\dim\Omega_N\ \text{(interpolation resolution tied to the actual parameter dimension)},\quad \frac{\dim\Omega_N\ln N}{N}\to0,\quad \delta_N\to0\ \text{with}\ \tfrac{\ln(1/\delta_N)}{N}\to0.$$
Then with probability $\ge1-\delta_N\to1$,
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(L_p\,\varepsilon_{\mathrm{approx}}(r_N)+2\Gamma_N(\dim\Omega_N,\delta_N)+\gamma^{\mathrm{opt}}_N\big)\ +\ 2h,$$
with $\varepsilon_{\mathrm{approx}}(r_N)\to0$ (TI-7, tied to $\dim\Omega_N$ via node-count), $\Gamma_N=C\bar L\sqrt{(\dim\Omega_N\ln N+\ln(1/\delta_N))/N}\to0$ (schedule), $\gamma^{\mathrm{opt}}_N\to0$ (declared). Hence the interior of $\Phi$ $\to0$ and, $\Phi$ being continuous at $0$,
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \longrightarrow\ 2h\qquad(\text{the fixed, declared design floor}).$$
At fixed mesh the estimator converges to $g^\star_\mu$ **up to the declared floor $2h$** — one target, no changing space, every symbol $g^\star_\mu$. The audit's two unstated typings are closed: $r_N$ is tied to $\dim\Omega_N$ by node-count, and $\delta_N\to0$ is declared. $\square$

## 2. Driving the floor to zero (optional, with the limiting space typed)

**Theorem TI-10 (mesh refinement with a defined limiting target and space). [declared / proved]**
To send $2h\to0$ one refines the output mesh $h_N\to0$. This changes the value-space grid, hence the Route-B law space, so a **common limiting operator space and a coherent target family must be typed** (the audit's first missing typing). Declare it as follows:
- **Nested meshes:** $h_{N+1}$ refines $h_N$ (dyadic), so the CDF-band value grids are nested and the finite-grid law classes embed isometrically into $(\Delta(V),W_1)$ — the **common limiting space is $(\Delta(V),W_1)$ itself** (compact), into which every mesh-$h_N$ operator value maps by the canonical inclusion of its band class.
- **Coherent target family:** $g^\star_{\mu,h_N}$ is the target at mesh $h_N$; the **single limiting target** is $g^\star_{\mu,0}=\arg\min_p[L_0(z,Bp)+\tfrac\mu2\|p\|^2]$ evaluated with the *continuum* value functional, and $d_{\mathbb M}\big(K_{h_N}(\beta),K_0(\beta)\big)\le 2h_N\to0$ for every band (the mesh-floor bound is exactly the distance between the grid-quantized and continuum band classes). So $g^\star_{\mu,h_N}\to g^\star_{\mu,0}$ in the common space at rate $2h_N$.
Then, along the joint schedule of TI-9 with $h_N\to0$,
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_{\mu,0})\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi(\cdots)+2h_N\ \longrightarrow\ 0.$$
**Still one target per statement:** TI-9's target is $g^\star_{\mu,h}$ (fixed mesh); TI-10's target is $g^\star_{\mu,0}$ (continuum), and the two are connected by an explicit, declared convergence $g^\star_{\mu,h_N}\to g^\star_{\mu,0}$ — not conflated, but related by a proved rate in a named common space. Each theorem names its own single target unambiguously. $\square$

**Declaration TI-11 (which target a deployment uses). [declared]**
A fixed-mesh deployment uses $g^\star_{\mu,h}$ (TI-9); a mesh-refining deployment uses $g^\star_{\mu,0}$ (TI-10). Both are regularized ($\mu$-ridge in barycentric coordinates) risk-optimal targets, never the unregularized Bayes target (TI-4). Within any single deployment there is exactly one target. The mesh index $h$ is part of the deployment declaration, exactly like $\mu$ and $z_H^0$.
