# Symbol Index (self-containment audit aid)

> **Status:** Phase-26.1, 2026-08-03. Every symbol in the two retained theorems (AL-5 calibration, AL-12 consistency) with an in-folder definition. No entry reads "external". The audit's undefined-symbol list ($\Xi_N$, A-STAT, A-CONT, C-IID, DE-T3, interpolation theorem, $I(S)$, conf, rung, external $\varepsilon_{\mathrm{approx}}$) is closed by removal or in-folder definition, as noted.

| Symbol | Meaning | Defined in | Audit item |
|---|---|---|---|
| $\mathcal D=(z_H^0,B(\cdot),\Delta_m,\mu,h)$ | fixed deployment | AL-2 | — |
| $\mathsf A(F,z)=K(B(z)F(z))$ | **declared operator output (no support restriction)** | AL-1 | 1 (Route A) |
| $I(S)$, conf, rung | — | **removed from output/theorems** (AL-1) | 1 |
| $Z,d_Z,\kappa,C_\kappa$ | statistic domain, metric, context map | AL-2 | closure |
| $V,D_V^{\mathrm{val}},(\Delta(V),W_1)$ | value interval, diameter, law space | AL-2 | — |
| $\Delta_m$ | coefficient simplex | AL-2 | — |
| $B(z),\beta_k(z),\kappa_B$ | band rule, columns, assembly-norm const. | AL-2 | 3 (notation) |
| $\mathbb B,K(\cdot)$ | band polytope, class map | AL-2 | — |
| $L,\bar L,L_{\mathrm{Lip}},L_0$ | loss, bound, band-Lipschitz, base cond. risk | AL-2, AL-3 | 1 (S-CONT) |
| (S-IID),(S-CONT),(S-GRID) | the three in-folder assumptions | AL-3 | replaces A-STAT/A-CONT |
| C-IID, DE-T3 | — | **removed** (AL-10) | 5 |
| $\zeta,\mu_\zeta$ | induced statistic, its law | AL-4 | — |
| $J_\mu,g^\star_\mu$ | operative risk, target | AL-4 | 2 (target) |
| $R_\mu,\mathcal E_\mu$ | population risk, excess risk | AL-4 | 4 |
| $d_{\mathbb M},d_H^{W_1}$ | operator metric (between $\mathsf A$'s), Hausdorff-$W_1$ | AL-4 | 1,4 |
| $D_V=D_V^{\mathrm{val}}\kappa_B,\Phi$ | transfer const. (incl. $\|B(z)\|$), calibration fn | AL-4 | 4 |
| $h$; floor $2h$ | fixed output mesh | AL-2 | — |
| $P_T,T_i=(S_i,Q_i,Y_i)$ | task law, **typed supervised sample** | AL-9 | 3 |
| $\Omega_N=(\Delta_m)^{\mathcal N_N},D_N,\nu_N,r_N$ | **one** parameter domain, dim, node count, resolution | AL-6 | 2 ($\Xi_N$ resolved) |
| $G_N,\phi_\nu,F_\omega,\mathcal H_N$ | realization, basis, coeff map, class | AL-6 | 2 |
| $\varepsilon_{\mathrm{approx}}(N)$ | approximation term, **derived in-folder** (witness $\in\mathcal H_N$) | AL-7 | 2 |
| $\widehat R_{\mu,N},\hat\omega_N,\gamma^{\mathrm{opt}}_N$ | empirical risk, estimator, tolerance | AL-9 | 3 |
| $\Lambda,C_0,\Gamma_N,\delta_N$ | param-Lipschitz const ($N$-indep), abs. const, gen. term, confidence | AL-11 | 5 (uniformity) |
| $L_p$ | coefficient-loss Lipschitz constant | AL-4/AL-12 | — |

No row reads "external"; every symbol resolves within `26_1_final_theorem_alignment/`. Density ($\varepsilon_{\mathrm{approx}}\to0$) is **derived** from the witness being a class member (AL-7), not imported.
