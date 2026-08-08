# Symbol Index (self-containment audit aid)

> **Status:** Phase-25.1, 2026-08-03. Every symbol appearing in any retained theorem, with its definition location **inside this folder**. An auditor confirms self-containment by checking that no cell reads "external".

| Symbol | Meaning | Defined in |
|---|---|---|
| $\mathcal D=(z_H^0,B(\cdot),\Delta_m,\mu,h)$ | fixed deployment | FC-1 |
| $z_H^0$ | frozen deployment state | FC-1 |
| $Z,\ d_Z$ | statistic domain and its metric | FC-2 |
| $z=z(S,Q,\gamma)$ | support/query-conditioned statistic | FC-2 |
| $\kappa,\ C_\kappa$ | context map, finite context set | FC-2 |
| $V,\ D_V^{\mathrm{val}}$ | value interval, its diameter | FC-2 |
| $(\Delta(V),W_1)$ | law space, Wasserstein-1 | FC-2 |
| $C=\Delta_m,\ \|\cdot\|$ | coefficient simplex, Euclidean norm | FC-3 |
| $B(z),\ \beta_k(z),\ \kappa_B$ | band rule, columns, assembly-norm constant | FC-4 |
| $\mathbb B,\ K(\cdot)$ | band polytope, class map | FC-5 |
| $\mathbb M$ | operator value space | FC-5 |
| $L,\ \bar L,\ L_{\mathrm{Lip}}$ | loss, bound, band-Lipschitz constant | FC-6 |
| $L_0(z,\beta)$ | base conditional risk | FC-6 |
| $J_\mu(z,p),\ \mu$ | operative regularized risk, ridge modulus | FC-6 |
| $g^\star_\mu$ | the single target | FC-6 |
| $L_p$ | coefficient-loss Lipschitz constant | FC-6 / FC-18 |
| $R_\mu(F)$ | population regularized risk | FC-7 |
| $\mu_\zeta$ | law of $\zeta$ on $Z$ | FC-7 |
| $\mathcal E_\mu(F)$ | regularized excess risk | FC-8 |
| $d_{\mathbb M},\ d_H^{W_1}$ | operator metric, Hausdorff-$W_1$ | FC-9 |
| $C_{\mathrm{cls}},\ D_V=D_V^{\mathrm{val}}\kappa_B$ | class-stability const., transfer const. (incl. $\|B(z)\|$) | FC-9 |
| $\Phi(t)=D_V\sqrt{2t/\mu}$ | calibration function | FC-10 |
| $h$ | fixed output mesh; floor $2h$ | FC-1 / FC-11 |
| $P_T,\ T_i=(S_i,Q_i),\ A_i$ | task law, task, identified target | FC-12 |
| $\Omega_N,\ D_N,\ G_N,\ F_\omega,\ \mathcal H_N$ | parameter space, dim, realization, coeff map, class | FC-13 |
| $\widehat R_{\mu,N}$ | empirical regularized risk | FC-14 |
| $\hat\omega_N,\ F_{\hat\omega_N}$ | estimator, learned map | FC-15 |
| $\varepsilon_{\mathrm{approx}}(N)$ | approximation term | FC-16 |
| $\Gamma_N,\ \delta_N$ | generalization term, confidence level | FC-17 |
| $\gamma^{\mathrm{opt}}_N$ | optimization tolerance | FC-18 |

No row reads "external": every symbol in the retained calibration theorem (FC-11) and consistency theorem (FC-19) is defined within `25_1_final_theory_closure/`.
