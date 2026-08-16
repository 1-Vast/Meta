# Frozen Theory: Symbol Index

The table lists active symbols only. Rejected symbols occur solely in `08_FAILURE_HISTORY.md`.

| Symbol | Definition | First occurrence | Dependency source |
|---|---|---|---|
| $\mathcal D$ | Fixed tuple $(z_H^0,B(\cdot),\Delta_m,\mu,h)$ | `01_FOUNDATIONS.md` section2 | AL-2, FD-3 |
| $z_H^0$ | Frozen deployment state | `01_FOUNDATIONS.md` section2 | DT-A, AL-2 |
| $Z,d_Z$ | Compact statistic domain and metric | `01_FOUNDATIONS.md` section2 | A-SKEL, AL-2 |
| $S,Q,\gamma$ | Support set, query, specification | `01_FOUNDATIONS.md` section1 | AL-2, AL-9 |
| $z(S,Q,\gamma)$ | Observable support/query statistic | `01_FOUNDATIONS.md` section1 | AL-2 |
| $\kappa,C_\kappa$ | Context map and finite context set | `01_FOUNDATIONS.md` section2 | AL-2 |
| $V$ | Compact affinity interval $[a_{\min},a_{\max}]$ | `01_FOUNDATIONS.md` section2 | AL-2 |
| $D_V^{\rm val}$ | Value-space diameter | `01_FOUNDATIONS.md` section2 | AL-2 |
| $h$ | Fixed output-grid mesh | `01_FOUNDATIONS.md` section2 | S-GRID, AL-2 |
| $\Delta_m$ | Coefficient simplex in $\mathbb R^{m+1}$ | `01_FOUNDATIONS.md` section2 | RP-1, AL-2 |
| $\mu$ | Positive ridge modulus | `01_FOUNDATIONS.md` section2 | MR-1 |
| $\mathbb B$ | Compact convex valid CDF-band polytope | `01_FOUNDATIONS.md` section3 | AL-2 |
| $B(z)$ | Deployment band-assembly matrix rule | `01_FOUNDATIONS.md` section3 | FC-4, AL-2 |
| $\beta_k(z)$ | Columns of $B(z)$; population column and anchors | `01_FOUNDATIONS.md` section3 | FC-4, AL-2 |
| $\kappa_B$ | $\sup_z\|B(z)\|_{\rm op}$ | `01_FOUNDATIONS.md` section3 | FC-4, AL-2 |
| $T=(S,Q,Y)$ | Observable supervised task | `01_FOUNDATIONS.md` section4 | AL-9 |
| $P_T$ | IID observable task law | `01_FOUNDATIONS.md` section5 | S-IID, AL-9 |
| $Y$ | Point-valued affinity target in $V$ | `01_FOUNDATIONS.md` section4 | AL-9 |
| $L$ | Bounded convex band loss | `01_FOUNDATIONS.md` section4 | A-LOSS, AL-2 |
| $\bar L$ | Uniform loss bound | `01_FOUNDATIONS.md` section4 | A-LOSS, AL-11 |
| $L_{\rm Lip}$ | Band-argument Lipschitz constant of $L$ | `01_FOUNDATIONS.md` section4 | A-LOSS, AL-2 |
| $\zeta$ | Random statistic $z(S,Q,\gamma)$ | `01_FOUNDATIONS.md` section4 | AL-4 |
| $\mu_\zeta$ | Law of $\zeta$ on $Z$ | `01_FOUNDATIONS.md` section4 | AL-4 |
| $L_0$ | Conditional base risk $\mathbb E[L(\beta,Y)\mid\zeta=z]$ | `01_FOUNDATIONS.md` section4 | S-CONT, AL-3 |
| $\varpi_\ell$ | Uniform continuity modulus of $L_0$ in $z$ | `01_FOUNDATIONS.md` section5 | S-CONT, AL-3 |
| $\Delta(V)$ | Borel probability laws on $V$ | `02_TARGET_AND_OPERATOR.md` section1 | AL-2 |
| $W_1$ | Wasserstein-1 metric on $\Delta(V)$ | `02_TARGET_AND_OPERATOR.md` section1 | AL-2 |
| $K(\beta)$ | Nonempty compact law class induced by band $\beta$ | `02_TARGET_AND_OPERATOR.md` section1 | AL-1, AL-2 |
| $\mathcal K(\Delta(V))$ | Nonempty compact subsets of the law space | `02_TARGET_AND_OPERATOR.md` section1 | Operator codomain |
| $\mathsf A(F,z)$ | Sole operator $K(B(z)F(z))$ | `02_TARGET_AND_OPERATOR.md` section2 | AL-1 |
| $d_H^{W_1}$ | Hausdorff distance induced by $W_1$ | `02_TARGET_AND_OPERATOR.md` section2 | AL-4 |
| $d_{\mathbb M}$ | Pointwise metric between retained operator outputs | `02_TARGET_AND_OPERATOR.md` section2 | AL-4 |
| $D_V$ | Transfer constant $D_V^{\rm val}\kappa_B$ | `02_TARGET_AND_OPERATOR.md` section2 | AL-4, FC-9 |
| $J_\mu$ | Local regularized risk $L_0(z,B(z)p)+\mu\|p\|^2/2$ | `02_TARGET_AND_OPERATOR.md` section3 | AL-4, MR-1 |
| $g_\mu^\star$ | Unique minimizer of $J_\mu(z,\cdot)$ | `02_TARGET_AND_OPERATOR.md` section3 | MR-1 |
| $R_\mu$ | Population regularized risk | `03_STRONG_CONVEXITY_AND_REGULARITY.md` section3 | AL-4 |
| $r_N$ | Coefficient interpolation resolution | `04_META_LEARNING_FORMULATION.md` section2 | AL-6 |
| $\mathcal N_N,\nu_N$ | Mesh node set and node count | `04_META_LEARNING_FORMULATION.md` section2 | AL-6 |
| $\Omega_N$ | Product-simplex parameter space | `04_META_LEARNING_FORMULATION.md` section2 | AL-6 |
| $D_N$ | Ambient parameter dimension $(m+1)\nu_N$ | `04_META_LEARNING_FORMULATION.md` section2 | AL-6 |
| $\phi_\nu$ | Nonnegative multilinear basis weight | `04_META_LEARNING_FORMULATION.md` section2 | AL-6 |
| $G_N$ | Multilinear realization map | `04_META_LEARNING_FORMULATION.md` section2 | AL-6 |
| $F_\omega$ | Coefficient hypothesis $G_N(\omega,\cdot)$ | `04_META_LEARNING_FORMULATION.md` section2 | AL-6 |
| $\mathcal H_N$ | Hypothesis class $\{F_\omega:\omega\in\Omega_N\}$ | `04_META_LEARNING_FORMULATION.md` section2 | AL-6 |
| $\widehat R_{\mu,N}$ | Empirical regularized risk | `04_META_LEARNING_FORMULATION.md` section3 | AL-9 |
| $\hat\omega_N$ | Measurable empirical near-minimizer | `04_META_LEARNING_FORMULATION.md` section3 | AL-9 |
| $\gamma_N^{\rm opt}$ | Optimization tolerance | `04_META_LEARNING_FORMULATION.md` section3 | AL-9 |
| $\varepsilon_{\rm approx}(N)$ | Best sup-norm coefficient approximation error | `05_APPROXIMATION_THEORY.md` section1 | AL-7 |
| $\omega_{g_\mu^\star}$ | Continuity modulus of the target | `05_APPROXIMATION_THEORY.md` section3 | AL-8 |
| $L_{\rm base}$ | Base-risk coefficient Lipschitz constant | `05_APPROXIMATION_THEORY.md` section4 | CL-1 |
| $L_p^\star$ | Sole regularized coefficient-to-risk constant | `05_APPROXIMATION_THEORY.md` section4 | CL-1, CL-2 |
| $\mathcal E_\mu$ | Excess risk $R_\mu(F)-R_\mu(g_\mu^\star)$ | `06_CALIBRATION_AND_GENERALIZATION.md` section1 | AL-4 |
| $\Phi$ | Calibration map $D_V\sqrt{2t/\mu}$ | `06_CALIBRATION_AND_GENERALIZATION.md` section2 | AL-4, AL-5 |
| $\Lambda$ | Sieve-uniform parameter-loss Lipschitz constant | `06_CALIBRATION_AND_GENERALIZATION.md` section3 | AL-11 |
| $C_0$ | Absolute covering/concentration constant | `06_CALIBRATION_AND_GENERALIZATION.md` section3 | AL-11 |
| $\delta_N$ | Failure probability at level $N$ | `06_CALIBRATION_AND_GENERALIZATION.md` section3 | AL-11, CL-6 |
| $\Gamma_N$ | Uniform generalization bound | `06_CALIBRATION_AND_GENERALIZATION.md` section3 | AL-11, CL-7 |

## Active dependency chain

$$
(S\text{-IID},S\text{-CONT},S\text{-GRID},\mathcal D)
\Longrightarrow g_\mu^\star,\mathsf A,\mathcal H_N
\Longrightarrow \varepsilon_{\rm approx},\Gamma_N
\Longrightarrow \mathcal E_\mu(F_{\hat\omega_N})
\Longrightarrow d_{\mathbb M}\text{ calibration up to }2h.
$$
