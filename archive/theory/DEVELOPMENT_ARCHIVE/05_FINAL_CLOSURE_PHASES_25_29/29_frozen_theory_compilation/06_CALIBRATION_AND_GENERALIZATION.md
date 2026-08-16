# Frozen Theory: Calibration and Generalization

## 1. Risk and excess risk

For measurable $F:Z\to\Delta_m$,

$$R_\mu(F)=\mathbb E_\zeta[J_\mu(\zeta,F(\zeta))],$$

$$\mathcal E_\mu(F)=R_\mu(F)-R_\mu(g_\mu^\star)\ge0.$$

## 2. Calibration theorem

Strong convexity gives

$$\frac\mu2\|F-g_\mu^\star\|_{L^2(\mu_\zeta)}^2
\le\mathcal E_\mu(F).$$

Define

$$D_V=D_V^{\rm val}\kappa_B,$$

$$\Phi(t)=D_V\sqrt{\frac{2t}{\mu}}.$$

The fixed-grid operator transfer and Minkowski's inequality yield

$$\boxed{
\|d_{\mathbb M}(F,g_\mu^\star)\|_{L^2(\mu_\zeta)}
\le\Phi(\mathcal E_\mu(F))+2h.
}$$

Every term refers to the same operator $\mathsf A(F,z)=K(B(z)F(z))$. The additive $2h$ is a fixed design floor, not a vanishing term.

## 3. Uniform generalization

Use the block-$\ell^\infty$ parameter metric on $\Omega_N$. Multilinear interpolation is uniformly $1$-Lipschitz from this metric to coefficient Euclidean distance. The regularized per-task loss is bounded by $\bar L+\mu/2$ and has the sieve-uniform parameter Lipschitz constant

$$\Lambda=L_{\rm Lip}\kappa_B
+\mu\operatorname{diam}(\Delta_m).$$

For confidence level $\delta_N\in(0,1)$, define

$$\Gamma_N
=C_0\left(\bar L+\frac\mu2\right)
\sqrt{\frac{D_N\log(\Lambda N)+\log(1/\delta_N)}{N}},$$

where $C_0$ is the absolute constant in the retained covering/concentration bound. Under (S-IID), with probability at least $1-\delta_N$,

$$\sup_{\omega\in\Omega_N}
|\widehat R_{\mu,N}(\omega)-R_\mu(F_\omega)|
\le\Gamma_N.$$

## 4. Excess-risk decomposition

Uniform deviation, empirical near-optimality, and approximation give

$$\mathcal E_\mu(F_{\hat\omega_N})
\le2\Gamma_N+\gamma_N^{\rm opt}
+L_p^\star\varepsilon_{\rm approx}(N).$$

Substitution into calibration gives, with probability at least $1-\delta_N$,

$$\boxed{
\|d_{\mathbb M}(F_{\hat\omega_N},g_\mu^\star)\|_{L^2(\mu_\zeta)}
\le\Phi\left(
2\Gamma_N+\gamma_N^{\rm opt}
+L_p^\star\varepsilon_{\rm approx}(N)
\right)+2h.
}$$

## 5. Consistency schedule and probability statements

The high-probability result tending to confidence one requires all of

$$\operatorname{mesh}(r_N)\to0,$$

$$\frac{D_N\log(\Lambda N)}{N}\to0,$$

$$\delta_N\to0,$$

$$\frac{\log(1/\delta_N)}{N}\to0,$$

$$\gamma_N^{\rm opt}\to0.$$

Under this schedule, $\Gamma_N\to0$, $\varepsilon_{\rm approx}(N)\to0$, the bound holds with probability $1-\delta_N\to1$, and its right-hand side tends to the fixed floor $2h$.

An almost-sure eventual statement is retained only when, in addition,

$$\sum_{N=1}^\infty\delta_N<\infty.$$

Together with the logarithmic rate condition, Borel-Cantelli then gives

$$\limsup_{N\to\infty}
\|d_{\mathbb M}(F_{\hat\omega_N},g_\mu^\star)\|_{L^2(\mu_\zeta)}
\le2h
\qquad\text{almost surely}.$$

No almost-sure conclusion follows from $\delta_N\to0$ alone.

## Provenance

This chapter consolidates AL-5, AL-11, AL-12, CL-1 through CL-8, with the final probability schedule and sole constant $L_p^\star$.
