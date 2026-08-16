# Frozen Theory: Meta-Learning Formulation

## 1. Task law and support/query structure

On one common probability space, let

$$T_i=(S_i,Q_i,Y_i)\sim P_T,\qquad i=1,2,\ldots,$$

be IID observable tasks. Here $S_i$ is the finite support set, $Q_i$ is the query, and $Y_i\in V$ is the identified point-valued affinity target available during meta-training. The task statistic is

$$\zeta_i=z(S_i,Q_i,\gamma)\in Z.$$

At deployment, the learned map is evaluated at the current task's support and query through the same statistic. Thus the formulation is support/query conditioned; it does not take a hidden task identity or query label as an inference input.

## 2. Parameter space and realization

For sieve level $N$, choose a finite mesh of $Z$ at resolution $r_N$, with node set $\mathcal N_N$ and node count $\nu_N$. Write $\operatorname{mesh}(r_N)$ for the maximum $d_Z$-diameter of a mesh cell.

The parameter space is

$$\Omega_N=(\Delta_m)^{\mathcal N_N}
=\{\omega=(\omega_\nu)_{\nu\in\mathcal N_N}:\omega_\nu\in\Delta_m\},$$

a compact convex subset of $(\mathbb R^{m+1})^{\nu_N}$ with ambient dimension

$$D_N=(m+1)\nu_N.$$

Let $\phi_\nu(z)$ be the nonnegative piecewise-multilinear basis functions, satisfying $\sum_\nu\phi_\nu(z)=1$. Define

$$G_N(\omega,z)=\sum_{\nu\in\mathcal N_N}\phi_\nu(z)\omega_\nu,$$

$$F_\omega(z)=G_N(\omega,z),\qquad
\mathcal H_N=\{F_\omega:\omega\in\Omega_N\}.$$

Because $F_\omega(z)$ is a convex combination of simplex points, every hypothesis is $\Delta_m$-valued. Its operator output is always

$$\mathsf A(F_\omega,z)=K(B(z)F_\omega(z)).$$

## 3. Empirical and population objectives

The empirical regularized risk is

$$\widehat R_{\mu,N}(\omega)
=\frac1N\sum_{i=1}^N\left[
L(B(\zeta_i)F_\omega(\zeta_i),Y_i)
+\frac\mu2\|F_\omega(\zeta_i)\|^2
\right].$$

Its expectation under (S-IID) is the population regularized risk

$$R_\mu(F_\omega)=\mathbb E_\zeta[J_\mu(\zeta,F_\omega(\zeta))].$$

An exact empirical minimizer exists by continuity on compact $\Omega_N$. The retained statement also allows a measurable approximate minimizer $\hat\omega_N$ satisfying

$$\widehat R_{\mu,N}(\hat\omega_N)
\le\inf_{\omega\in\Omega_N}\widehat R_{\mu,N}(\omega)
+\gamma_N^{\rm opt},$$

where $\gamma_N^{\rm opt}\ge0$ is the declared optimization tolerance.

The learned coefficient map is $F_{\hat\omega_N}$, and the learned operator on a current task is

$$\mathsf A(F_{\hat\omega_N},z(S,Q,\gamma)).$$

## Provenance

This chapter consolidates AL-6, AL-9, AL-10, MR-1, and the retained fixed-deployment task typing.
