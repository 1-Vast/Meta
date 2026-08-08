# Frozen Theory: Foundations

## 1. Problem

For each task, an observable support set `S`, a query `Q`, and a declared specification `gamma` determine a statistic

$$z=z(S,Q,\gamma)\in Z.$$

The learned object is a measurable coefficient map $F:Z\to\Delta_m$. Its output is decoded into a valid set of probability laws for a continuous scalar affinity. Historical tasks are used only to estimate $F$; the deployment state is frozen.

## 2. Fixed deployment

The complete deployment declaration is

$$\mathcal D=(z_H^0,B(\cdot),\Delta_m,\mu,h).$$

- $z_H^0$ is the frozen historical/deployment state.
- $(Z,d_Z)$ is a compact metric statistic domain, represented by a finite union of compact cubes.
- $\kappa:Z\to C_\kappa$ is a measurable context map with finite codomain.
- $V=[a_{\min},a_{\max}]\subset\mathbb R$ is the compact affinity range, with diameter $D_V^{\rm val}=a_{\max}-a_{\min}$.
- $h>0$ is the fixed output-grid mesh. It is not refined by any theorem in this package.
- $\Delta_m=\{p\in\mathbb R^{m+1}:p_k\ge0,\ \sum_{k=0}^m p_k=1\}$ is the compact convex coefficient simplex with Euclidean norm.
- $\operatorname{diam}(\Delta_m)=\sup_{p,q\in\Delta_m}\|p-q\|$ is its Euclidean diameter.
- $\mu>0$ is the fixed ridge modulus.

## 3. Band assembly

Let $\mathbb B$ be the compact convex polytope of valid lower/upper CDF-band vectors on the fixed mesh. The deployment determines the matrix-valued rule

$$B(z)=[\beta_0(z)\mid\beta_1\mid\cdots\mid\beta_m],\qquad
\beta_0(z)=b^{\rm pop}_{\kappa(z)},$$

where every column lies in $\mathbb B$ and $\beta_1,\ldots,\beta_m$ are fixed anchors. For each fixed $z$, assembly is linear:

$$p\longmapsto B(z)p=\sum_{k=0}^m p_k\beta_k(z)\in\mathbb B.$$

Define the finite assembly norm

$$\kappa_B=\sup_{z\in Z}\|B(z)\|_{\rm op}<\infty,$$

where the operator norm maps Euclidean coefficient distance to the band sup norm.

## 4. Data and loss

A task is

$$T=(S,Q,Y),\qquad Y\in V,$$

where $Y$ is the observable identified point target. The loss

$$L:\mathbb B\times V\to[0,\infty)$$

is convex in its band argument, bounded by $\bar L$, and $L_{\rm Lip}$-Lipschitz in the band sup norm, uniformly in $Y$.

For $\zeta=z(S,Q,\gamma)$, let $\mu_\zeta$ denote the law of $\zeta$ on $Z$ and define the conditional base risk

$$L_0(z,\beta)=\mathbb E[L(\beta,Y)\mid\zeta=z].$$

## 5. Assumptions

The frozen theory uses exactly these assumptions:

1. **(S-IID)** Meta-training tasks and the current task are IID draws from one observable task law $P_T$.
2. **(S-CONT)** $L_0$ has an everywhere-defined version satisfying, for a declared modulus $\varpi_\ell$,
   $$\sup_{\beta\in\mathbb B}|L_0(z,\beta)-L_0(z',\beta)|\le\varpi_\ell(d_Z(z,z')),
   \qquad \varpi_\ell(t)\to0\ \text{as }t\downarrow0.$$
3. **(S-GRID)** The Route-B affinity grid and mesh $h$ are fixed components of $\mathcal D$.

No conditional-IID branch, distribution-shift theorem, ranking assumption, continuum-mesh assumption, or varying-$z_H$ assumption is retained.

## Provenance

This chapter consolidates the validated definitions AL-2, AL-3, AL-9, MR-1, and CL-1 without changing their mathematical content.


# Frozen Theory: Target and Operator

## 1. Probability-law output space

Let $\Delta(V)$ be the Borel probability laws on $V$, equipped with Wasserstein-1 distance

$$W_1(P,P')=\int_V|F_P(v)-F_{P'}(v)|\,dv.$$

For a valid CDF band $\beta\in\mathbb B$, let $K(\beta)\subset\Delta(V)$ be the nonempty compact, convex, $W_1$-closed class of laws satisfying that band. Thus

$$K:\mathbb B\to\mathcal K(\Delta(V)),$$

where $\mathcal K(\Delta(V))$ denotes nonempty compact subsets of $(\Delta(V),W_1)$.

## 2. The sole operator

For every measurable $F:Z\to\Delta_m$, the retained output operator is exactly

$$\boxed{\mathsf A(F,z)=K(B(z)F(z)).}$$

No support intersection, identified-set restriction, confidence coordinate, rung coordinate, or alternate output is part of this frozen operator.

The pointwise operator metric is

$$d_{\mathbb M}(F,G)(z)
=d_H^{W_1}(\mathsf A(F,z),\mathsf A(G,z)),$$

where $d_H^{W_1}$ is Hausdorff distance induced by $W_1$. Its population norm is

$$\|d_{\mathbb M}(F,G)\|_{L^2(\mu_\zeta)}
=\left(\mathbb E_{\zeta}\,d_{\mathbb M}(F,G)(\zeta)^2\right)^{1/2}.$$

The fixed-grid class stability bound is

$$d_H^{W_1}(K(\beta),K(\beta'))
\le D_V^{\rm val}\|\beta-\beta'\|_{\mathbb B}+2h.$$

Combining it with assembly gives

$$d_{\mathbb M}(F,G)(z)
\le D_V\|F(z)-G(z)\|+2h,
\qquad D_V:=D_V^{\rm val}\kappa_B.$$

## 3. The sole target

Define the operative local risk

$$J_\mu(z,p)=L_0(z,B(z)p)+\frac\mu2\|p\|^2.$$

The unique learning target is

$$\boxed{g_\mu^\star(z)=\operatorname*{arg\,min}_{p\in\Delta_m}J_\mu(z,p).}$$

The positive ridge is part of the target-defining objective. It is also present in the empirical and population risks. The theory makes no claim that $g_\mu^\star$ equals an unregularized Bayes target.

The target operator is the same retained operator evaluated at this coefficient map:

$$\mathsf A(g_\mu^\star,z)=K(B(z)g_\mu^\star(z)).$$

## 4. Positive-ridge justification

For each $z$, $p\mapsto L_0(z,B(z)p)$ is convex because $B(z)p$ is linear in $p$ and the band loss is convex. Adding $\frac\mu2\|p\|^2$ makes $J_\mu(z,\cdot)$ at least $\mu$-strongly convex. The negative quadratic is not an equivalent notation: it would remove this guarantee and is not retained.

## Provenance

This chapter consolidates AL-1, AL-4, AL-5, MR-1, MR-2, and the validated fixed-grid stability statement.


# Frozen Theory: Strong Convexity and Regularity

## 1. Existence and uniqueness

For every $z\in Z$, $J_\mu(z,\cdot)$ is continuous and at least $\mu$-strongly convex on the compact convex simplex $\Delta_m$. Therefore it attains a unique minimizer $g_\mu^\star(z)$.

Assumption (S-CONT) supplies an everywhere-defined continuous version of $L_0$. Consequently $J_\mu$ is measurable in $z$, continuous in $p$, and the unique argmin map $g_\mu^\star:Z\to\Delta_m$ is measurable.

## 2. Square-root continuity modulus

Let

$$p=g_\mu^\star(z),\qquad q=g_\mu^\star(z').$$

Strong convexity and optimality at $z$ give

$$J_\mu(z,q)-J_\mu(z,p)\ge\frac\mu2\|q-p\|^2.$$

Strong convexity and optimality at $z'$ give

$$J_\mu(z',p)-J_\mu(z',q)\ge\frac\mu2\|p-q\|^2.$$

Adding these inequalities yields

$$\mu\|p-q\|^2
\le [J_\mu(z,q)-J_\mu(z',q)]
   +[J_\mu(z',p)-J_\mu(z,p)].$$

The ridge is independent of $z$, so (S-CONT) bounds each bracket by $\varpi_\ell(d_Z(z,z'))$. Hence

$$\boxed{
\|g_\mu^\star(z)-g_\mu^\star(z')\|
\le\sqrt{\frac{2\varpi_\ell(d_Z(z,z'))}{\mu}}.
}$$

Thus $g_\mu^\star$ is uniformly continuous on compact $Z$. If $\varpi_\ell(t)=O(t)$, the proved conclusion is generally Holder-$1/2$, not Lipschitz. No gradient-regularity strengthening is retained.

## 3. Regularized Bayes optimality

For any measurable $F:Z\to\Delta_m$, define

$$R_\mu(F)=\mathbb E_\zeta[J_\mu(\zeta,F(\zeta))].$$

Pointwise optimality gives

$$R_\mu(F)\ge R_\mu(g_\mu^\star),$$

with equality only when $F=g_\mu^\star$ $\mu_\zeta$-almost everywhere. This is optimality for the regularized risk only.

## Provenance

This chapter consolidates TI-2, TI-5, TI-6, AL-4, AL-8, and MR-1 without adding a stronger regularity claim.


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


# Frozen Theory: Approximation

## 1. Approximation error

For the hypothesis class $\mathcal H_N$, define

$$\varepsilon_{\rm approx}(N)
=\inf_{F\in\mathcal H_N}\sup_{z\in Z}
\|F(z)-g_\mu^\star(z)\|.$$

## 2. Interpolation witness

Assign the target value to every mesh node:

$$\omega_\nu^\star=g_\mu^\star(\nu),\qquad \nu\in\mathcal N_N.$$

Because $g_\mu^\star(\nu)\in\Delta_m$, the parameter $\omega^\star=(\omega_\nu^\star)_\nu$ lies in $\Omega_N$. Therefore the multilinear witness

$$F_{\omega^\star}(z)=\sum_\nu\phi_\nu(z)g_\mu^\star(\nu)$$

belongs to $\mathcal H_N$.

## 3. Density bound

Write

$$\omega_{g_\mu^\star}(t)
=\sup\{\|g_\mu^\star(z)-g_\mu^\star(z')\|:d_Z(z,z')\le t\}$$

for the target's continuity modulus. Multilinear interpolation and this modulus give

$$\begin{aligned}
\varepsilon_{\rm approx}(N)
&\le\sup_{z\in Z}\|F_{\omega^\star}(z)-g_\mu^\star(z)\|\\
&\le\omega_{g_\mu^\star}(\operatorname{mesh}(r_N))\\
&\le\sqrt{\frac{2\varpi_\ell(\operatorname{mesh}(r_N))}{\mu}}.
\end{aligned}$$

Consequently,

$$\operatorname{mesh}(r_N)\to0
\quad\Longrightarrow\quad
\varepsilon_{\rm approx}(N)\to0.$$

Only this forward implication is retained. Approximation follows because the witness is a member of $\mathcal H_N$ and the node mesh refines; parameter-dimension growth by itself is not an approximation theorem.

## 4. Approximation contribution to risk

Define the base coefficient Lipschitz constant

$$L_{\rm base}
=\sup_{z\in Z}\sup_{p_1\ne p_2}
\frac{|L_0(z,B(z)p_1)-L_0(z,B(z)p_2)|}
{\|p_1-p_2\|}
\le L_{\rm Lip}\kappa_B.$$

The sole regularized coefficient-to-risk constant is

$$\boxed{L_p^\star=L_{\rm base}+\mu\operatorname{diam}(\Delta_m).}$$

It satisfies

$$|J_\mu(z,p_1)-J_\mu(z,p_2)|
\le L_p^\star\|p_1-p_2\|,$$

and hence

$$|R_\mu(F)-R_\mu(G)|
\le L_p^\star\|F-G\|_{L^\infty}.$$

Here $\|F-G\|_{L^\infty}=\sup_{z\in Z}\|F(z)-G(z)\|$.

Therefore

$$\inf_{F\in\mathcal H_N}R_\mu(F)-R_\mu(g_\mu^\star)
\le L_p^\star\varepsilon_{\rm approx}(N).$$

No alternate coefficient-to-risk constant is retained.

## Provenance

This chapter consolidates AL-7, AL-8, CL-1, CL-2, CL-3, and the forward-only correction recorded in the Phase-28.1 stopping criterion.


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


# Frozen Theory: Scope and Limitations

## Supported scope

The frozen theory supports exactly:

1. **One fixed deployment**
   $$\mathcal D=(z_H^0,B(\cdot),\Delta_m,\mu,h).$$
   All population bands, anchors, ridge modulus, and output mesh are fixed.

2. **Support-conditioned meta-learning**
   Historical IID tasks train $F_{\hat\omega_N}$. At inference, the current support and query enter through $z(S,Q,\gamma)$.

3. **Continuous point-valued affinity regression**
   The supervised target is $Y\in V\subset\mathbb R$. The emitted object is a fixed-grid CDF-band class of probability laws over $V$.

4. **Fixed-resolution guarantees**
   Approximation resolution for the coefficient map may refine, but the output-grid mesh $h$ remains fixed. Calibration and consistency retain the design floor $2h$.

## Not claimed

The frozen theory provides no theorem for:

- pairwise, listwise, or metric ranking;
- coherent joint-order learning;
- derivation of ranking from affinity regression;
- continuum output-mesh refinement or a zero-mesh target;
- convergence between mesh-indexed targets;
- a model conditioned on varying $z_H$;
- transport to an undeclared task distribution;
- conditional-IID fibers or missing-fiber generalization;
- optimization efficiency or an architecture-specific training algorithm;
- equality between the regularized target and an unregularized Bayes target;
- a support-intersected output class.

## Interpretation of the result

The theorem guarantees that empirical learning approaches the single regularized target operator in the declared operator metric up to the fixed $2h$ floor, with high probability under the full schedule in `06_CALIBRATION_AND_GENERALIZATION.md`. It does not claim exact zero-error recovery at fixed output resolution.

## Final status

Within these explicit limitations, the retained mathematical theory is frozen. The limitations are part of the contract and must not be removed by interpretation.


# Frozen Theory: Rejected Formulations

This chapter records failed formulations for provenance only. None is part of the active definitions or theorems.

## 1. Negative ridge

**Rejected form**

$$L_0(z,B(z)p)-\frac\mu2\|p\|^2.$$

**Reason rejected:** subtracting the quadratic does not supply $\mu$-strong convexity. When the base risk is flat in a direction, the negative quadratic is strictly concave there. Existence may survive by compactness, but uniqueness, continuous argmin control, coefficient calibration, and the retained consistency chain do not follow.

**Retained form:** the positive ridge in $J_\mu$.

## 2. Old coordinate target

**Rejected form:** applying an Euclidean ridge to the earlier coordinates $(\lambda,w)$ and identifying its minimizer with the barycentric-$p$ minimizer.

**Reason rejected:** Euclidean ridge is not invariant under $p_0=1-\lambda$, $p_j=\lambda w_j$. Equal assembled-band images do not imply equal regularized objectives or minimizers.

**Retained form:** one target defined directly in $p\in\Delta_m$.

## 3. Support-intersection output

**Rejected form:** emitting the class $K(B(z)F(z))$ after an additional support/identified-set intersection while calibrating the unrestricted class.

**Reason rejected:** the declared output and calibrated object were different. Intersection need not preserve nonemptiness or be nonexpansive in Hausdorff distance.

**Retained form:** $\mathsf A(F,z)=K(B(z)F(z))$ only.

## 4. Continuum mesh convergence

**Rejected form:** inferring convergence of mesh-dependent risk minimizers to a continuum target from a same-band discretization estimate.

**Reason rejected:** a bound comparing discretized and continuum classes for the same band does not compare different minimizers. That would require a typed objective family and an argmin-stability theorem not present in the retained scope.

**Retained form:** fixed output mesh $h$ and additive floor $2h$.

## 5. Approximation from dimension growth alone

**Rejected form:** concluding $\varepsilon_{\rm approx}\to0$ merely because parameter dimension tends to infinity.

**Reason rejected:** a growing family may still contain only constant or otherwise nondense maps.

**Retained form:** the hypothesis class is the multilinear witness family itself, and density follows from node refinement plus target continuity.

## 6. Canonical/risk-optimal target switching

**Rejected form:** using a computable canonical operator for approximation or imitation while calling the task-risk minimizer the same target.

**Reason rejected:** the two objects need not coincide, and their supervision and guarantees are different.

**Retained form:** one positive-ridge regularized risk target $g_\mu^\star$ throughout.

## 7. Weak probability schedule

**Rejected forms:** using only $\delta_N\to0$, or only $\log(1/\delta_N)/N\to0$, and claiming both confidence tending to one and $\Gamma_N\to0$; claiming an almost-sure result without summability.

**Reason rejected:** the two high-probability conditions have independent roles, and Borel-Cantelli requires summable failure probabilities.

**Retained form:** both high-probability conditions, plus summability for the optional almost-sure conclusion.


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
