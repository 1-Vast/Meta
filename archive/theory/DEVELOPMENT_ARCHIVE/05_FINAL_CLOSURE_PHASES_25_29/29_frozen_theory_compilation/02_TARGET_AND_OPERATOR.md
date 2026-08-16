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
