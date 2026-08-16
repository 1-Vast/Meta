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
