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
