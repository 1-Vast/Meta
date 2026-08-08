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
