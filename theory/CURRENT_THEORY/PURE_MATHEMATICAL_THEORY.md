# Pure Mathematical Theory of Quotient-Preserving Few-Shot Prediction

## 1. Purpose and scope

This document develops a self-contained mathematical theory for few-shot prediction under four constraints:

1. only a finite support set is observed for a new task;
2. nuisance main effects must be removed without destroying interaction information;
3. task adaptation must be identifiable from the support set;
4. source-task learning must be separated from deployment shift.

The theory is independent of any particular neural architecture, dataset, or scientific application.

## 2. Probability space and task sampling

Let $(\Omega,\mathcal F,\mathbb P)$ be a probability space. A dependency component is a random object
$C$ taking values in a measurable space $\mathcal C$. The conditional laws of $\tau\mid C$ and
$(S,Q)\mid\tau,C$ are model primitives given by Markov kernels. The support and query measurement
units are disjoint, and $1\le|Q|<\infty$ almost surely.

For a predictor $h$ and bounded loss $\ell\in[0,M]$, define the component loss

$$
\bar\ell_h(C)
=\mathbb E_{\tau\mid C}
\mathbb E_{S,Q\mid\tau,C}
\left[
\frac1{|Q|}\sum_{q\in Q}\ell(h(q\mid S),Y_q)
\right].
$$

For a domain $d\in\{\mathrm{src},\mathrm{dep}\}$ with component law $\Pi_d$, define

$$
R_d(h)=\mathbb E_{C\sim\Pi_d}\bar\ell_h(C).
$$

The independent sample size is the number of independent components, not the number of rows,
queries, pairs, or resampled episodes within a component.

## 3. Information limits

### 3.1 Minimax indistinguishability

Let $V=[v_-,v_+]$ and $D_V=v_+-v_-$. Let $W$ denote all information visible to an estimator,
including the source archive, covariates, support observations, and all legal side information.

An admissible randomized estimator is one fixed Markov kernel from $W$ to its action space; its
internal random seed is independent of the data-generating mechanism.

**Theorem 1 (finite-information lower bound).** Suppose two admissible mechanisms $P_0,P_1$
satisfy $P_0^W=P_1^W$, while the same query truth equals $v_-$ under $P_0$ and $v_+$ under
$P_1$. Then, over all measurable randomized estimators,

$$
\inf_{\widehat y}\max_{i\in\{0,1\}}
\mathbb E_{P_i}|\widehat y-Y_q|
\ge\frac{D_V}{2},
$$

and

$$
\inf_{\widehat y}\max_{i\in\{0,1\}}
\mathbb E_{P_i}(\widehat y-Y_q)^2
\ge\frac{D_V^2}{4}.
$$

For ranking, restrict the action space to $\{>,<,=\}$. A correct strict action has loss $0$, the wrong
strict action has loss $1$, and a predicted tie against a strict truth has loss $1/2$. If the two
mechanisms induce opposite strict rankings and the visible-information laws are identical, then the
minimax pairwise ranking loss is at least $1/2$.

**Proof.** The estimator has the same distribution under both mechanisms. For every realized value
$t$, the triangle inequality gives
$|t-v_-|+|t-v_+|\ge D_V$, and the parallelogram identity gives
$(t-v_-)^2+(t-v_+)^2\ge D_V^2/2$. Averaging and taking the larger risk proves the two bounds.
For ranking, every strict or tied output incurs average loss at least $1/2$ across the two opposite
orders. $\square$

### 3.2 Compression risk

Let $X$ be the complete legal information and let $Z=z(X)$ be a measurable representation. Define

$$
R_X^*=\inf_{G\in L^2(\Omega,\sigma(X),\mathbb P)}\mathbb E(Y-G)^2,
\qquad
R_Z^*=\inf_{H\in L^2(\Omega,\sigma(Z),\mathbb P)}\mathbb E(Y-H)^2.
$$

**Theorem 2 (Bayes compression identity).** If $Y\in L^2$, then

$$
R_Z^*-R_X^*
=\mathbb E\operatorname{Var}
\left(\mathbb E[Y\mid X]\mid Z\right)\ge0.
$$

The gap is zero if and only if $\mathbb E[Y\mid X]$ has a $\sigma(Z)$-measurable version.

**Proof.** Orthogonal projection in $L^2$ gives the Bayes predictors
$\mathbb E[Y\mid X]$ and $\mathbb E[Y\mid Z]$. Apply the conditional variance decomposition and
use $\sigma(Z)\subseteq\sigma(X)$. $\square$

## 4. Difference operators and additive quotients

Let $P$ and $L$ be nonempty sets and let $f:P\times L\to\mathbb R$. Define

$$
D_Lf(p;l,l')=f(p,l')-f(p,l),
$$

$$
D_P\delta(p,p';l,l')
=\delta(p;l,l')-\delta(p';l,l'),
$$

and the rectangle operator

$$
\mathcal Rf(p,p';l,l')
=f(p,l')-f(p,l)-f(p',l')+f(p',l).
$$

Thus $\mathcal R=D_PD_L$.

### 4.1 Characterization of the additive kernel

**Theorem 3 (rectangle kernel).** On the complete Cartesian domain $P\times L$,

$$
\mathcal Rf\equiv0
\quad\Longleftrightarrow\quad
f(p,l)=u(p)+v(l)
$$

for some functions $u:P\to\mathbb R$ and $v:L\to\mathbb R$.

**Proof.** The reverse implication follows by cancellation. For the forward implication, fix
$(p_0,l_0)$ and set

$$
u(p)=f(p,l_0),
\qquad
v(l)=f(p_0,l)-f(p_0,l_0).
$$

The identity $\mathcal Rf(p,p_0;l_0,l)=0$ gives $f(p,l)=u(p)+v(l)$. $\square$

For a sparse domain, the same reconstruction holds on a reference star if a point $(p_0,l_0)$
exists such that, for every controlled $(p,l)$, all four cells
$(p,l),(p,l_0),(p_0,l),(p_0,l_0)$ are observed and the corresponding reference rectangle vanishes.

### 4.2 Antisymmetry and integrability

The ligand difference is antisymmetric under $l\leftrightarrow l'$. The rectangle is separately
antisymmetric under $p\leftrightarrow p'$ and $l\leftrightarrow l'$.

Let $G=(L,E)$ be a connected graph and let $\delta:E^{\mathrm{or}}\to\mathbb R$ be an oriented
edge field.

**Theorem 4 (discrete potential criterion).** There exists $s:L\to\mathbb R$ such that

$$
\delta(l,l')=s(l')-s(l)
$$

if and only if $\delta$ is antisymmetric and its circulation along every cycle is zero. The potential
is unique up to an additive constant.

**Proof.** A potential difference is antisymmetric and telescopes along cycles. Conversely, fix a
root and define $s(l)$ as the path integral from the root to $l$. Zero circulation makes the value
path independent. Connectedness gives uniqueness up to the root value. $\square$

Antisymmetry alone therefore does not imply integrability, and neither antisymmetry nor integrability
removes a second-coordinate-only difference. The rectangle quotient is required for that purpose.

### 4.3 Deterministic quotient reconstruction

**Theorem 5 (reconstruction modulo additive gauge).** Let $f,g:P\times L\to\mathbb R$, fix
$(p_0,l_0)$, and write $e=f-g$. Define

$$
u(p)=e(p,l_0),
\qquad
v(l)=e(p_0,l)-e(p_0,l_0).
$$

Then

$$
e(p,l)-u(p)-v(l)=\mathcal Re(p,p_0;l_0,l).
$$

Consequently, if all reference rectangles satisfy
$|\mathcal R(f-g)|\le\varepsilon_R$, then $f-g$ is within $\varepsilon_R$ of an additive function
pointwise on the reference star.

An average rectangle error does not imply a pointwise scalar error without an additional norm
equivalence or spectral-gap condition, reference-star coverage, and an anchor for the additive gauge.

## 5. Affine task sections

Let $\phi:X\to\mathbb R^r$ and $b:X\to\mathbb R$ be shared functions, and fix
$R\in[0,\infty)$. A task-specific member is

$$
f_c(x)=b(x)+\phi(x)^\top c,
\qquad
\|c\|_2\le R.
$$

For support points $x_1,\ldots,x_k$ with point observations $y_i$, define

$$
r_S=(y_i-b(x_i))_{i=1}^k,
\qquad
\Phi_S=(\phi(x_i)^\top)_{i=1}^k.
$$

### 5.1 Exact identifiability

For a feasible coefficient $c^\dagger$, the support-equivalent coefficients are

$$
\mathcal F_S(c^\dagger)
=\{c^\dagger+n:n\in\ker\Phi_S, \|c^\dagger+n\|_2\le R\}.
$$

A query $q$ is identifiable if and only if $\phi(q)^\top n=0$ for every feasible null direction.
In the nondegenerate interior case $\|c^\dagger\|_2<R$, this is equivalent to

$$
\phi(q)\in\operatorname{row}(\Phi_S).
$$

Complete coefficient recovery occurs if and only if the feasible fiber is a singleton. In the
nondegenerate interior case, this requires $\operatorname{rank}(\Phi_S)=r$.

### 5.2 Centered adaptation

Suppose an unknown scalar level is a nuisance parameter. For $k\ge1$, define

$$
H_k=I_k-\frac1k\mathbf1\mathbf1^\top.
$$

The centered ridge state is

$$
\widehat c(S)
=\arg\min_{c\in\mathbb R^r}
\left\{
\frac1k\|H_k(r_S-\Phi_Sc)\|_2^2+\lambda\|c\|_2^2
\right\},
\qquad\lambda>0.
$$

Let $A=H_k\Phi_S$ and $d=H_kr_S$. Then

$$
\widehat c(S)
=\left(\frac1kA^\top A+\lambda I\right)^{-1}
\frac1kA^\top d.
$$

**Proposition 1 (existence, stability, and information rank).** The centered objective is
$2\lambda$-strongly convex, so the minimizer exists and is unique. Moreover,

$$
\operatorname{rank}(H_k\Phi_S)
\le\min\{k-1,r\}.
$$

If $A$ and $d$ are $C^1$ functions of shared parameters, then $\widehat c$ is $C^1$ and

$$
\frac{d\widehat c}{d\theta}
=-[\nabla_{cc}^2\mathcal L_{\mathrm{in}}]^{-1}
\nabla_{\theta c}^2\mathcal L_{\mathrm{in}},
$$

with

$$
\|[\nabla_{cc}^2\mathcal L_{\mathrm{in}}]^{-1}\|_{\mathrm{op}}
\le\frac1{2\lambda}.
$$

The zero-support branch is defined separately by $S=\varnothing$ and $c=0$; $H_0$ is not used.
The rank ceiling concerns only the centered SAR state. It is not an endpoint-risk optimality
statement and does not justify discarding a support-free query predictor or the task-level channel.

### 5.3 Random-intercept level calibration

Let support residuals satisfy the conditional Gaussian random-intercept model

$$
A\sim\mathcal N(m,\tau^2),\qquad r_i=A+\varepsilon_i,
\qquad \varepsilon_i\overset{\mathrm{iid}}\sim\mathcal N(0,\sigma^2),
$$

where $\sigma^2>0$ and $\tau^2\ge0$. Then, for $k\ge1$,

$$
\mathbb E[A\mid r_1,\ldots,r_k]
=m+\frac{k\tau^2}{\sigma^2+k\tau^2}(\overline r-m).
$$

**Proposition 2 (Bayes level shrinkage).** Under squared loss and the stated model, the posterior
mean above is the Bayes estimator of the task level. Relative to the raw residual mean, its Bayes
risk reduction is

$$
\frac{\sigma^2}{k}-\frac{\tau^2\sigma^2}{\sigma^2+k\tau^2}
=\frac{\sigma^4}{k(\sigma^2+k\tau^2)}\ge0.
$$

When the hyperparameters are estimated from source tasks, this becomes a plug-in empirical-Bayes
rule; the finite-sample inequality is not automatically inherited under misspecification.

### 5.4 Interval observations and exact sections

Let each support observation be a nonempty closed interval $O_i$. Define

$$
\mathcal C_S
=\{c:\|c\|_2\le R, b(x_i)+\phi(x_i)^\top c\in O_i\ \forall i\}.
$$

If $\mathcal C_S\ne\varnothing$, it is compact and convex. For a query $q$, define

$$
I_q(S)=
\left[
\min_{c\in\mathcal C_S}f_c(q),
\max_{c\in\mathcal C_S}f_c(q)
\right].
$$

The midpoint $m_q^{\mathrm{sec}}$ is the conditional Chebyshev predictor under absolute loss, and
the half-width $r_q^{\mathrm{sec}}$ is the exact conditional section radius.

For feasible noiseless point support, let $c^\dagger$ be the minimum-norm solution of
$\Phi_Sc=r_S$, with $\|c^\dagger\|_2\le R$. Then

$$
r_q^{\mathrm{sec}}
=\sqrt{R^2-\|c^\dagger\|_2^2}
\|P_{\ker\Phi_S}\phi(q)\|_2.
$$

The formula includes the boundary case $\|c^\dagger\|_2=R$, where the radius is zero. If the
coefficient set is empty, the section is undefined and the procedure must abstain.

## 6. Certificates and selective decisions

Let $\widehat f_q$ be a trained point predictor, which need not equal $m_q^{\mathrm{sec}}$. Define

$$
r_q^{\mathrm{center}}=|\widehat f_q-m_q^{\mathrm{sec}}|.
$$

Let $f_q^*$ be the true latent query value on the same scale as the section. Assume that, on one
joint event $\mathcal E_q$, the representation, pointwise deployment-shift, and observation bounds
satisfy

$$
\operatorname{dist}(f_q^*,I_q(S))
\le r_q^{\mathrm{repr}}+r_q^{\mathrm{trans}}+r_q^{\mathrm{obs}}.
$$

Then a valid total radius around $\widehat f_q$ is

$$
r_q^{\mathrm{tot}}
\ge r_q^{\mathrm{center}}+r_q^{\mathrm{sec}}
+r_q^{\mathrm{repr}}+r_q^{\mathrm{trans}}+r_q^{\mathrm{obs}}.
$$

Indeed, on $\mathcal E_q$,

$$
|\widehat f_q-f_q^*|\le r_q^{\mathrm{tot}}.
$$

If the component bounds are only marginal, a simultaneous calibration rule or an explicit union
bound is required.

For two queries $q,q'$ on the same truth scale, if both certificates hold on one joint event and

$$
|\widehat f_q-\widehat f_{q'}|>r_q^{\mathrm{tot}}+r_{q'}^{\mathrm{tot}},
$$

then the predicted strict order is correct on that event. Coverage of the acceptance event and
conditional ranking error must be controlled separately.

## 7. Component-level learning

Let $\mathcal H$ be a fixed pointwise measurable and separable hypothesis class. Assume
$\bar\ell_h(C)\in[0,M]$. Given IID source components $C_1,\ldots,C_N$, define

$$
\widehat R_{\mathrm{src}}(h)
=\frac1N\sum_{i=1}^N\bar\ell_h(C_i),
$$

and the expected Rademacher complexity

$$
\mathfrak R_N(\bar\ell\circ\mathcal H)
=\mathbb E_{C,\sigma}
\sup_{h\in\mathcal H}
\frac1N\sum_{i=1}^N\sigma_i\bar\ell_h(C_i),
$$

where the $\sigma_i$ are IID Rademacher signs independent of the components. Let the training
algorithm output $\widehat h$ with empirical regret

$$
\gamma
=\widehat R_{\mathrm{src}}(\widehat h)
-\inf_{h\in\mathcal H}\widehat R_{\mathrm{src}}(h)\ge0.
$$

Define

$$
\Delta_{\mathcal H}
=\sup_{h\in\mathcal H}|R_{\mathrm{dep}}(h)-R_{\mathrm{src}}(h)|.
$$

**Theorem 6 (component generalization with shift).** With probability at least $1-\delta$,

$$
R_{\mathrm{dep}}(\widehat h)
-\inf_{h\in\mathcal H}R_{\mathrm{dep}}(h)
\le
4\mathfrak R_N(\bar\ell\circ\mathcal H)
+2M\sqrt{\frac{\log(2/\delta)}{2N}}
+\gamma+2\Delta_{\mathcal H}.
$$

**Proof.** The standard expected-Rademacher uniform deviation bound controls source excess risk by
the first three terms. Apply
$|R_{\mathrm{dep}}(h)-R_{\mathrm{src}}(h)|\le\Delta_{\mathcal H}$ to $\widehat h$ and to an
$\varepsilon$-optimal deployment comparator, then let $\varepsilon\downarrow0$. $\square$

## 8. Transport decomposition

Let $(\mathcal D,d)$ be a Polish metric space with its Borel sigma-field, and let
$D(C)\in\mathcal D$ be a measurable descriptor. Write $\Pi_s^D$ and $\Pi_t^D$ for the source and
deployment Borel probability pushforwards and assume finite first moments. Fix versions

$$
m_h^a(z)=\mathbb E_a[\bar\ell_h(C)\mid D=z],
\qquad a\in\{s,t\}.
$$

Assume every $m_h^s$ has a specified, everywhere-defined $L$-Lipschitz version on
$\operatorname{supp}(\Pi_s^D)\cup\operatorname{supp}(\Pi_t^D)$, or use a fixed Lipschitz extension.
Define

$$
\Delta_{\mathrm{mech}}
=\sup_{h\in\mathcal H}
\mathbb E_{D\sim\Pi_t^D}|m_h^t(D)-m_h^s(D)|.
$$

**Proposition 2 (transport and mechanism decomposition).**

$$
\Delta_{\mathcal H}
\le L W_1(\Pi_s^D,\Pi_t^D;d)+\Delta_{\mathrm{mech}}.
$$

**Proof.** Add and subtract the integral of the specified source conditional-risk version under the
deployment descriptor law. Bound the marginal transport term by Kantorovich-Rubinstein duality and
the conditional-mechanism term by its definition. $\square$

Let $R_{\mathrm{dep}}^*$ be the Bayes risk over all legal deployment predictors. Assume every member
of $\mathcal H$ belongs to that comparison class, and define

$$
\varepsilon_{\mathrm{app}}^{\mathrm{dep}}
=\inf_{h\in\mathcal H}R_{\mathrm{dep}}(h)-R_{\mathrm{dep}}^*\ge0.
$$

Combining Theorem 6 and Proposition 2 yields

$$
R_{\mathrm{dep}}(\widehat h)-R_{\mathrm{dep}}^*
\le U_{\mathrm{gen}}
+2\left[LW_1(\Pi_s^D,\Pi_t^D;d)+\Delta_{\mathrm{mech}}\right]
+\varepsilon_{\mathrm{app}}^{\mathrm{dep}},
$$

where

$$
U_{\mathrm{gen}}
=4\mathfrak R_N(\bar\ell\circ\mathcal H)
+2M\sqrt{\frac{\log(2/\delta)}{2N}}+\gamma.
$$

This result does not imply that the transport, mechanism, approximation, or complexity terms are
small. Lipschitz continuity of a predictor with respect to its inputs does not by itself imply
Lipschitz continuity of the conditional component risk.

Theorem 6 and Proposition 2 are expectation-risk results. They do not by themselves construct the
pointwise quantity $r_q^{\mathrm{trans}}$ required by a query certificate.

## 9. Logical dependency of the theory

The valid implication chain is

```text
visible information
    -> representation sufficiency
    -> quotient coverage
    -> scalar integrability
    -> finite-support section
    -> trained point prediction
    -> joint certificate
    -> component generalization
    -> explicit transport and approximation terms
```

No arrow may be reversed without an additional theorem. In particular:

- predictive performance does not identify a mechanism;
- antisymmetry does not identify cross-variable interaction;
- average quotient accuracy does not imply pointwise scalar accuracy;
- a ridge point is not automatically a section midpoint;
- IID episodes do not replace independent components;
- representation distance does not prove conditional-mechanism invariance;
- a marginal uncertainty object does not automatically imply ranking validity.
