# Convex Parameterization Repair (Item 2)

> **Status:** Phase-22.1, 2026-08-03. Repairs the audit's finding 3: the coefficient $c=(\lambda,w)$ made the assembly **bilinear** ($\lambda w_j$ terms), so convexity of the band loss did not transfer to $c$ and adjoining $\tfrac\mu2\|c\|^2$ did not yield $\mu$-strong convexity (audit's $(\lambda w-a)^2$ indefinite-Hessian witness). The repair reparameterizes to make the assembly **linear**, after which the ridge is genuinely strongly convex with modulus exactly $\mu$. New results **RP-1–RP-3**, tagged **[proved] / [declared]**.

---

## 1. The linear reparameterization

Fix $z_H$ (deployment constant, Item 1), so the basis vectors $\beta_0=b^{\mathrm{pop}}_{\kappa(z)},\beta_1=b_1,\dots,\beta_m=b_m\in\mathbb B$ are **fixed** (the anchors are SKEL constants; $b^{\mathrm{pop}}$ is a fixed function of the fixed $z_H$ and the observable context $\kappa(z)$).

**Definition RP-1 (barycentric coefficient). [declared]**
Replace $c=(\lambda,w)\in[0,1]\times\Delta_{m-1}$ by the **direct barycentric weight**
$$p=(p_0,p_1,\dots,p_m)\in\Delta_m\ \ (\text{the }m\text{-simplex in }\mathbb R^{m+1}),\qquad\text{assembly}\quad \mathsf{asm}(p)=\sum_{k=0}^m p_k\beta_k\ =\ Bp,$$
where $B=[\beta_0\,\cdots\,\beta_m]$ is the fixed band-matrix. The old map $(\lambda,w)\mapsto(1-\lambda)\beta_0+\lambda\sum_jw_j\beta_j$ is exactly $p_0=1-\lambda,\ p_j=\lambda w_j$ — so **$\Delta_m$ is the image of the old coefficients, with no loss** (surjective onto the same $\operatorname{conv}\{\beta_k\}=\mathbb B$), and $p\mapsto Bp$ is **linear**, eliminating the bilinearity at its source. The coefficient space is now $C=\Delta_m$: compact convex.

## 2. Genuine strong convexity

**Theorem RP-2 (the operative risk is $\mu$-strongly convex in $p$ — with modulus exactly $\mu$). [proved]**
Under (A-LOSS) [band-score $L(\cdot,a)$ convex in the band vector] and (A-STAT):
$$\ell_0(z,p)=\mathbb E\big[L(Bp,A_T)\mid\zeta=z\big]\ \text{is convex in }p;\qquad \ell(z,p)=\ell_0(z,p)+\tfrac\mu2\|p\|^2\ \text{is }\mu\text{-strongly convex in }p.$$
*Proof.* $p\mapsto Bp$ is linear; $L(\cdot,a)$ is convex in the band vector (interval score: width $u-l$ affine, violation $\mathrm{dist}(y,[l,u])$ convex); composition of a convex function with a linear map is convex, and expectation preserves it — so $\ell_0(z,\cdot)$ is convex (finite, bounded loss). Adding $\tfrac\mu2\|p\|^2$ (a $\mu$-strongly convex function) to a convex function yields a $\mu$-strongly convex function, **and the modulus is exactly $\mu$** because the base is genuinely convex (nonnegative curvature everywhere) — the ridge no longer has to dominate negative curvature, which was precisely the audit's objection to the bilinear case. The $(\lambda w-a)^2$ counterexample is dissolved: in barycentric coordinates that term is $(p_j-a)^2$-type, convex in $p$. $\square$

**Corollary RP-3 (well-definedness, now with the right hypotheses — repairs the PT-2 measurability overstatement). [proved]**
Under (A-STAT, A-LOSS, A-SC, **A-CONT**): (A-CONT) selects a continuous — hence measurable and **everywhere-defined** — version of $z\mapsto\ell_0(z,\cdot)$ (the previous PT-2 omitted A-CONT while asserting a function on every $z$; A-CONT is now included). For that version, $\ell(z,\cdot)$ is continuous and $\mu$-strongly convex on the compact convex $\Delta_m$ for every $z\in Z$, so
$$g^\star(z)=\operatorname*{arg\,min}_{p\in\Delta_m}\ \ell(z,p)$$
is a **single-valued function defined on all of $Z$** (existence by Weierstrass, uniqueness by strong convexity). The target of Phase 21 stands, now on a genuinely convex footing. $\square$

## 3. What is unchanged

The target is still the single risk-optimal map (Phase-21 Route B); only its coordinates changed ($\Delta_m$ in place of $[0,1]\times\Delta_{m-1}$), removing the bilinearity. Bayes-optimality over measurable maps (PT-3) transfers verbatim (same image $\mathbb B$, same risk). No operator redesign, no new generality — a coordinate change plus the correct hypothesis list.
