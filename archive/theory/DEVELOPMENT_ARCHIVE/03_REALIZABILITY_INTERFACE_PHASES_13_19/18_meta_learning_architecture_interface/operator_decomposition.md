# Meta-Operator Decomposition (§2)

> **Status:** Phase-18, 2026-08-03. The minimal factorization $A_\theta(S)=D_\theta(r(S))$ — mathematical objects only. New results **MI-4–MI-6**, tagged **[proved] / [declared]**.

---

## 1. The representation

**Definition MI-4 (canonical representation).**
$$r(S)\ =\ \big(\ b_{\mathrm{can}}(S),\ \ z(S)\ \big)\ \in\ \mathbb B\times Z,$$
- $b_{\mathrm{can}}(S)$: the canonical band vector — the frozen forced/compatible construction with its exact count-dependent margins, per declared index; a closed-form function of $S$ (order statistics of the observations enter here: e.g. per-location interval endpoints are max/min functionals — see MI-8);
- $z(S)$: the compact stratified statistic — context label $\kappa(S)$ and (at meta-training) fiber frequencies/counts; for the current task alone, $z$ reduces to $\kappa(S)$ plus the declared index arguments.
Both components are computable from observables; $\mathbb B\times Z$ is a fixed finite-dimensional compact set under the declared skeleton.

## 2. The decoder

**Definition MI-5 (decoder).**
$$D_\theta(b_{\mathrm{can}},z)\ =\ \Big(\ K\big((1-\lambda)\,b_{\mathrm{can}}+\lambda\textstyle\sum_{j=1}^m\varphi_j(z)\,b_j\big),\ \ \text{confidence}(z),\ \ \text{rung}(z)\ \Big),\qquad \theta=(\lambda,b_1,\dots,b_m)\in[0,1]\times\mathbb B^m,$$
with $K(\cdot)$ the fixed band-to-class denotation, $\varphi$ the declared partition of unity, and the side channels computed canonically from the counts in $z$. $D_\theta$ is affine in (reparameterized) $\theta$, Lipschitz in $r$ (Hoffman / $W_1$-stability constants), and valid for every $(\theta,r)$ (MR-9).

## 3. The two theorems that make this *the* decomposition

**Theorem MI-6 (sufficiency and family-minimality of $r$). [proved]**
(i) *Sufficiency:* every member of the family factors through $r$: $A_\theta=D_\theta\circ r$ for all $\theta$ — by construction (MR-3), since the decoder reads $S$ only through $(b_{\mathrm{can}},z)$.
(ii) *Minimality (relative to the family — the only well-posed sense):* the equivalence kernel of $r$ equals the joint equivalence kernel of the family: if $r(S)=r(S')$ then $A_\theta(S)=A_\theta(S')$ for **every** $\theta$ (immediate from (i)); conversely if $r(S)\ne r(S')$ then some family member separates them — a difference in $b_{\mathrm{can}}$ is exposed at $\lambda=0$ (the canonical member emits different bands), and a difference in $z$ is exposed either by a $\varphi$-weight difference at some $\lambda>0$ with distinct anchors, or by the side channels (confidence/rung read $z$). Hence no coarser representation is sufficient, and every finer one carries information the family provably cannot use. Minimal sufficiency is *family-relative* by necessity: a richer declared family would demand a richer $r$ — this scoping is honest and is re-examined in the failure audit. $\square$

**Remark (why the decomposition is forced, not chosen).** The split reproduces the program's central separation one more time: $r$'s first component is the **identification channel** (certificate-carrying, $\theta$-invariant), its second the **population channel** (preference-carrying, learnable); the decoder is the convex assembly that provably cannot mix them into a false claim (MR-13(iii)). Any decomposition that merged the channels inside a single opaque representation would violate the channel-typing rule MI-3.4 — so the minimal decomposition is also the only admissible one up to reparameterization. $\square$
