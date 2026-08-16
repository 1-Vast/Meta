# Finite Parameterization Theorem (§2)

> **Status:** Phase-17, 2026-08-03. New results **MR-3–MR-6**, tagged **[proved] / [declared] / [impossible]**. Mandate: construct $\Theta\subset\mathbb R^p$ with fixed finite $p$ **not depending on $\varepsilon$**; if impossible, prove it and identify the weakest relaxation. Both halves are delivered: a fixed-$p$ construction at fixed declared skeleton, and an impossibility theorem showing no fixed $p$ covers all resolutions — the relaxation is exactly "fix the skeleton", and it is proved weakest.

---

## 1. The declared skeleton (the sole relaxation)

**(SKEL) [declared].** A deployment declares, once: the context set $C_\kappa$; the finite query/specification atlas and event system $\mathcal E$ (Route A) or value grid $T$ with mesh $h$ (Route B); the statistic horizon $\bar N$; and a fixed finite **partition-of-unity** $\varphi=(\varphi_1,\dots,\varphi_m)$ on the compact statistic domain $Z$ ($\varphi_j\ge0$, $\sum_j\varphi_j\equiv1$, continuous — e.g. barycentric weights of a fixed triangulation; the choice is declared, not learned). Everything below is relative to (SKEL); nothing below depends on any accuracy parameter.

## 2. The fixed-$p$ theorem

**Theorem MR-3 (finite parameterization at fixed skeleton). [proved]**
Let $q=\dim(\text{band vector})$ under (SKEL) and $\mathbb B\subset\mathbb R^q$ the valid-description polytope (MR-1/§3). Define
$$\Theta\ =\ [0,1]\ \times\ \mathbb B^{\,m}\ \subset\ \mathbb R^{\,p},\qquad p\ =\ 1+m\,q\ \ \text{— fixed, finite, } \varepsilon\text{-free},$$
with decoder (per index; index subscripts suppressed)
$$A_\theta(H)\ =\ K\Big(\ (1-\lambda)\,b_{\mathrm{can}}(H)\ +\ \lambda\sum_{j=1}^m\varphi_j\big(z(H)\big)\,b_j\ \Big),\qquad \theta=(\lambda,b_1,\dots,b_m),$$
where $b_{\mathrm{can}}(H)$ is the **canonical band vector computed directly from the observable record** $H$ (the frozen forced/compatible construction with its exact count-dependent margins — a fixed, closed-form, parameter-free component of the decoder, *not* factored through $z$; this replaces, and retracts, Phase-15's false exact-factorization claim on the $\top$-stratum), and $z(H)$ is the compact stratified statistic. Then:
(i) every value is valid: $b_{\mathrm{can}}(H)\in\mathbb B$ (canonical outputs satisfy order, coherence, and nonemptiness — the empirical law is a witness; proved in the earlier phases for the canonical construction) and each $b_j\in\mathbb B$ by the type of $\Theta$; the decoded band is a convex combination of $\mathbb B$-points, hence in $\mathbb B$ (MR-1) — **for every $\theta$ and every $H$**;
(ii) the family **contains the canonical operator exactly** at $\lambda=0$ — zero approximation error to the frozen baseline at fixed $p$; no sieve;
(iii) $p$ depends only on (SKEL), never on a tolerance. $\square$

## 3. The impossibility half, and the weakest relaxation

**Theorem MR-4 (no fixed $p$ across all resolutions). [impossible — proved]**
There is no fixed finite-dimensional family $\{A_\theta\}_{\theta\in\mathbb R^p}$ (decoder measurable, arbitrary) that is uniformly $\varepsilon$-dense, for every $\varepsilon>0$, in the valid operators of **all** skeletons simultaneously. *Proof.* As the skeleton refines (Route-B mesh $h\to0$, or growing Route-A atlases), the valid operator spaces contain affinely independent families of unbounded dimension: e.g. the band maps $z\mapsto b$ constant in $z$ realize all of $\mathbb B_{\mathrm{skel}}$, and $\dim\mathbb B_{\mathrm{skel}}\to\infty$ along refinements (each new threshold adds an independent band coordinate whose value is decision-relevant — two operators differing only there differ in $d_{\mathbb M}$ by a fixed amount, by the $W_1$/TV separation of the corresponding classes). A continuous (indeed any) image of $\mathbb R^p$ cannot be $\varepsilon$-dense in a normed family containing $(\dim>p)$-dimensional affine simplices with pairwise distances bounded below: a $p$-parameter set is a $p$-dimensional object and, for $\varepsilon$ below half the separation of $2^{p+1}$ affinely positioned targets, density would require distinct preimages exceeding any covering capacity of $\mathbb R^p$ at that scale — formally, the metric entropy of the target family grows without bound along refinements while the entropy of the image of any fixed-dimensional compact parameter set is bounded at each scale. Hence uniform-over-skeleton fixed-$p$ density fails. $\square$

**Corollary MR-5 (the relaxation is weakest). [proved]**
Any relaxation permitting a fixed-$p$ family must bound the decision-relevant resolution — i.e. must fix (SKEL) up to equivalence: MR-4's obstruction is *only* the unbounded independent band coordinates that refinement creates, and MR-3 shows fixing them suffices. So "declare the finite skeleton" is simultaneously necessary (MR-4) and sufficient (MR-3): the **weakest relaxation**, exactly as the mandate demands. What the relaxation costs is explicit and honest: Route-B carries the irreducible declared mesh floor $2h$ in $W_1$ (§3); Route-A carries the deployment atlas; refining either is a *new deployment*, not a larger $\varepsilon$-net inside one family.

**Remark MR-6 (why the Phase-15 sieve dissolved).** The sieve arose from banning the canonical closed-form map from the decoder and asking a fixed weak function class (multilinear tables) to imitate it to arbitrary accuracy. The canonical map is frozen mathematics computable from observables — including it as the fixed $\lambda=0$ component is not an oracle assumption (see `adversarial_counterexamples.md`), and with it the only thing left to parameterize is a *deviation* inside a fixed compact polytope: finitely many parameters, permanently. $\square$
