# Parameter Space (Deliverable 1)

> **Status:** Phase-20 (trainable operator foundation), 2026-08-03. Phases 0–19.5 unmodified. Audit of record: `../20_final_theory_audit/BLOCKING_ISSUES.md` — the minimal obstruction is that $F_\omega:Z\to C$ was introduced with $\Omega$, $\mathcal H$, and the empirical objective undefined, and C3 imposed rather than derived. This phase defines all of them, architecture-free. New results carry **TF-** numbers, tagged **[proved] / [declared] / [impossible]**. No architecture, no implementation vocabulary, no application.

---

## 1. The fixed context (carried, not re-derived)

Under a declared deployment skeleton (SKEL, Phase 17): the statistic domain $Z$ is a finite union of compact cubes (compact metric); the **coefficient set**
$$C\ =\ [0,1]\ \times\ \Delta_{m-1}\qquad(\text{mixing weight }\lambda\ \times\ \text{partition weights }w\text{ over the }m\text{ declared anchors})$$
is compact convex; the fixed **assembly** $\mathsf{asm}:C\times\mathbb B^m\times\mathbb B\to\mathbb B$, $\mathsf{asm}\big((\lambda,w),(b_j),b_{\mathrm{pop}}\big)=(1-\lambda)b_{\mathrm{pop}}+\lambda\sum_jw_jb_j$, lands in the valid-description polytope $\mathbb B$ for every argument (convexity, MR-9) and is Lipschitz (band coordinates affine). The learned object is exactly the map $Z\to C$ plus the anchor choice — nothing else.

## 2. The parameter space $\Omega$

**Definition TF-1 (implementation parameter space). [declared, generic]**
$$\Omega\ =\ \Xi\ \times\ \mathbb B^{m},\qquad \Xi\subset\mathbb R^{D}\ \text{nonempty compact},\ D\in\mathbb N\ \text{declared},$$
with the Euclidean metric and Borel σ-algebra. $\omega=(\xi,b_1,\dots,b_m)$: $\xi$ parameterizes the coefficient map; $b_1,\dots,b_m\in\mathbb B$ are the anchor bands (finite by SKEL). $\Omega$ is a compact subset of $\mathbb R^{D+m\dim\mathbb B}$ — **finite-dimensional, compact, explicit**; the only architecture-dependent quantity is $D$, and it enters solely as a dimension.

## 3. The realization map (the architecture slot, typed)

**Definition TF-2 (realization map and its regularity — the sole interface to "an architecture"). [declared]**
An **admissible realization** is a map
$$G:\ \Xi\times Z\ \to\ \mathbb R^{\dim C},\qquad(\xi,z)\mapsto G(\xi,z),$$
required to satisfy **only** (checkable, architecture-neutral) conditions:
- **(G1) joint continuity** on the compact $\Xi\times Z$ (hence uniform continuity and joint measurability);
- **(G2) parameter-Lipschitz:** $\|G(\xi,z)-G(\xi',z)\|\le L_\xi\|\xi-\xi'\|$ uniformly in $z$, declared $L_\xi<\infty$;
- **(G3) input-regularity:** a declared modulus of continuity $\varpi_G$ with $\|G(\xi,z)-G(\xi,z')\|\le\varpi_G(d_Z(z,z'))$ uniformly in $\xi$.
Any construction meeting (G1)–(G3) is admissible — the definitions and theorems below quantify over *all* such $G$; a specific dense witness family is exhibited in `approximation_theorem.md` (TF-9) to prove the class is non-empty and the approximation achievable. **No architecture is named; (G1)–(G3) are the entire contract an architecture must meet.**

**Proposition TF-3 (projection to the coefficient set is free). [proved]**
Let $\pi_C:\mathbb R^{\dim C}\to C$ be the Euclidean metric projection onto the compact convex $C$. Then $\pi_C$ is well-defined and $1$-Lipschitz (projection onto a convex set is nonexpansive), hence continuous and measurable. Define the **coefficient map**
$$F_\xi\ :=\ \pi_C\circ G(\xi,\cdot)\ :\ Z\to C.$$
$F_\xi$ inherits (G1)–(G3) with the same constants (composition with a $1$-Lipschitz map): jointly continuous in $(\xi,z)$, $L_\xi$-Lipschitz in $\xi$, modulus $\varpi_G$ in $z$, and **$C$-valued for every $\xi$ and every raw $G$** — so no feasibility obligation is placed on the architecture; validity of the coefficient is a property of the fixed projection, not of the learned weights. $\square$

This is the first half of the audit's repair: $\Omega$ exists, is compact and finite-dimensional, and the map $F_\omega:Z\to C$ is now a defined object (TF-3) for any admissible realization, with its regularity constants declared and its codomain guaranteed by construction. The hypothesis class and objective follow in the next files.
