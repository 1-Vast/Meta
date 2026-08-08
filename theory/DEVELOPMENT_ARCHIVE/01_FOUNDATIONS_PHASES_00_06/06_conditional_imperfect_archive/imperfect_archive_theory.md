# Imperfect-Archive Theory (Part III)

> **Status:** Phase-6, 2026-08-03. Replaces exact historical values by bounded ones $|y_{ai}-f_a(x_{ai})|\le\delta_a$. Sources: frozen corpus (Theorem 1, CP-2; the archive-noise quantitative constant was the corpus's flagged **OPEN** item, handoff §9.1). New results **CI-C**; tags **[proved] / [conditional] / [impossible] / [open]**. Adversarially refereed (**confirmed**, with hypotheses added); the corrections — the $\delta$-smallness threshold in CI-C2, the no-coupling/well-specification hypotheses in CI-C1, the surrogate-modulus specification, and the operational-metric reading of CI-C3 — are incorporated.

**The central question.** How does archive uncertainty propagate: archive error $\to$ window-system error $h$ $\to$ query-ambiguity inflation? Is there a uniform bound, and if not, what is the strongest impossibility?

---

## 1. The union-family reduction — archive noise needs no new optimality theory

**Definition.** For a declared class and archive data with bounds $(\delta_a)$, let the **consistent set**
$$\mathcal W_{\mathrm{arch}}=\{\text{candidate window systems }W\text{ in the declared class}:\ \text{the archive is }(\delta_a)\text{-consistent with }W\},\qquad \mathcal F_{\mathrm{union}}=\bigcup_{W\in\mathcal W_{\mathrm{arch}}}\mathcal F_W.$$

**Theorem CI-C1 (union reduction). [proved under well-specification + no-coupling]**
(i) Under **well-specification** (the truth's window system lies in the declared class and the noise bounds hold), $\mathcal F_{\mathrm{union}}$ is nonempty and contains the truth (the truth is $(\delta_a)$-consistent with itself).
(ii) Applying **Theorem 1 verbatim** to $\mathcal F_{\mathrm{union}}$: the canonical operator on the union — report the hull of the union of sections, center + radius — is exactly minimax-optimal, **conditional on the observed archive**, against the joint adversary (choice of consistent window system *and* member *and* noise). The equivalence "joint adversary $=$ member-adversary over $\mathcal F_{\mathrm{union}}$" is an exact projection identity **under the no-coupling hypothesis**: for each $W$, the admissible new member ranges over $\mathcal F_W$ independently of the archive members' realizations (automatic in the exactly-$d$ linear class, where $\mathcal F_W=W$ and members are uncoupled elements). Sections commute with unions, so the operator is literally hull-of-union-of-sections; risk depends on $(W,f)$ only through $f$, and the truth-set depends on the archive only through $\mathcal W_{\mathrm{arch}}$ (CP-2-style sufficiency — no extra exploitable information).
(iii) Its certificates are **valid outer certificates** by construction (the union covers every consistent candidate). **[proved]**
(iv) **Hence archive noise requires no new optimality theory, and the OPEN sharp perturbation constant is bypassed for validity purposes** — it remains open only as a question about the *size* $h$ of the union, never about correctness. The empty-union case (well-specification violated) triggers the misspecification-detection clause, not an estimate. **[proved]**

This is the central move of Part III: **imperfect archives are a special case of Theorem 1 applied to a larger family.**

---

## 2. The conditional quantitative perturbation theorem

**Theorem CI-C2 (perturbation, with the required smallness threshold). [conditional]**
Exactly-$d$ linear class; archive pattern satisfying the exact-case sufficient conditions (F17-type locally + unisolvent-overlap chaining, CI-B5); all $|y_{ai}-f_a(x_{ai})|\le\delta$; values bounded by $M$; **conditioning hypothesis**: every linear system solved in the identification (local cores, per-point systems, gluing steps) has smallest singular value $\ge\sigma_0>0$. Then, **provided $\delta\le\delta_0(\sigma_0,M,d,\text{transport})$** (small enough that every perturbed system keeps $\sigma_{\min}\ge\sigma_0/2$ by Weyl, and the F17 rank certificates persist), the consistent set has operational diameter
$$h\;\le\;C\cdot\delta,\qquad C=C(\sigma_0,M,d,\text{transport structure}),$$
via the classical per-solve bound $\|\Delta x\|\le(\|e\|+\|E\|\,\|x\|)/(\sigma_0-\|E\|)$ composed along the finite transport, and the final query inflation is
$$\le\ \tfrac12\,\omega^{\mathrm{surr}}_{x,D}(2\varepsilon+2h)+h$$
where $\omega^{\mathrm{surr}}$ is the modulus of the **surrogate (union/identified) family** (refereed specification; using the true modulus adds $O(h)$ terms).
*Refereed provisos (all load-bearing, none decorative):* (a) the smallness threshold is required — **for $\delta\ge M$ the zero member $g_a=0$ lies in every subspace and matches all data, so $\mathcal W_{\mathrm{arch}}=\operatorname{Gr}(d,N)$ and $h$ is unbounded** even though the true instance is well-conditioned; (b) $M$ is what closes the large-$\delta$ regime up to the threshold; (c) $C$ is **exponential in transport depth** ($C\sim C_{\mathrm{step}}^{L}$ over $L$ gluing steps) — legitimately inside "depends on the transport structure," and forced by the CI-C3 lower bound; (d) the quantitative statement is restricted to the covered region (at uncovered $x$ the union modulus is already large at $\delta=0$, and Theorem 1 still returns the honest radius, so no error).

The schematic chain **archive error $\to$ window error $h\le C\delta$ $\to$ query inflation $\le\tfrac12\omega(2\varepsilon+2h)+h$** is therefore sound with explicit, if regime-restricted, constants.

---

## 3. Impossibility of a uniform (conditioning-free) bound

**Theorem CI-C3 (no uniform bound). [impossible]**
No bound $h\le C(\delta)$ independent of conditioning exists. Two refereed witnesses:
- **Near-degenerate core.** $d=2$, core value-vectors at angle $\sim\sigma_{\min}$ (explicit: core values $(1,0),(1,\gamma)$, third-point values $(0,\gamma)$, all F17 rank hypotheses holding; $\sigma_{\min}\sim\gamma$). A $\delta$-shift $\gamma+\delta$ is matched by a subspace whose off-core predictions swing by $\Theta(\delta/\gamma)=\Theta(\delta/\sigma_{\min})$ — **unbounded at fixed $\delta$ as $\sigma_{\min}\to0$**, in the **operational/value-diameter metric** (the Grassmannian *angle* saturates on the compact space — this is the correct reading; the unbounded quantity is the value/row diameter, which is what the operator reports). Taking $\sigma_{\min}=\delta^2$ kills any modulus $C(\delta)\to0$.
- **Chain amplification.** Transport along $L$ near-degenerate unisolvent overlaps realizes multiplicative growth $\kappa^L$ (already at $d=1$: gluing at near-zero overlap values gives per-step relative factor $\kappa=(v+\delta)/(v-\delta)>1$). The growth is **realized, not an upper-bound artifact**, so no polynomial-in-$L$ uniform bound exists even with per-step condition number fixed at $\kappa>1$.

Both are correctly framed as impossibility of **uniform** bounds only: each fixed well-posed instance retains a finite bound via its own $\sigma_{\min}$ (consistent with CI-C2), so $h(\delta)\to0$ as $\delta\to0$ remains possible per instance.

**The additional regularity required (to escape CI-C3):** exactly the $\sigma_0$-conditioning hypothesis of CI-C2 — a floor on the smallest singular value of every solve in the transport, plus the $\delta\le\delta_0$ threshold. This is the archive-side analogue of the corpus's stability requirement (a positive modulus floor); without it the propagation is genuinely unbounded, which is the strongest impossibility available and answers the mandate's "if no uniform bound exists, provide the impossibility and identify the missing regularity."

---

## 4. $\delta$-monotonicity

**Proposition CI-C4. [proved]** $\mathcal W_{\mathrm{arch}}$ is nondecreasing (setwise) in each $\delta_a$ (constraint relaxation), hence union sections are nested and certified radii are nondecreasing in $\delta$ — the archive-noise analogue of the established $\varepsilon$-monotonicity. Centers may move (not claimed monotone), mirroring the estimate/guarantee split of the base theory (OP-8).

---

## 5. Summary of Part III

| Result | Statement | Tag |
|---|---|---|
| CI-C1 | union-family reduction: hull-of-union operator is minimax-optimal + valid outer certificate; OPEN constant bypassed for validity | **proved** (well-spec + no-coupling) |
| CI-C2 | $h\le C\delta$ under $\sigma_0$-conditioning and $\delta\le\delta_0$; $C\sim C_{\mathrm{step}}^L$; inflation $\le\tfrac12\omega^{\mathrm{surr}}(2\varepsilon+2h)+h$ | **conditional** |
| CI-C3 | no uniform conditioning-free bound (value-diameter swing $\Theta(\delta/\sigma_{\min})$; $\kappa^L$ chains) | **impossible** |
| CI-C4 | radii nondecreasing in $\delta$ | **proved** |
| — | sharp perturbation constant $C$ (exact, not schematic) | **open (non-blocking)** |

**Net answer to the mandate.** The chain *archive error $\to$ window error $h$ $\to$ query inflation* **holds with explicit constants under $\sigma_0$-conditioning + $\delta\le\delta_0$** (CI-C2), and **provably cannot hold uniformly without conditioning** (CI-C3). Crucially, **validity never needed the constant**: the union reduction (CI-C1) makes the honest certificate correct at *any* $\delta$, demoting the corpus's OPEN archive-noise item from a correctness gap to a mere *tightness* question about the size of the union.
