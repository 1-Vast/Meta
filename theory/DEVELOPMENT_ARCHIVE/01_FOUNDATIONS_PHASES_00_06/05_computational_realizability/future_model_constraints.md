# Future Model Constraints

> **Status:** Phase-5 interface, 2026-08-02. This file answers exactly one question: **"What mathematical properties must a future computational model satisfy?"** It does not answer how the model should be implemented. Every item cites its theorem. Restructured per the refereed plan audit: the monotonicity item is split (radii required monotone, center exempt); the branch/flag item is the certified *disjunction*, not a prescription; declared-class items are added; the MAY clause is scoped so it cannot be read against the MUSTs.

---

## MUST (each violation is certified error)

| # | Property | Certificate |
|---|---|---|
| M1 | **Declare its class, tameness route, and stability certificate**: the closure class it assumes (no closure assumption ⇒ no valid radius), whether tameness is by definability or closed form, and where the modulus is locally bounded | handoff §6; DM-6/CR-1; CR-7 |
| M2 | Approximate the **four objects** O1–O4 — class representation, task adaptation, query readout, coverage/validity decision — with O4 separate and one-sided (false inclusion = fabrication; false exclusion = conservative) | Phase-4 interface; CR-8 |
| M3 | Output the **pair** (center, radius), radius in $[0,+\infty]$ with $+\infty$ a lawful value; envelopes **outer** (one-sided); loss of any kind priced outward by inflation, never absorbed | IB-7; CR-6.2 |
| M4 | Take $\varepsilon$ as an input; radii nondecreasing in $\varepsilon$ and nonincreasing under support refinement — **the center is exempt from all monotonicity** | A4/A5; OP-8 (split per audit) |
| M5 | Satisfy the composite invariances: permutation of support pairs; joint affine equivariance including reflections; gauge invariance (nothing reported may depend on latent coordinates or basis choices) | A1/A2/A3; CR-5 |
| M6 | Respect the capacity ceilings: at most $k$ continuous task dimensions; the size-$(k{+}1)$ window truncation; latent dimension $\ge$ class dimension (continuous **and** definable categories, with topological excess under the full-class hypothesis; equivariant representations budgeted for their extra width) | F20/CP-3; DM-3; DM-2 |
| M7 | Handle certified discontinuities by the **disjunction**: piecewise structure aligned with the (lower-dimensional, definable) transition sets, **or** declared localized transition-band error / certificate inflation $\sim J/2$ — soft or randomized selection does not create a third option | IB-9; MP-4; CR-3 |
| M8 | Guard composition by the **margin lemma**: selector decisions trusted only at margin $\ge2h$; the collar flagged one-sidedly | CR-4 |
| M9 | Surface the three partiality flags — unrealizable support (misspecification detector; silent projection incurs the doubling cost), off-coverage, unbounded section — which exhaust the undefined locus | IB-8; CR-8 |
| M10 | Carry **data-dependent certificate content** outside the surjective-trace stratum; treat DM-5 stratum membership as a flagged piecewise object; not defer $\varepsilon$ to the readout outside strata (i)/(iii) | DM-5; AP-5 |
| M11 | Reproduce exact identifiable data with zero error; keep the reported center $1$-Lipschitz in the reported envelopes; satisfy the sensitivity-sum on constants-including classes | A7; IB-10; IB-11 |
| M12 | Be evaluated only by the worst-case protocols (P1–P10, NP-1–NP-6) and compared to other systems only in the induced-operator metric on the double gauge quotient | F3/IB-14; DM-9/NP-6 |

## MUST NOT

- Emit finite radii at unstable-local configurations while carrying any representation error (CR-7(c): pointwise void — unconditional), or anywhere off its declared coverage without declaring the member-level assumption (F18).
- Attach semantics to latent coordinates, compare latents across systems, or claim a globally continuous basis-style canonical form on topologically nontrivial classes (CR-5).
- Use symmetric error notions for envelopes (false certificates), suppress the misspecification flag, operate $\varepsilon$-free, or impose monotone center updates.
- Claim more than $k$ task dimensions, structure beyond the $(k{+}1)$-truncation, generic identification of nonlinear classes with $d\ge3$ at $k\le5$ without per-configuration verification, or any average-case vindication of a worst-case obligation.

## MAY (unconstrained beyond the listed extensional constraints)

Choice of gauge locally (chart-wise bases are fine — the obstruction is global); choice of branch representation and of any internal structure whatsoever; smooth surrogates anywhere the class is stable and the margin is guarded — provided the extensional constraints above hold. The mathematics is exact about *what* the model computes and silent about *how*; this clause may not be read against M6's floors or M7's disjunction.

---

*One sentence:* a valid future model is any computational object that — for a declared, tame, stability-certified class — computes outer envelope pairs and flags satisfying M1–M12; everything else about it is mathematically free.
