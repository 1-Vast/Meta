# Operator Approximation Theory

> **Status:** Phase-5 derivation, 2026-08-02. Sources: frozen corpus + Phase 4 (DM-, IB-, O-numbers). New results carry CR-numbers; all were adversarially refereed and the corrections are incorporated (notably: the two-routes scoping of the normal form, the compactified-value completion of the regularity corollary, and the formalized soft-selection dichotomy). No parameterization is chosen anywhere.

**Question.** Given $\mathbb A(\mathbb T,S,x,\varepsilon)$ and the Phase-4 factorization $\Phi/U/R$: what is the most general differentiable approximation form, and what must be continuous / may be discontinuous / needs piecewise representation / needs explicit validity flags?

---

## 1. Two normal forms, not one (CR-1)

There is no single "most general differentiable form"; the corpus's two tameness routes (DM-6) yield **two normal forms**, refereed as genuinely distinct:

**(A) Definable route.** For a uniformly definable class, for every finite $r$:
$$\mathbb A\;=\;\sum_{b=1}^{B}\mathbf 1_{C_b}\cdot A_b\;\ \oplus\ \text{flags},$$
with finitely many definable cells $C_b$, each branch $A_b$ of class $C^r$ on its cell (full-dimensional cells are their own neighborhoods; lower-dimensional cells extend via definable tubular neighborhoods — extension to cell *closures* can fail by derivative blow-up, a hedge doing real work), the selector $\mathbf 1_{C_b}$ definable and discrete, and the $\pm\infty$ cells routed to the flag/constant branch, never to the smooth machinery. Finiteness of $B$ is an *achievability* guaranteed by o-minimal $C^r$ cell decomposition, not a necessity. The decomposition may be taken **jointly in $(\theta;S,x,\varepsilon)$** — the class parameter included — so meta-level (archive-side) transitions (NP-3) are cells of the same decomposition, and the DM-5 $\varepsilon$-placement strata (rank conditions on $G$ are definable) are distinguished cells of it: stratum membership is itself a flagged, piecewise object, not a side condition.

**(B) Closed-form route.** For Lipschitz classes and RKHS balls (not definable families), the envelopes are explicit finite formulas (McShane envelopes; Golomb–Weinberger interval): continuous **on their domain of definition** for each fixed $k$, with all discontinuity carried by the *domain boundary* (the moving realizability threshold) — i.e., by the flag, not the value. *Subsumption remark (hedged):* for specific kernels definable in a richer o-minimal structure (e.g. Gaussian, definable in $\mathbb R_{\exp}$) the closed-form route embeds into route (A) at each fixed $k$; the routes differ in uniformity in $k$, which is neutralized within scope by $k\le5$: all cell/branch counts may be tabulated over $k\in\{0,\dots,5\}$.

---

## 2. Regularity across branches: jumps only, for tame classes (CR-2)

**Corollary CR-2 (jump-type regularity; refereed — a corollary of DM-6(a) plus the o-minimal Monotonicity Theorem, not a new theorem).** Route (A): along any *definable* curve $\gamma(t)\to p$ in the operator's domain, the envelope composition is definable in one variable; the $+\infty$ locus along the curve is definable, so near $0^+$ the composition is either identically $+\infty$ (limit trivially exists) or real-valued, where the Monotonicity Theorem's endpoint addendum gives one-sided limits in $\mathbb R\cup\{\pm\infty\}$. **Every discontinuity is jump-type along definable curves** — no oscillation. Scope: per-definable-curve (non-definable spirals may still oscillate); route (B) has continuous envelopes on-domain, so its only "jumps" are flag flips.

**Necessity of tameness (refereed witnesses).** A merely continuous family can oscillate: $\Theta=\overline{\{(t,\sin(1/t)):t\in(0,1]\}}$ (topologist's sine curve, compact), $f_{(t,s)}(x_1)=t$, $f_{(t,s)}(x)=s$ — the upper envelope at data $y$ is $\sin(1/y)$ for $y>0$ and $1$ at the realizable datum $y=0$: $\limsup=1$, $\liminf=-1$. Tameness is exactly what buys jump-only discontinuity structure.

---

## 3. The continuity taxonomy (mandate items, answered)

| | Object | Certificate |
|---|---|---|
| **Must be continuous** | each branch $A_b$ on its cell (indeed $C^r$, any finite $r$); the center *given* the envelopes ($1$-Lipschitz, constant sharp); route-(B) envelopes on their domain | DM-6(b); OP-5/IB-10 |
| **May (must be allowed to) be discontinuous** | the branch selector; the flags; the radius at validity boundaries (jump to $+\infty$; codomain compactified to $[0,+\infty]$); the operator across genuine section-topology transitions | MP-4; DM-6(d); O4 |
| **Needs piecewise representation** | the envelope pair across transitions (route A); the DM-5 strata; the meta-level (archive-side) transition structure | DM-6; NP-3; CR-1 |
| **Needs explicit validity flags** | the three partiality sources — unrealizable support (misspecification detector; silent projection costs the doubling penalty), off-coverage, unbounded sections — plus the coverage/validity decision object with its one-sided semantics | IB-8; O4; F18 |

**Lemma CR-8 (exhaustiveness of the flags).** For a declared tame class, the locus where the branch form emits no finite value is exactly partitioned by the three flags: empty section (unrealizable), unknown window (off-coverage — *epistemic*, an assertion about the declared class's knowledge, with one-sided semantics), and nonempty unbounded section (structural). Every undefined point falls under exactly one, each with its own one-sided error semantics; there is no fourth failure mode of the value map.

---

## 4. The selection cannot be smoothed away (CR-3)

**Proposition CR-3 (soft-selection no-escape; refereed with the formalization it needs).** Formalize "end-to-end differentiable" as: *all non-flag outputs — including any selection distribution — depend continuously on the input.* Then at a genuine jump of size $J$:
- a **soft selector** composed with continuous branches is a continuous map, so MP-4 applies directly: sup-error $\ge J/2$ near the boundary;
- a **continuously randomized selector** has mixing weight passing through $\tfrac12$, so worst-case expected error $\ge\sim J/2$;
- a **discontinuously randomized selector** exposes a discrete object in distribution — i.e., it *is* the flag horn.

Hence the dichotomy is exhaustive: **expose the selector as a discrete output, or accept a localized half-jump error band** — with the honest third reading of the second horn being *certificate inflation*: a system may widen its reported radius by $\sim J/2$ on the band and remain valid (the error is priced into the certificate rather than hidden).

---

## 5. Composition stability: the margin lemma (CR-4)

Branch-wise accuracy does not compose for free: upstream error $h$ (in $\Phi$ or $U$) perturbs envelope values and can **flip the selector near envelope ties**, causing $O(\text{local jump})$ composite error at points *off* the transition set.

**Lemma CR-4 (selector margin with flagged collar).** Define the margin $\mu(\text{input})$ as the relevant envelope gap governing branch selection. If all factor-wise sup-errors are $\le h$, then on $\{\mu\ge2h\}$ the selector is unflipped and the composite error is bounded by the branch error; on the **definable collar** $\{\mu<2h\}$ the system must emit the one-sided uncertainty flag (or inflate the radius by the local jump). On compact subsets of open cells the collar shrinks to the transition set as $h\to0$, so factorized approximants compose convergently. This is the bridge between per-factor approximation and whole-domain certified behavior; without it, the Phase-4 factorization's error analysis is incomplete.

---

## 6. Summary

The most general differentiable approximation form compatible with the corpus is: **finitely many $C^r$ branches on definable cells (jointly in class, support, query, and $\varepsilon$), a discrete selector, one-sided flags exhausting the undefined locus, margin-guarded composition, and certificate inflation as the lawful alternative to discreteness** — in route (A); with route (B)'s explicit continuous-on-domain envelopes as the parallel form for the nonparametric regimes. Anything smoother is provably wrong at the jumps (CR-3); anything less structured discards regularity the mathematics guarantees for free (CR-2).
