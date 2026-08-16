# Final Theory Stopping Criterion (Part VI)

> **Status:** Phase-6 terminal decision, 2026-08-03. Determines whether the mathematical program is sufficient for handoff, after the three final additions (auxiliary information, irregular archives, imperfect archives). Sources: Phase-6 files (CI-A/B/C/D) and the frozen corpus. The decision is audited against every open item the phase produced; the mandate forbids opening research directions unless one blocks the core result.

---

## The decision

$$\boxed{\textbf{THEORY\_COMPLETE\_FOR\_HANDOFF}}$$

The mathematical program is sufficient for handoff. Each of the three additions was **reduced to the frozen base theory**, introducing no new primitive and no blocking gap:

- **Auxiliary information** reduces to **conditioning on a fiber** (CI-A1): the base canonical operator on $\mathcal F_{c_b}$, with the exact useful-iff theorem (CI-A3) — *auxiliary information is useful iff it changes the joint window* — and the impossibility that its harm cannot be graded by a distance on $C$ without declared structure (CI-A5(iii)).
- **Irregular archives** reduce to a **sheaf-like consistency system with a gluing law** (CI-B): a counting necessity $\sum(k_a-d)_+\ge d(N-d)$ (CI-B1, real constant-rank proof), a per-point necessity (CI-B2), a decisive proof that connectivity is *not* the invariant (CI-B3), an exact gluing lemma across unisolvent overlaps (CI-B4a) with an explicit holonomy obstruction when overlaps degenerate (CI-B4b), and a chaining sufficiency (CI-B5) — necessary and sufficient conditions that **sandwich** the one open combinatorial item.
- **Imperfect archives** reduce to **Theorem 1 on a union family** (CI-C1): the hull-of-union operator is minimax-optimal with valid outer certificates, a conditional $O(\delta)$ perturbation bound under $\sigma_0$-conditioning (CI-C2), and a matching impossibility of any uniform bound without conditioning (CI-C3).

- **Combined** (CI-D): a single justified conditional operator $A_c(\text{archive},c_b,S_b,x,\varepsilon)$ = the base operator on the union-fiber family, optimal by composition of the three reductions, with an exact five-way information ledger and five decisive substitution tests (T1–T5).

**Why no missing theorem blocks the core.** The base theory (Phases 0–5) was already handoff-ready for declared, tame, stable classes with exact archives. Phase 6 relaxes the three simplifying assumptions, and in each case the relaxation is absorbed by an existing base theorem applied to a modified family — fiber, sheaf-glued window system, or union — so **correctness and optimality never required a result the corpus lacks.**

---

## Open items, each certified non-blocking

The mandate forbids opening research directions unless one blocks the core. None does; recorded for completeness.

| Open item | Why it does not block |
|---|---|
| **Exact combinatorial characterization of *unique* completability for general $d$** (CI-B8) | **Non-blocking.** CI-B1 (necessity) and CI-B5 (sufficiency) sandwich it; the operator uses whichever windows are *identified*, and honestly reports residual ambiguity (radius $=$ realized affine/gauge freedom) where they are not. The literature (Pimentel-Alarcón–Boston–Nowak) exactly characterizes *finite* completability and gives only sufficient conditions for *unique* — so this is open in the field, not merely in this corpus, and the operator is correct regardless of which side of the open line a pattern falls: it never over-claims. |
| **Sharp archive-noise perturbation constant** $C$ (CI-C2; the corpus's original OPEN item) | **Non-blocking.** The union reduction (CI-C1) makes the honest certificate *valid at any $\delta$* without the constant; the constant governs only the *tightness* (size $h$ of the union), not correctness. Demoted from a correctness gap to a quantitative-refinement question, exactly as CI-C1(iv) states. |
| **Exact optimality under member/archive coupling** (general non-linear classes, CI-D2) | **Non-blocking.** The union operator remains a **valid outer certificate** under coupling; only the exact-minimaxity claim weakens, and the exactly-$d$ linear class (the corpus's identification setting) is uncoupled by construction. |
| Inherited Phase-5 items (distributional guarantees; unstable classes) | Unchanged: distributional theory is deliberately out of scope (would need new axioms); unstable classes are settled *negatively* (CR-7), a prohibition, not a gap. |

None of these is a theorem whose absence prevents constructing, certifying, or falsifying the operator. Each is a refinement that would tighten a bound or widen a scope already correctly handled by conservative reporting.

---

## What handoff now delivers

The complete theoretical object handed off is: the conditional operator $A_c$, expressible as the base canonical operator (center $+$ compactified one-sided radius, with partiality flags) on the **union-fiber window system** — auxiliary information as fiber, irregular archive as sheaf-glued identification with certified holonomy, imperfect archive as union — under the Phase-5 realizability constraints (branch normal form, M1–M12) and the Phase-6 information ledger, falsifiable by P1–P10, NP-1…NP-6, and the new substitution tests T1–T5.

**Distinguishing the four epistemic statuses across the whole final phase:**
- **Proved:** CI-A1, CI-A2, CI-A3; CI-B1, CI-B3, CI-B4a, CI-B4b, CI-B6, CI-B7; CI-C1, CI-C4; CI-D1, CI-D2, T1–T5.
- **Conditional on stated assumptions:** CI-A (class-conditional augmentation transfer, linear label); CI-B2 (non-pivot genericity), CI-B5 (cover + local F17 + unisolvent overlaps); CI-C2 ($\sigma_0$-conditioning, $\delta\le\delta_0$); CI-D2 (well-specification, no-coupling, augmented-class $\mathcal W_{\mathrm{arch}}$).
- **Impossible:** CI-A5(iii) (no distance-graded harm without structure on $C$); CI-C3 (no uniform conditioning-free perturbation bound).
- **Open (non-blocking):** exact unique-completability characterization for general $d$; sharp constant $C$; exact optimality under coupling.

---

## Closing

The program that began with a pure impossibility theorem — *with no assumptions the data determine nothing off the design* — ends by absorbing its three most realistic complications into the same single primitive, the trace modulus, via three reductions to one minimax identity. Auxiliary information conditions it, irregular archives glue it, imperfect archives union it; the operator, its certificate, its ledger, and its falsification tests survive each relaxation intact. The theory is complete for handoff, and the boundary it draws between what remains provable and what is provably impossible is, itself, one of its results.

**Verdict: `THEORY_COMPLETE_FOR_HANDOFF`.**
