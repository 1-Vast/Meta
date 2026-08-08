# Decision-Extension Stopping Criterion (Part XV)

> **Status:** Phase-7 terminal decision, 2026-08-03. Determines whether the decision extension is complete. Sources: the nine Phase-7 files (DE-S, DE-P, DE-H, DE-J, DE-R, DE-L, DE-T, DE-U, DE-O numbers) and the frozen corpus (cited only). The mandate's constraint is honored throughout: the identification theory was not modified, reinterpreted, or extended; every Phase-7 object consumes the frozen $I(O)$ (or its Phase-6 conditional/union form) as an opaque input.

---

## The decision

$$\boxed{\textbf{DECISION\_EXTENSION\_COMPLETE}}$$

The phase's single question — *given an honestly identified set, what additional mathematical object permits a learnable and scientifically honest decision under residual ambiguity?* — is answered, with proofs, and the answer is closed under its own audit:

**The object is a monotone total-preorder completion of the dominance order on loss profiles (DE-P2/P3) — concretely generated, when history is to supply it, by an ambiguity class of population laws conditioned inside the identified set (DE-H5, DE-L1), legitimate exactly under declared (EXCH/IID) + (LIK) + (COV) (+ transport class under shift), learnable at explicit $n^{-1/2}$ rates where historical members' own values are identified and up to an $n$-irreducible censoring width where they are not, and collapsing conservatively to the frozen minimax operator when nothing is declared.**

---

## Audit against the five mandated theorems (XII)

| Required | Delivered | Status |
|---|---|---|
| 1. Impossibility: the identified set alone does not determine a non-minimax decision | DE-S3 (two-context underdetermination), DE-S4 (equivariance: intrinsic rules are pinned to the Chebyshev center on symmetric sets), DE-S5 (no canonical weighting; gauge obstruction) | **proved** |
| 2. Sufficiency: under explicit added decision information an optimal action exists | DE-O2 (existence + undominated selection; via DE-T4(i) l.s.c./compactness), DE-P3(i) (Szpilrajn existence at the preorder level) | **proved** |
| 3. Learning theorem/bound for the decision object from historical members | DE-L3 (Hoeffding/DKW + BV risk transfer, explicit constants), DE-L4 (censored-history interval bound; systematic width theorem) | **proved, conditional on declared (IID), (COV)** |
| 4. Robustness / impossibility under population shift | DE-T2 (Λ-ratio and TV-ball degradation formulas; robust ranking threshold), DE-T3 (undeclared shift: history worthless, committal rules harmful), DE-T4 (Γ-minimax well-posedness; frozen endpoint) | **proved** |
| 5. Joint-query theorem for ranking decisions | DE-J6/DE-J8 (marginal insufficiency; joint pushforward necessary), DE-R1–R6 (sign trichotomy; comparative probability $p_{ab}$ as the exact minimal object; three-tier hierarchy; second-order partial identification) | **proved** |

Further mandated determinations, all discharged: formal separation $I$ vs $D$ (Part I file); minimal primitive without presupposing any candidate (Part II — the candidates are *located* as completions, none privileged without proof); feasibility/frequency separation with the no-shrinkage prohibition as a theorem (DE-H1–H3); single-law-vs-set with logical-status comparison and the two independent sources of set-valuedness (DE-H4/H5); criterion comparison with divergence witnesses and the prediction collapse of regret (DE-H6); five joint-loss regimes with the boundary theorem (DE-J2–J8); error decomposition with the which-term-shrinks corollary (DE-U7/U8); the four-quantity ledger with pairwise non-determination (DE-U1–U6); the honest operator with H1–H6 and the reduction, inertness, and consistency theorems (DE-O1–O4). No computational design appears anywhere in the phase (mandate XIII).

**Epistemic classification across the phase.**
- **Proved:** DE-S1–S5, DE-P1–P6, DE-H1–H6, DE-J1–J8, DE-R1–R6, DE-L1, DE-L2, DE-L5, DE-T1, DE-T3, DE-T4, DE-U1–U8 (decomposition conditional as tagged), DE-O1–O4.
- **Conditional on declared axioms:** DE-L3/L4 (IID, COV), DE-T2 (declared transport class), DE-U7 (the full declaration stack), all Tier-2/3 ranking outputs.
- **Impossible (proved impossibilities):** canonical non-minimax selection from the set alone (DE-S4/S5); a single posterior from the frozen noise model (DE-H4); frequency-driven shrinkage of $I(O)$ (DE-H2/H3); decision learning under undeclared shift (DE-T3); reconstruction of any ledger row from the others (DE-U1–U6); unconditional sub-radius guarantees (DE-L5/DE-O4).

---

## Open items, each certified non-blocking

| Open item | Why it does not block |
|---|---|
| Exact axiomatics of **menu-dependent** coherent selectors (minimax regret's tier, DE-P4): a characterization theorem for the (preorder, reference) pairs | The phase's primitive covers all menu-independent coherent selection; regret is correctly *priced* (one extra declared object) and usable as-is; only its abstract classification is open — a taxonomy refinement, not a correctness gap. |
| Sharp constants and finite-$n$ exchangeable (non-IID) predictive bounds (DE-T1(i) beyond the de Finetti regime) | The IID rates (DE-L3/L4) and the exchangeable qualitative statement are proved; finite-exchangeable sharpening tightens constants only; conservative reporting keeps every emitted bound valid. |
| Characterization of decision-robust (Tier-2) regions for **general non-separable joint losses** beyond ranking | Tier logic (DE-R6) and the general Γ-minimax machinery (DE-T4) are proved; per-loss geometric characterizations are conveniences. The operator computes robustness by direct test ($\tfrac12\notin$ interval, or argmin stability over $\mathcal Q$), which is always available. |
| Optimal choice of transport class (Λ vs TV vs $f$-divergence) for a given declared shift mechanism | The mandate requires declaring *some* class and knowing its price (DE-T2 does both); optimality among declarations is a modeling question outside the information ceiling, exactly analogous to the frozen theory's declared closure classes. |
| Infinite query sets / non-compact action spaces in DE-O2 | Compact-case existence is proved; extensions are standard l.s.c./coercivity refinements. All mandated decisions (finite $Q$, mandate V) are covered. |

None of these is a theorem whose absence prevents defining, certifying, auditing, or falsifying the decision operator. Each would tighten a bound, refine a taxonomy, or widen scope already handled by conservative reporting under H2/H6.

---

## Closing

Phase 7 leaves the identification theory exactly where Phases 0–6 froze it, and builds the one thing it provably could not contain: a selection principle. The separation is now itself a theorem chain — the set determines dominance and nothing finer (DE-S2/P1); selecting requires a completion (DE-P3); completions are either canonical-and-minimax or declared (DE-P6/S4); declarations can be learned from history exactly when a cross-member axiom is declared (DE-L, DE-T3), only tilt the inside of the identified set (DE-H2), and die gracefully into the frozen operator when withdrawn (DE-T4/O1). The decision layer selects; it never identifies. The ledger keeps the two forever distinct.

**Verdict: `DECISION_EXTENSION_COMPLETE`.**
