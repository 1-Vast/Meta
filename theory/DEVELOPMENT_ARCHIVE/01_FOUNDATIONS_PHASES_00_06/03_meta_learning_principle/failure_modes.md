# Failure Modes — Where the Principle Cannot Work

> **Status:** Phase-3 catalogue, 2026-08-02. Every mode below is *certified*: it carries a reference to a theorem or counterexample of the corpus (F-, C-, CP-, OP-, MP-numbers), and an **observable signature** — what a system that has hit the mode will exhibit. The catalogue is exhaustive relative to the corpus; modes are grouped by the layer at which the principle breaks. All statements are worst-case (the standing scope of the theory).

Format per mode: *Description* → *Certificate* → *Signature*.

---

## A. Structural absence — nothing is shared

**A1. No shared structure at all.** The task family imposes no constraint linking support to query (trace sets are full cylinders). Meta-learning is vacuous: minimax error $+\infty$ at every off-support query; information exists exactly on the support.
*Certificate:* F4 / C1. *Signature:* performance indistinguishable from a rule that never saw the archive; error unbounded under the F4 two-member protocol.

**A2. Unbounded prior with no binding constraint at the query.** Even with genuine shared structure, a specific query can sit where the family imposes nothing (section $=\mathbb R$): the center is undefined, radius $+\infty$.
*Certificate:* F6 ($\phi(x)\notin\operatorname{row}(G)$); §2.3 of `operator_specification.md`. *Signature:* any emitted point-value at such a query is fabrication; adversarial members realize arbitrary values with identical support data.

---

## B. Support-configuration failures — the evidence is placed wrongly

**B1. Rank failure at the actual support.** Identifiability is a joint property of family and configuration; the multivariate world guarantees bad configurations exist: on any triod-containing domain, every $d\ge2$-dimensional space of continuous functions has singular size-$d$ designs.
*Certificate:* F10(ii) (Mairhuber–Curtis–Sieklucki, size-$d$ scoping); F6. *Signature:* abrupt loss of adaptation quality under small support-placement changes; rank drop of the induced evaluation map. *Scope:* with $k>d$ points design-independence can return (F10(iii) triod example) — the mode is about size-$d$ configurations.

**B2. Support size below the task dimension.** For families containing a continuously parametrized $d$-cell, $k<d$ makes exact recovery impossible for every configuration.
*Certificate:* F14(i). *Scope (refereed):* constrained coefficient sets can be identifiable below full rank ($\ker G\cap(C-C)=\{0\}$ with $\operatorname{rank}G<d$, F6(iii)) — the mode applies to $d$-cell families. *Signature:* residual ambiguity of dimension $\ge d-k$ no matter where the support is placed.

**B3. Support size between $d$ and $2d$: global ambiguity at every configuration.** Nonlinear families can defeat *every* support of size $k=d$ (and more): for $\sin(\theta x)$ tasks, $k=1$ and $k=2$ fail at every configuration — the explicit collision pair $\theta=\frac{\pi}{2b}-\frac{\pi}{a}$, $\theta'=\frac{\pi}{2b}+\frac{\pi}{a}$ produces identical supports with different query values; $k=3=2d+1$ with pairwise-irrational ratios resolves it.
*Certificate:* F14(ii) / C4; protocol P7. *Signature:* bimodal ambiguity — two well-separated task hypotheses fit the support exactly; a support-size cliff at $2d+1$.

**B4. Non-quasianalytic task families: flat directions.** For merely smooth (non-quasianalytic) families, agreement sets can have nonempty interior: *positive-measure* sets of support configurations fail, for every $k$, and no genericity statement of any kind holds.
*Certificate:* F14(iii) / C5 ($\theta\,e^{-1/x}\mathbf1_{x>0}$). *Signature:* entire regions of support placements carry zero task information; random placement does not escape (the failure set is not null).

**B5. Saturation by noise.** When $\omega_{x,D}(2\varepsilon)=\operatorname{diam}T_x$ (bounded-prior regime), the support provably adds nothing at $x$ in the worst case: adaptation gain zero.
*Certificate:* MP-2. *Signature:* posterior certificate equals prior certificate; any apparent adaptation benefit at such queries does not survive the adversarial protocol.

---

## C. Coverage and archive failures — the meta-level evidence is insufficient

**C1. Off-coverage queries.** At locations never observed in any previous task, family structure is unknowable: within the linear exactly-$d$ class the consistent value set is all of $\mathbb R$ (or forced-zero $\{0\}$, a degenerate exception that vanishes for $\varepsilon>0$). No volume of archive at *other* locations helps.
*Certificate:* F18 / C13; protocol P6. *Signature:* the F18 modification changes the truth arbitrarily while every observable is unchanged — claimed off-coverage accuracy is prior smuggling; the smuggled member-level assumption (continuity, analyticity, kernel) must be declared to be legitimate.

**C2. Task-diversity rank failure.** Fewer than $d$ sufficiently diverse previous tasks ($\operatorname{rank}A<d$): the window is identified only up to an affine ambiguity of dimension $d-\operatorname{rank}A$, and *every* point of the ambiguity is realized by a genuine candidate family.
*Certificate:* F17 (necessity direction). *Signature:* consistent multi-model disagreement that no additional per-task data resolves; $n\ge d$ tasks are necessary, and a **single-task archive** with $d\ge2$ leaves ambiguity of dimension $\ge d-1$.

**C3. Unknown or misspecified task dimension.** Archive rank certifies $\dim\ge d$ only, never $\le d$: choosing the class dimension too small biases every window; too large re-opens B2-type ambiguity. The identification theorem is an iff *within the class of exactly-$d$-dimensional families with $d$ known*.
*Certificate:* F17 scoping; `theorem_summary.md` §8. *Signature:* systematic residuals at core configurations (under-dimensioned) or unresolvable coset ambiguity (over-dimensioned).

**C4. Noisy archive.** With inexact archive observations, only the factorization of the total operator and a first-order subspace-perturbation expectation are established; the quantitative constant is **open**.
*Certificate:* `operator_formulation.md` §6 (open flag). *Signature:* certificates computed from noisy-archive windows are not yet backed by a theorem; the interface requires the flag to propagate.

**C5. Consistent-but-unrealizable meta-objects.** A learned window system can satisfy every projective-consistency check and still correspond to no family at all (uncountable coverage, non-compact windows).
*Certificate:* MP-1(iii) (Waterhouse). *Signature:* a "family model" passing all pairwise/finite consistency audits while its global object is empty; countable coverage or compact (closed, bounded) windows restore realizability (MP-1(i),(ii)).

---

## D. Stability failures — identifiable, yet unusable

**D1. Identifiable but noise-fragile families.** Exact adaptation possible at $\varepsilon=0$, worthless at every $\varepsilon>0$: $\omega(0)=0$ with $\omega(t)=\infty$ for all $t>0$.
*Certificate:* F3 / C3 (tanh family); protocol P8. *Signature:* performance cliff between clean and noisy evaluation, error diverging along the family's saturating direction.

**D2. Misspecified noise level.** The optimal rule depends on $\varepsilon$ (F1 Rem. 1.4). Under-reported $\varepsilon$: realizability failures (empty sections) on genuine data — the operator's misspecification detector fires on valid inputs, or, if silenced, certificates become false. Over-reported $\varepsilon$: inflated radii, adaptation gain destroyed as far as saturation (B5).
*Certificate:* F1 Rem. 1.4; MP-2; §2 of the specification. *Signature:* asymmetric — spurious "inconsistent data" flags (under) versus systematically vacuous certificates (over).

**D3. Discontinuity of optimal adaptation.** The true operator jumps as support values cross realizability boundaries of competing members; continuous surrogates carry irreducible sup-error $\ge$ half the jump, localized at the transitions.
*Certificate:* MP-4; protocol P9. *Signature:* error spikes concentrated on narrow bands of support values; unremovable by capacity increases while the surrogate stays continuous.

**D4. Misspecified data (model error).** Support inconsistent with every member at level $\varepsilon$ (distance $\eta>0$ to the realizable set): the operator is undefined; the recorded projection convention costs the doubling $\tfrac12\omega(2\varepsilon+2\eta)$ plus the family's approximation error at the query.
*Certificate:* treatise §10.4; OP definitions. *Signature:* the realizability check fails — which is a feature: with $k>d$ consistent overdetermination, the same check is the misspecification detector.

---

## E. Representation failures — the demanded object does not exist

**E1. No finite continuous task state.** Infinite-dimensional families (the 1-Lipschitz class) admit *no* finite-dimensional continuous summary: any architecture-agnostic "task vector" of fixed finite dimension provably cannot represent the task. The principle still operates — via sections — but only with interval outputs and query-coupled information.
*Certificate:* F19(v) / C7. *Signature:* task-vector approaches plateau at an error floor equal to the width the discarded information carries; the floor moves with the query (min-distance geometry).

**E2. Dimension excess of continuous representations.** Even finite-dimensional task spaces can require more coordinates than their dimension: tripod-parametrized families need $m_{\min}=2>d=1$; only the sandwich $d\le m_{\min}\le2d+1$ is guaranteed.
*Certificate:* F19(ii),(iii) / C6. *Signature:* topological obstructions — no continuous injective $d$-dimensional representation exists although the family "is $d$-dimensional".

**E3. Residual gauge at small $k$.** With too few support points, a symmetry of the family survives: only invariants (differences, ratios, cross-ratios — per the gauge taxonomy's exact observation counts) are identifiable; absolute values output at that stage are fabricated.
*Certificate:* F16; F15. *Signature:* systems agree across runs on relative quantities and disagree arbitrarily on absolute ones; the disagreement is exactly along the stabilizer orbit.

**E4. Dimension claims without continuity.** Any claim that a $1$-dimensional (or any fixed-dimensional) code "suffices" is vacuously true measurably and false continuously; dimension talk is meaningful only in the continuous/stable category.
*Certificate:* F19(iv) / C8 (Borel collapse). *Signature:* compression results that do not survive a stability/continuity audit.

---

## F. Evaluation and aggregation pitfalls — misreading the theory

**F1. Scalar design scores.** Two supports with equal radius at a query need not be informationally equivalent; the complete comparison is the section-containment (Blackwell-type) preorder, of which the radius is a strictly coarser scalar quotient.
*Certificate:* C12; treatise §6. *Signature:* support-selection procedures optimizing a single scalar mis-rank supports that the preorder distinguishes.

**F2. Adaptive-versus-nonadaptive gaps.** For convex balanced families, adaptive support selection gains nothing (so adaptive machinery is provably wasted there); for nonconvex (multi-branch) families it can gain unboundedly (so nonadaptive evaluation understates achievable performance).
*Certificate:* C14; treatise §6. *Signature:* regime-dependent reversal of the adaptivity ablation.

**F3. Average-case refutations of worst-case claims.** Every floor and ceiling of the theory is minimax. A benign task distribution can beat the radius floor on average without contradicting anything. Predictions are falsifiable only under the adversarial protocols (P1–P10); benchmark averages neither confirm nor refute the theory.
*Certificate:* scope declaration, `meta_learning_abstraction.md` §0. *Signature:* apparent "violations" that vanish under the two-member midpoint-data protocol.

**F4. Multi-query loss mismatch.** The exact constant $\tfrac12$ is scalar-target; only the sup-loss multi-query bridge is exact. Joint losses (e.g. $\ell_2$ over queries) carry radius-versus-diameter slack; evaluating the canonical per-query rule under a joint loss and penalizing the slack misattributes a theory property to the system.
*Certificate:* F1 Rem. 1.2 / Rem. 1.1. *Signature:* per-query optimal systems appearing suboptimal under joint metrics by exactly the radius/diameter gap.

**F5. Monotone-update priors.** Enforcing that more evidence moves estimates monotonically toward anything contradicts optimal adaptation: guarantees are monotone, estimates are not.
*Certificate:* OP-8; protocol P10. *Signature:* monotonicity-constrained systems deviate from conditional minimaxity by a computable margin on the nested-section protocol.

---

## G. Degenerate cases (conventions, not pathologies)

| Case | Status | Certificate |
|---|---|---|
| $k=0$ (no support) | Operator $=$ baseline; minimax $=\tfrac12\operatorname{diam}T_x$; only gauge-invariants of F16-type families are task-adapted "for free" | MP-2; F16 |
| Query in support ($x=x_i$) | Floor is $\varepsilon$, not $0$, for rich families — memorization has a worst-case price under noise | F4 |
| Singleton family | Everything identifiable at $k=0$; $\Delta\equiv0$; pure-baseline regime | trivial from F2 |
| Empty family | Excluded by convention (SA1): no realizable data, operator nowhere defined | SA1 |
| Repeated/coincident support points | Rank drop of the evaluation map; effective $k$ decreases | F6 |
| Overdetermined consistent support ($k>d$) | Redundancy $=$ free misspecification detector (realizability check) | D4 above |
| Forced-zero archive case | The single exception to off-coverage unconstrainedness; vanishes at $\varepsilon>0$ | F18 / C13 |

---

## Closing statement

None of these modes is an implementation defect; each is a theorem about the problem. A future system does not avoid them by architecture or scale — it avoids them by (i) operating inside the identifiable region, (ii) declaring the assumptions that extend that region, and (iii) shipping the radius with every prediction so that the modes, when hit, are visible rather than silent.
