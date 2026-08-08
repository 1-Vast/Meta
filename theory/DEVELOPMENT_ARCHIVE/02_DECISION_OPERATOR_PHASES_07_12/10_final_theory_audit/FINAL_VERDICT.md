# Phase 8.1 Repair Verification Audit: Consolidated Final Report

## Executive judgment

The repair successfully corrects the original central direction error: the exact-set minimax value is the deterministic information floor, while the value on an outer envelope is a conservative upper robust surrogate and never a lower floor. It also correctly adds explicit decision indexing to \(M_\phi\), replaces exchangeability-based concentration by IID/C-IID or a declared concentration condition, supplies a finite union bound, and makes off-coverage infinity loss-dependent.

The repaired interface is still invalid for model compilation. At least three load-bearing statements fail mathematically, and the composition theorem claims a conditioned population-law guarantee that its learning theorem does not establish.

## Repair matrix

| Required repair | Verdict |
|---|---|
| Outer envelope is upper surrogate, exact \(J\) retains the floor | `PASS` |
| Argmin set or declared \(\tau\), no hidden measure | `PASS` |
| Abstention consistently determined by \((\mathcal A,L)\) | `FAIL` |
| \(M_\phi\) explicitly indexed by context/query/decision specification | `PASS` |
| No concentration from bare exchangeability | `PASS` |
| Simultaneous finite-family event confidence | `PASS` |
| Complete joint-law object from those bounds | `FAIL` |
| Off-coverage infinity is loss-typed | `PASS` |
| DR-S4 scope correction is a valid theorem | `FAIL` |

## Blocking findings

### 1. The repaired selector-discontinuity theorem is false as stated

DR-S4-R concludes that endpoint minimizers in separated closed sets force every selector to jump. Its proof silently adds the condition that every intermediate argmin lies in the union of those two sets.

Let \(\mathcal A=[0,1]\), \(\rho_t(a)=(a-t)^2\), \(\mathcal A_0=\{0\}\), and \(\mathcal A_1=\{1\}\). The written hypotheses hold, but the unique selector \(s(t)=t\) is continuous. The missing union/separated-components hypothesis must be added to the theorem, not only to its proof.

### 2. The terminal abstention rule can select a strictly worse action

With both pairwise signs feasible, strict 0-1 actions have robust value 1. Give abstention cost \(c=2\) and tolerance \(T=1/2\). Then the repaired contract's \(G_{\rm cert}=1>T\), so Failure 2 mandates abstention, although abstention has loss 2 and is not criterion-optimal.

When no action meets tolerance, the contract must flag infeasibility. Abstention can be selected only when its declared loss/criterion makes it optimal, or under a separately declared refusal semantics. The current clause contradicts the repair's own statement that abstention is an ordinary action rather than a meta-rule.

### 3. Simultaneous event intervals are mis-typed as a joint law

The Hoeffding/union-bound calculation gives valid simultaneous intervals. But a probability law “on \(S\)” is not the right object when \(S\subsetneq S_m\): a one-element \(S=\{123\}\) would force probability one on 123. Overlapping pairwise events are not categorical outcomes either.

The confidence class must live on the full jointly realizable outcome space, or add a residual outcome and the required joint constraints. Until then the ranking population codomain is incomplete.

### 4. The composition theorem does not establish current-observation conditioning

DR-L3-R estimates transported population event frequencies. \(V_M\) and DR-M1-R then refer to the “true conditioned law,” but no repaired theorem conditions that law on current observations \(O\). A marginal law of \(g(f)\) does not determine \(P(g(f)\mid O)\).

The contract must explicitly choose one proved route: learn a joint pushforward of the conditioning statistics and decision target, declare a likelihood and conditioning theorem, or declare a context/fiber design under which DR-L3-R directly samples the current conditional population.

## Further required scoping

- \(R_{\rm set}\) is a deterministic floor. The assertion that “no rule” beats a witness value must use \(R_{\rm rand}\) for randomized policies.
- \(G_{\rm cert}\) is an achieved guarantee only under argmin existence. For an \(\eta\)-argmin realization, the guarantee is the selected action's outer risk or carries the appropriate \(+\eta\) slack.
- The archive feasibility channel must be formally invariant to task order and duplicate multiplicity so the declared set-of-traces channel cannot become a frequency channel in implementation.

## Existence, identification, learning

- **Existence:** conditional on compactness/lower-semicontinuity, or approximate with explicit slack; not unconditional as the floor file currently phrases attainment.
- **Identification:** the set-valued decision is determined by the exact joint object and declared decision specification; this part passes.
- **Learning:** simultaneous event bounds are proved under correct concentration assumptions, but the complete conditioned joint decision law and end-to-end certified realization are not.

## Model-compilation decision

A separate engineering agent cannot yet construct a trainable model without making theorem-level choices about conditioning, joint probability structure, abstention failure semantics, randomized versus deterministic floors, and approximate attainment. These are not selectable deployment assumptions inside an otherwise complete contract.

PROCESS_VERDICT: `THEORY_FREEZE_CONFIRMED`

MODEL_COMPILATION_VERDICT: `DECISION_OPERATOR_REPAIR_FAILED`
