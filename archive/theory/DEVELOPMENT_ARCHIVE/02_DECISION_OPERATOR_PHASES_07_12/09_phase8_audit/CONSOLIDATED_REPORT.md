# Phase 8 Final Decision Operator Audit: Consolidated Report

## Executive judgment

Phase 8 contains several mathematically sound components: it preserves the identification/decision separation, replaces the scalar half-diameter by a loss-typed exact minimax value, requires honest set-valued or explicitly tie-broken selection, and uses joint query/order objects rather than independent marginals.

The compilation claim nevertheless fails. DR-F4 reverses the certificate meaning of outer approximation. For \(J\subseteq\widehat J\), monotonicity gives

\[
R_{\mathrm{set}}(J,\mathcal A,L)
\le
R_{\mathrm{set}}(\widehat J,\mathcal A,L).
\]

The right side is an upper bound on the exact minimax value and a conservative surrogate robust value. It is not a larger valid lower information floor for the true problem. Phase 8 uses it as such in the validity predicate and terminal contract. This is a false central theorem, not a deployment assumption.

Additional theorem and typing defects occur in the abstention specialization, exchangeability-based rates, simultaneous listwise confidence, the domain of the population operator, and the off-coverage failure clause.

## 1. Theory freeze

**Result:** `THEORY_FREEZE_CONFIRMED`.

- All 26 Phase 0-6 files match the SHA-256 snapshot recorded by the preceding in-scope freeze audit; mismatches: `0`.
- Phase 7 was written from `2026-08-03T08:10:05Z` through `08:18:51Z`. Phase 8 begins at `08:37:39Z` and adds seven Markdown theory/interface files in its own directory.
- No Phase 0-7 file was overwritten by Phase 8.
- All explicit local Markdown links in Phase 8 resolve; broken links: `0`.
- Phase 8 adds theorem, interface, and stopping-criterion text only. No model architecture, experiment, or dataset contaminates the theory extension.

There is no independently signed pre-Phase-8 Phase-7 manifest in the allowed tree. The Phase-7 conclusion therefore rests on chronology and the current snapshot. No freeze violation is evidenced.

## 2. Identification layer

**Result:** pass.

The exact current-task objects remain

\[
\mathcal I(O),\qquad
J_Q(O)=\{(f(x_1),\ldots,f(x_m)):f\in\mathcal I(O)\}.
\]

The learned feasibility output \(\widehat J\) is expressly an outer representation under a declared closure class. Population learning produces a law class on a decision-relevant pushforward. The interface's hard no-feedback rule prevents this frequency channel from shrinking current-task feasibility.

Losses, criteria, preferences, and tie-breaks act only after identification. No historical frequency is promoted to a current-member observation. The later error about the certificate meaning of \(\widehat J\) occurs at the decision layer and does not redefine \(\mathcal I(O)\).

## 3. Decision layer

**Result:** pass with scoping corrections.

Phase 8 requires the loss and criterion to be declared. It permits a single-valued action only when uniqueness is proved or a tie-break \(\tau\) is explicit. An implicit full-support reference measure is correctly reclassified as a declared preference and otherwise forbidden.

The same identified set can and should yield different actions under different declared criteria. For \(J=\{0,1\}\) and \(\mathcal A=\mathbb R\):

- minimax absolute loss selects \(1/2\);
- Bayes squared loss with declared \(P(V=1)=0.9\) selects \(0.9\);
- Bayes 0-1 classification under the same declared law selects \(1\).

None changes \(J\). They differ only through declared decision information.

DR-S4 is too broad. A jump is forced for a branch-switching tie between separated components, such as the two discrete rankings. It is not forced along every continuous problem path that merely contains a nonunique argmin.

## 4. Robust information floor

**Result:** fail.

### Exact-set theorem

The definition

\[
R_{\mathrm{set}}(J,\mathcal A,L)
=\inf_{a\in\mathcal A}\sup_{v\in J}L(a,v)
\]

is the correct deterministic minimax value for the exact joint identified object. It is compatible with arbitrary declared loss and structured action spaces. Its randomized version is likewise correctly typed.

For scalar absolute loss it gives

\[
R_{\mathrm{set}}(J,\mathbb R,|a-v|)
=\tfrac12(\sup J-\inf J),
\]

so the frozen absolute-error theorem is recovered rather than universalized. For pairwise 0-1 ranking with both signs feasible and no abstention, the deterministic and randomized values are \(1\) and \(1/2\).

### Outer-envelope counterexample

Take \(J=\{0\}\), \(\widehat J=\{0,100\}\), \(\mathcal A=\mathbb R\), and \(L(a,v)=|a-v|\). Then

\[
R_{\mathrm{set}}(J)=0,
\qquad
R_{\mathrm{set}}(\widehat J)=50.
\]

Action \(a=0\) has true worst-case loss zero. Hence 50 cannot be a valid lower floor on achievable risk for the true problem. DR-F4(i)'s monotonicity calculation is correct, but the sentence interpreting the larger number as a valid conservative information floor is false.

The proper type discipline is:

- \(R_{\mathrm{set}}(J)\): exact minimax value and true information lower bound;
- \(R_{\mathrm{set}}(\widehat J)\) for outer \(\widehat J\): upper bound on the true minimax value and conservative robust surrogate;
- \(R_{\mathrm{set}}(\widetilde J)\) for inner \(\widetilde J\): possibly weak lower bound, although \(\widetilde J\) is not a valid feasibility certificate.

The distinction between a feasible-set certificate and a lower-bound certificate is essential. Phase 8 conflates them.

### Abstention counterexample

DR-F1(iii) asserts the pairwise values \(1\) and \(1/2\) “with or without abstention.” If abstention is available with declared constant loss \(c<1/2\), both minimax values are at most \(c\). There is no loss-independent abstention result; the abstain action and its loss must be included in \((\mathcal A,L)\).

## 5. Selection operator

**Result:** pass with DR-S4 narrowed.

The honest outputs are exactly the appropriate two forms:

\[
\mathcal A^*=\operatorname*{argmin}_{a\in\mathcal A}\rho(a)
\]

or \(\tau(\mathcal A^*)\) for an explicit admissible \(\tau\). The symmetry proof for two freely swapped rankings correctly rules out a canonical equivariant singleton. A hidden measure, prior, or preference is therefore neither required nor permitted.

The approximate-risk result is valid: uniform risk error at most \(\eta/2\) makes a minimizer of the approximate risk \(\eta\)-optimal for the true risk. Set-valued \(\eta\)-argmins are an honest realization device.

## 6. Joint decision object

**Result:** pass.

Phase 8 correctly distinguishes marginals, pairwise differences, jointly realizable orders, and the full joint value set.

For a direct counterexample, let

\[
J_{\mathrm{diag}}=\{(0,0),(1,1)\},\qquad
J_{\mathrm{anti}}=\{(0,1),(1,0)\}.
\]

Both have marginal intervals \([0,1]\) in each coordinate. The diagonal object always ties the two queries; the anti-diagonal object admits either strict order and never ties them. Independent marginals therefore do not determine ranking behavior.

The DR-J3 reversal witness also correctly proves that separately realizable pairwise signs can over-admit listwise orders. Its exact order set has Kendall minimax value 2, while the all-orders outer proxy has value 3. The latter is a conservative surrogate robust value, not a true information floor of 3.

## 7. Meta-learning interface

**Result:** fail as a complete compilation contract.

The three roles are conceptually well separated:

1. \(I_\theta\): outer feasibility representation under declared closure assumptions;
2. \(M_\phi\): outer confidence/ambiguity class for population decision information;
3. \(D_\psi\): declared criterion evaluation and honest selection/abstention with a ledger.

The contract is not closed for three reasons.

First, \(V_D\) and DR-M1 invoke the false proposition that the minimax value on outer \(\widehat J\) is a valid lower information floor. The claimed end-to-end validity theorem therefore does not follow.

Second, the stated map

\[
M_\phi:\mathcal W\times\mathrm{tags}\to\mathfrak Q_g\times(0,1]
\]

does not take \(g\), \(Q\), or an equivalent decision context as input, although its codomain is a class of laws on the corresponding \(g\)-space. The index must be an explicit input, or the output must be a universally indexed family.

Third, the failure contract declares \(R_{\mathrm{set}}=\infty\) for off-coverage queries. Unbounded value feasibility does not force infinite risk under bounded 0-1, ranking, or abstention losses. The fallback must remain loss-typed.

## 8. Learnability

**Result:** incomplete and partly incorrect.

The three required notions separate as follows.

- **Existence:** established under the compactness/lower-semicontinuity assumptions for the argmin, or through a declared approximate action set.
- **Identification:** set-wise from \(J_Q(O)\), the declared decision context, and the population ambiguity object; a canonical singleton is intentionally absent at undeclared ties.
- **Finite-history estimation:** DR-L3 gives a real forced/compatible interval construction for pairwise signs and listwise orders, but its advertised coverage theorem is not valid as stated.

Mere exchangeability does not yield Hoeffding or DKW concentration. If all task variables equal a shared Bernoulli \(Z\), the sequence is exchangeable but the empirical frequency is always 0 or 1 and does not concentrate around its marginal probability \(1/2\). The theorem needs IID, conditional IID within fibers, or a separately declared concentration condition.

The interval radius shown in DR-L3 controls the two bounds for a fixed order. Intersecting those intervals across \(m!\) orders does not automatically give simultaneous \(1-\delta\) coverage. A union allocation such as \(\delta/m!\), or a proved uniform multinomial bound, is required.

Finally, the claimed “if and only if” is too broad. Per-task set identification is sufficient for this distribution-free forced/compatible construction. It is not necessary over all declared statistical models, because a known likelihood or measurement-error model may identify a population aggregate without identifying every historical latent value. Necessity must be scoped to the frozen distribution-free information model.

## 9. Required corrections before compilation

The smallest theorem-level repair is:

1. retype outer-envelope minimax values as surrogate upper robust values, while retaining the exact-set value as the information floor;
2. recompute every abstention result from the declared abstention loss;
3. restate finite-history rates under IID/conditional IID or another explicit concentration assumption;
4. prove simultaneous coverage for the listwise law polytope;
5. index \(M_\phi\) by \(g\), \(Q\), or decision context;
6. make off-coverage behavior depend on \((\mathcal A,L)\);
7. scope selector discontinuity to separated branch-switching ties.

Until these corrections are proved and propagated through \(V_D\), DR-M1, the compilation contract, and the stopping criterion, a future model builder would have to invent or repair mathematics.

## Final classification

PROCESS_VERDICT: `THEORY_FREEZE_CONFIRMED`

MODEL_COMPILATION_VERDICT: `DECISION_OPERATOR_INVALID`
