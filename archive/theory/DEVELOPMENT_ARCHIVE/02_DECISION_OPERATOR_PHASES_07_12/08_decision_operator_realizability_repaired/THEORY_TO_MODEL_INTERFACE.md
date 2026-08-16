# THEORY_TO_MODEL_INTERFACE — REPAIRED (Part VI, Final Model Compilation Contract)

> **Status:** Phase-8.1 terminal contract, 2026-08-03. Supersedes `../08_decision_operator_realizability/THEORY_TO_MODEL_INTERFACE.md` in full. Every clause cites its (repaired) theorem. Certificate types are load-bearing throughout: **floor** = exact-set or witness-certified lower bound; **guarantee** = outer-envelope upper bound; the two are never interchangeable (DR-F4-R).

---

## INPUT — what the model receives

1. **Archive** $w$: historical tasks $\{(D_i,y_i,c_i)\}$ — two typed channels, feasibility → $I_\theta$, frequency → $M_\phi$; disjoint by theorem (DE-H1) and by construction.
2. **Current task**: support $S_t$ ($k\le5$) and $\varepsilon$ as an argument.
3. **Decision specification $\gamma\in\mathcal G$** — explicitly including the query set $Q$, the pushforward map $g$, the action set $\mathcal A$ **with any abstain action and its declared cost**, the loss $L$, criterion, confidence $\delta$, tolerance $\eta$, tie-break $\tau$ (if single-valued output demanded).
4. **Declarations**: closure class; population rungs — (IID) or (C-IID) or a declared concentration condition (bare exchangeability buys no rate — DR-L2-R(i)); transport class; complexity declaration (finite target family for union bounds, or a declared uniform-convergence condition).

## LEARNED — what is allowed to be learned

1. **$I_\theta$**: an outer envelope $\widehat J\supseteq J_Q(O)$ (order projection $\widehat\Sigma\supseteq\Sigma$) under the declared closure class; optionally a certified **witness list** $\widetilde J\subseteq J_Q(O)$ (verified admissible members — the floor device, DR-F5-R). Objective: envelope tightness and witness richness (bracket width).
2. **$M_\phi(\cdot,\gamma)$**: the outer confidence class $\widehat{\mathcal Q}_\gamma$ — the DR-L3-R interval polytope with simultaneous coverage at the union-bound constant (or a declared uniform bound), fiber-relative under (C-IID).
3. **$D_\psi$**: criterion evaluation and $\eta$-argmin computation.
4. **Nothing else.** $\tau$, criterion, loss, abstention cost: declared, never learned.

## FORBIDDEN — what cannot be inferred (each a theorem or type rule)

1. Shrinking $I(O)$, $J_Q(O)$, $\Sigma$, $\Delta_{ab}$ by historical frequencies (DE-H2/H3); no $M_\phi\to I_\theta$ edge.
2. **Type-swapping certificates**: presenting an outer-envelope value $G_{\mathrm{cert}}(\widehat J)$ as an information floor (DR-F4-R: the $\{0\}$-vs-$\{0,100\}$ witness); presenting an inner witness set as a feasibility report; claiming an achievable unconditional worst case below any certified witness floor $R_{\mathrm{set}}(\widetilde J)$ (DR-F2 + DR-F5-R).
3. **Inner approximations as feasibility certificates** (DR-F4-R(c)); spuriously excluded orders in $\widehat\Sigma$ (false Tier-1 claims).
4. **Hidden measures**: undeclared tie-breaks at symmetric ties (DR-S3); posteriors without declared likelihood (DE-H4); undeclared reference weightings (DE-S5).
5. **Undeclared statistical strength**: rates from bare exchangeability (shared-Bernoulli witness, DR-L2-R); simultaneous coverage without the union/uniform declaration (DR-L3-R); cross-population transfer without a declared class (DE-T3); cross-fiber borrowing without a declared modulus (CI-A5(iii) lift).
6. **Loss-independent verdicts**: any abstention, off-coverage, or failure value not computed from the declared $(\mathcal A,L)$ (DR-F1-R(iii), T7); "optimal" where only "$\eta$-optimal" is certified (DR-S5); continuous single-valued output across a branch-switching tie (DR-S4-R).

## OUTPUT — what a prediction represents

$$\big(\ \mathcal A^*_\eta\ \text{or}\ \tau(\mathcal A^*_\eta)\ \text{or declared abstain}\ ;\ \ \mathrm{Ledger}\ \big)$$

The action output represents the criterion-optimal act under the declared information — never the true value, order, or member. Mandatory Ledger rows (mutually non-recoverable, DE-U1–U6 + DR-F4-R):
1. **identification row** (declaration-invariant, H1): $\widehat J$/envelope summary, flags;
2. **bracket row**: certified floor bound $R_{\mathrm{set}}(\widetilde J,\mathcal A,L)$ (or "no witnesses computed") $\big|$ certified guarantee $G_{\mathrm{cert}}(\widehat J,\mathcal A,L)$ — typed as lower/upper brackets of the uncomputable exact floor $R_{\mathrm{set}}(J_Q(O),\mathcal A,L)$;
3. **conditional rows** with axiom tags: class-risk of the selected action, regret over $\widehat{\mathcal Q}_\gamma$, ranking tier (1 identified / 2 decision-robust / 3 tie-broken) with the simultaneous $p$-intervals;
4. **tolerance row**: $\eta$ and bracket-width diagnostics (DR-M2-R);
5. **echo row**: every declaration consumed, including $\tau$, abstain cost, rungs, and complexity declaration (H6).

## FAILURE — when the model must abstain or flag (all loss-typed)

1. **Empty section** (unrealizable support under the declared class): misspecification flag; no action (frozen semantics).
2. **Off-coverage query**: report $G_{\mathrm{cert}}(\widehat J,\mathcal A,L)$ under the declared loss — $+\infty$ only for losses unbounded on the feasible set (absolute error on an unbounded section); finite for bounded losses ($\le1$ for $0$–$1$ ranking, $\le c$ with declared abstention). Abstain **iff** the certified guarantee exceeds the declared tolerance under $(\mathcal A,L)$ — never as a loss-independent reflex (T7).
3. **Symmetric tie with single-valued strict output demanded and no declared $\tau$**: return the action set; or the declared abstain action if its declared cost makes it criterion-optimal ($c\le$ the strict actions' guarantee, DR-F1-R(iii)); a strict pick here would be a hidden measure (DR-S3).
4. **Empty transported class** $\widehat{\mathcal Q}_\gamma=\emptyset$: declarations jointly inconsistent with history — surface the shift, do not silently widen.
5. **Validity self-check failure** (envelope containment or class coverage audit fails): certificates void; only flags and the DE-T4 minimax fallback under the declared $(\mathcal A,L)$ remain (DR-M1-R monotone degradation).
6. **Branch-switching tie proximity** with demanded continuous single-valued output: emit $\mathcal A^*_\eta$ or flag localized error (DR-S4-R); convex-prediction selections are exempt.

---

**Compilation criterion.** A system is a valid realization iff it factors as $D_\psi\circ(I_\theta\times M_\phi)$ with $M_\phi$ explicitly indexed by $\gamma$ (or a declared projectively-consistent family), each factor satisfies its repaired validity predicate ($V_I$, $V_M$, $V_D$), and its emissions instantiate this contract. Falsification protocols: the frozen suites (P1–P10, NP-1…6, T1–T5) plus the repaired Phase-8 checks — **type audits** (no outer value ever labeled "floor"; bracket rows correctly ordered $R_{\mathrm{set}}(\widetilde J)\le G_{\mathrm{cert}}$), witness-floor consistency (emitted guarantees dominate every certified witness floor), outer-envelope containment on constructed families (the DR-J3 reversal witness with its exact bracket $[2,3]$), no-hidden-measure symmetry tests (DR-S3/DE-S4), abstention-cost sensitivity (moving the declared $c$ across the strict-action guarantee must flip the abstain decision — a loss-typing probe), and the three-width ablations (fiber count, per-task coverage, transport radius must move exactly their own term, with the union-bound constant, DR-L4-R).
