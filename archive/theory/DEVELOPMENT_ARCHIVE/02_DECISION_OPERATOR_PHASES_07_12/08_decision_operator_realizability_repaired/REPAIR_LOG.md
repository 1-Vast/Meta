# Phase 8.1 Repair Log

> **Status:** Phase-8.1, 2026-08-03. This directory supersedes `../08_decision_operator_realizability` as the compilation source. Phases 0–7 untouched. Audit of record: `../09_phase8_audit/CONSOLIDATED_REPORT.md` (verdict `DECISION_OPERATOR_INVALID`). Repaired results carry the original **DR-x** numbers with suffix **-R**; unrepaired results are carried over verbatim by citation. No new theory is created: every repair is a re-typing, a re-scoping, or an elementary consequence (monotonicity, Hoeffding + union bound) of frozen material.

---

## Audit finding → repair map

| # | Audit finding | Repair | File |
|---|---|---|---|
| T1 | **DR-F4 reversed the certificate meaning of outer envelopes** ($R_{\mathrm{set}}(\widehat J)\ge R_{\mathrm{set}}(J)$ is an *upper* surrogate, not a valid floor; counterexample $J=\{0\},\widehat J=\{0,100\}$: true achievable worst case $0$, claimed "floor" $50$) | Three-type discipline (DR-F4-R): exact $R_{\mathrm{set}}(J)$ = information floor; outer $R_{\mathrm{set}}(\widehat J)$ = **certified guarantee value** $G_{\mathrm{cert}}$ (upper bound on achievable and on the floor); inner witness set $\widetilde J\subseteq J$ = **certified floor lower bound** $R_{\mathrm{set}}(\widetilde J)$. Honest report = the bracket $[R_{\mathrm{set}}(\widetilde J),\,G_{\mathrm{cert}}(\widehat J)]$. Propagated through $V_D$, DR-M1-R, contract, stopping criterion. | `loss_typed_information_floor.md` + all downstream |
| T2 | Hidden-measure risk in selection (already largely correct; DE-O2's $\mu_0$) | Carried over: Option A (argmin set) / Option B (declared $\tau$) exhaustive; DR-S3 unchanged (audit: pass). | `honest_selection_operator.md` |
| T3 | **Abstention constants loss-independent** (DR-F1(iii) claimed pairwise values $1,\tfrac12$ "with or without abstention") | DR-F1-R(iii): all values computed from declared $(\mathcal A,L)$ *including* the abstain action and its declared cost $c$: pairwise values $\min(1,c)$ deterministic, $\min(\tfrac12,c)$ randomized; abstention is an action, never a meta-rule. Propagated to failure clauses. | floor + contract files |
| T4 | **$M_\phi$ untyped in $g,Q$** (codomain referenced a $g$-space absent from the domain) | DR-M-R: $M_\phi$ takes the decision specification $(g,Q,\text{context})$ as explicit input — or outputs a universally indexed, projectively consistent family $\{\widehat{\mathcal Q}_{g,Q}\}$. | `meta_learning_interface.md` |
| T5 | **Exchangeability alone claimed to yield concentration** (audit counterexample: all tasks equal a shared Bernoulli $Z$ — exchangeable, empirical frequency never concentrates at the marginal $\tfrac12$) | DR-L2-R/L3-R: rates restated under declared **(IID)** or **(C-IID)** (conditional IID within fibers) or a separately declared concentration condition; exchangeability retained only as the *mixture* statement (targets the directing-measure-conditional law, a different object) and explicitly marked insufficient for rates. | `decision_information_learnability.md` |
| T6 | **No simultaneous coverage for the listwise polytope**; "iff" over-broad | DR-L3-R: union-bound allocation $\delta/(2|S|)$ over the queried order family $S\subseteq S_m$, $\eta_n=\sqrt{\ln(4|S|/\delta)/2n}$, with the finite-index complexity assumption declared (infinite families require a declared uniform-convergence condition); necessity of the per-task gate **scoped to the frozen distribution-free information model** (a declared likelihood model can identify aggregates without per-task identification). | same |
| T7 | **Off-coverage failure clause universal-$\infty$** (unbounded value set does not force infinite risk under bounded losses) | Failure clauses re-typed: off-coverage reports $G_{\mathrm{cert}}(\widehat J,\mathcal A,L)$ — infinite only for losses unbounded on the feasible set; finite (≤1, ≤c) for 0–1/ranking/abstention losses; abstain iff the certified guarantee exceeds declared tolerance under the declared $(\mathcal A,L)$. | contract |
| — | DR-S4 too broad (jump not forced at every non-unique argmin) | DR-S4-R: discontinuity forced exactly at **branch-switching ties between separated components** of $\mathcal A$ (e.g. discrete rankings); connected argmin variation excluded from the claim. | `honest_selection_operator.md` |
| — | DR-J3(iii) called the pairwise proxy value "a floor of 3" | Sentence re-typed: $R^{\mathrm{pair}}_{\mathrm{set}}=3$ is a conservative **surrogate guarantee**, the true floor is $2$. | `joint_query_decision_objects.md` |

## Unchanged (audit: pass)

Identification layer typing (outer $\widehat J$, no-feedback rule); the exact-set definition of $R_{\mathrm{set}}$ and the floor theorem DR-F2 *for the true set*; selection dichotomy DR-S1–S3, DR-S5; the joint-object lattice DR-J1–J4 (minus the one sentence above); honesty axioms H1–H6; the three-operator separation as concept.

## New declared assumptions introduced by the repair (none hidden)

(i) **Witness certificates** for floor lower bounds: finitely many explicitly constructed members verified admissible under the declared class (an elementary inner device, monotonicity only). (ii) **(IID)/(C-IID)** replacing (EXCH) wherever a rate is claimed. (iii) **Finite-index union bound** (or a declared uniform-convergence condition) for simultaneous coverage. Each appears in the axiom tags of the contract; nothing else was added.
