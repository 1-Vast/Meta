# Decision-Operator Repair — Stopping Criterion

> **Status:** Phase-8.1 terminal decision, 2026-08-03. Sources: the repaired files in this directory, the audit `../09_phase8_audit/CONSOLIDATED_REPORT.md`, and the frozen Phases 0–7 (cited only, unmodified). This file replaces `../08_decision_operator_realizability/realizability_stopping_criterion.md` as the verdict of record; the original Phase-8 verdict `DECISION_OPERATOR_COMPLETE` is **withdrawn** (the audit was correct: DR-F4's certificate reading was a false central claim).

---

## The decision

$$\boxed{\textbf{DECISION\_OPERATOR\_REPAIRED}}$$

## Condition 1 — all audit failures corrected

| Audit failure (consolidated report §) | Correction | Where |
|---|---|---|
| §4 outer-envelope value used as information floor (counterexample $J=\{0\},\widehat J=\{0,100\}$) | Three-type discipline: exact floor / outer $G_{\mathrm{cert}}$ guarantee / inner witness floor bound; Phase-8 interpretive sentence retracted; bracket $[R_{\mathrm{set}}(\widetilde J),G_{\mathrm{cert}}(\widehat J)]$ is the honest emission | DR-F4-R, DR-F5-R |
| §4 abstention values loss-independent | All values from declared $(\mathcal A,L)$ incl. abstain cost: $\min(1,c)$, $\min(\tfrac12,c)$ | DR-F1-R(iii); contract Failure 2–3 |
| §3/§5 DR-S4 too broad | Narrowed to separated branch-switching ties, with the connectedness proof; connected-argmin variation exempt | DR-S4-R |
| §6 proxy value "floor of 3" | Re-typed: exact floor $2$, proxy guarantee $3$, witness pair closes the bracket $[2,3]$ | DR-J3(iii)-R |
| §7 $M_\phi$ untyped in $g,Q$ | Explicit index $\gamma\in\mathcal G$ in the domain, or universally indexed projectively consistent family with joint confidence accounting | DR-M-R §2 |
| §7/§4 $V_D$, DR-M1 built on the false floor reading | $V_D$ re-stated (guarantee validity + floor honesty + loss-typed fallback); DR-M1-R proved without the retracted step | DR-M-R §3–4 |
| §7 off-coverage universal-$\infty$ | Loss-typed: $G_{\mathrm{cert}}(\widehat J,\mathcal A,L)$, finite for bounded losses; abstain iff guarantee exceeds declared tolerance | contract Failure 2 |
| §8 exchangeability $\Rightarrow$ concentration | Rates only under (IID)/(C-IID)/declared concentration; shared-Bernoulli witness adopted; (EXCH) demoted to mixture statements | DR-L2-R |
| §8 no simultaneous listwise coverage | Union-bound theorem with $\eta_n=\sqrt{\ln(4|S|/\delta)/2n}$ over the declared finite family; DKW route for scalar CDFs; complexity declaration for infinite families | DR-L3-R |
| §8 over-broad "iff" | Necessity scoped to the frozen distribution-free information model; declared-likelihood aggregate identification acknowledged as a permitted extra declaration | DR-L3-R |

## Condition 2 — no new hidden assumptions

Every assumption added by the repair is surfaced in `REPAIR_LOG.md` and in the contract's echo row: witness certificates (elementary: monotonicity + verified membership), (IID)/(C-IID)/declared-concentration rungs, finite-index union bound or declared uniform-convergence condition. Nothing enters silently; each is a tag the output must carry.

## Condition 3 — the interface is mathematically closed

$I_\theta$ (outer envelope + optional witness list; $V_I$) $\times$ $M_\phi$ (explicitly $\gamma$-indexed outer class; $V_M$ with simultaneous coverage) $\to$ $D_\psi$ (guarantee-valid, floor-honest, loss-typed selection; $V_D$): the composition theorem DR-M1-R is proved **without** the retracted step — every emission valid under its stated type, degradation monotone and loss-typed; the tightness calculus DR-M2-R gives the auditable bracket-width accounting; the minimality guard stands (three typed objects, no legal merger). No emission type lacks a validity proof; no validity proof appeals to a retracted claim; a future model builder inherits only theorems, types, and declared slots — nothing to invent or repair.

## Residual open items (tightness only, non-blocking — carried from Phase 8 with types corrected)

Discrete-loss tightness rates near order boundaries (validity unaffected — both bracket ends remain valid); exact computability of $\Sigma(J)$ beyond the convex regime ($\Sigma^{\mathrm{pair}}$ remains the valid outer proxy); fiber-smoothing rates under a declared modulus; sharp constants beyond the union bound. Each affects bracket *width*, never emission *validity*.

---

**Verdict: `DECISION_OPERATOR_REPAIRED`.**
