# Phase 8.2 Closure Log

> **Status:** Phase-8.2, 2026-08-03. This directory, together with the unrefuted parts of `../08_decision_operator_realizability_repaired`, is the compilation source of record. Phases 0–7 untouched. Audit of record: `../10_final_theory_audit/` (verdict `DECISION_OPERATOR_REPAIR_FAILED`, four blocking findings + three scoping corrections). New results carry **DC-** numbers. No new architectures; no networks; re-statements, proofs, and explicitly declared assumptions only.

---

## Audit finding → closure map

| # | Blocking finding (audit) | Closure | File |
|---|---|---|---|
| 1 | **DR-S4-R false as stated**: the union-of-branches hypothesis lived in the proof, not the statement; witness $\rho_t(a)=(a-t)^2$ on $[0,1]$ has continuous selector $s(t)=t$ | **DC-S1**: corrected theorem with the union hypothesis explicit; audit's witness adopted (DC-S2) + ramp witness (DC-S3) showing endpoint separation and connected $\mathcal A$ never suffice; positive side: Berge/Michael sufficiency for continuous selectors, and the unconditional discrete-action jump lemma that preserves the ranking warning | `selector_theorem_repair.md` |
| 2 | **Abstention clause selects strictly worse actions**: with strict robust value $1$, abstain cost $c=2$, tolerance $T=\tfrac12$, the 8.1 contract mandated abstention (loss $2$) | **DC-A**: abstention is an action, selected **iff criterion-optimal** (threshold theorem: abstain iff $\rho(a_{\mathrm{abs}})\le R_{\mathrm{set}}(J,\mathcal A_{\mathrm{strict}},L)$ — the threshold is *derived*, and equals the declared cost compared to the strict robust value, nothing else); tolerance violation is **infeasibility**, reported as failure, never converted into abstention; the audit's triple $(1,2,\tfrac12)$ is the defect witness of record | `abstention_semantics_repair.md` |
| 3 | **Simultaneous intervals mis-typed as a joint law** (a "law on $S\subsetneq S_m$" forces mass 1 on $S$; overlapping pairwise events are not outcomes) | **DC-R**: the population ranking object is a constraint class of laws on the **full** outcome space $\Omega_m=S_m$: $\{P\in\Delta(\Omega_m):\ell(E)\le P(E)\le u(E),\ E\in\mathcal E\}$ with events (orders, pairwise, top-$k$) as *constraints*, not outcomes; pairwise probabilities typed as marginals; linear-ordering-polytope consistency; decomposable-loss sufficiency theorem + equal-marginal witness for exact-match loss | `joint_ranking_object_repair.md` |
| 4 | **Current-observation conditioning unproved**: DR-L3-R estimates population marginals, DR-M1-R used them as $P(\cdot\mid O)$ | **DC-C**: route **B — declared conditional fiber construction** (with the sufficiency declaration explicit): a declared context map $\kappa$ computable from any task record; (C-IID-$\kappa$); the conditioning-ladder theorem giving exactly what each declaration buys (nothing → support-only class = frozen minimax endpoint; SUFF-$\kappa$ → fiber class; LIK → single posterior), so no conditional claim exceeds its declared assumptions | `conditional_population_learnability.md` |

## Scoping corrections (audit "further required scoping") — absorbed

| Correction | Where |
|---|---|
| "No rule beats the witness floor" must use $R_{\mathrm{rand}}(\widetilde J)$ for randomized policies ($R_{\mathrm{set}}$ for deterministic only) | `final_interface.md` §Ledger typing |
| $G_{\mathrm{cert}}$ is *achieved* only under argmin existence; the unconditionally valid guarantee row is the **selected action's own outer risk** $\sup_{v\in\widehat J}L(\hat a,v)$ (with $\eta$-slack claims separate) | `final_interface.md` §Ledger typing |
| Archive feasibility channel formally invariant to task order and duplicate multiplicity (set-of-traces semantics enforced as an axiom, so the feasibility channel cannot become a frequency channel) | `final_interface.md` §$V_I$ |

## Carried unrefuted (audit: PASS)

Exact-set floor vs outer-surrogate typing (DR-F4-R bracket); argmin-set / declared-$\tau$ selection with no hidden measures; $M_\phi$ decision-specification indexing; IID/C-IID/declared-concentration rungs; simultaneous finite-family confidence (as *event intervals* — now correctly consumed by DC-R as constraints); loss-typed off-coverage values.
