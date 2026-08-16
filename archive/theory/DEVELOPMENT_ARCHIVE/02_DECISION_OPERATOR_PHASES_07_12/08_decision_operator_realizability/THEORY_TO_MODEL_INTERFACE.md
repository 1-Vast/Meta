# THEORY_TO_MODEL_INTERFACE (Part VI — Final Model Compilation Contract)

> **Status:** Phase-8 terminal contract, 2026-08-03. This document is the **complete compilation target**: a future meta-learning system is correct relative to this program iff it implements the three-operator interface (`meta_learning_interface.md`) under the clauses below. Every clause cites its theorem; nothing here is advisory. Phases 0–7 frozen; Phase-8 results DR-F/S/J/L/M.

---

## INPUT — what the model receives

1. **Archive** $w\in\mathcal W$: historical tasks $\{(D_i,y_i,c_i)\}_{i\le n}$ — designs, observed values, optional auxiliary labels. Consumed through exactly two typed channels: feasibility (→ $I_\theta$) and frequency (→ $M_\phi$); the channels are disjoint by theorem (DE-H1) and must be disjoint in the implementation.
2. **Current task** $o\in\mathcal O$: support $S_t=\{(x_i,\tilde y_i)\}_{i\le k}$, $k\le5$, and the noise level $\varepsilon$ **as an argument** (frozen F1 Rem. 1.4).
3. **Query set** $Q\in\mathcal X^m$ (finite; ranking decisions consume $Q$ jointly, never per-query — DR-J).
4. **Declarations** (all explicit, all echoed in the output): family closure class (legitimizes outer envelopes — frozen); population axioms EXCH/C-EXCH + transport class + confidence $\delta$ (DR-L2); decision context $(\mathcal A,L,\text{criterion})$; tie-break $\tau$ if single-valued output is demanded (DR-S3); tolerance $\eta$ (DR-S5).

## LEARNED — what is allowed to be learned

1. **$I_\theta$ (feasibility):** an **outer envelope** of the window system / joint pushforward $J_Q(O)$ (and its order projection $\Sigma$) under the declared closure class. Objective: tightness. Constraint: $\widehat J\supseteq J_Q(O)$, $\widehat\Sigma\supseteq\Sigma$ ($V_I$; DR-F4, DR-J outer semantics).
2. **$M_\phi$ (frequency):** an **outer confidence class** $\widehat{\mathcal Q}$ of population laws on the decision-relevant pushforward — concretely the DR-L3 interval class (forced/compatible frequencies ± Hoeffding ± transport), at fiber-relative effective count under C-EXCH. Constraint: $V_M$ coverage.
3. **$D_\psi$ (criterion computation):** the conditioning, robust-risk evaluation, and $\eta$-argmin computation. Constraint: $V_D$ = H1–H6 + floor consistency + honest selection.
4. **Nothing else.** In particular the tie-break, the criterion, and the loss are *declared inputs*, never learned; learning them from data would be an undeclared preference (DE-P/DE-S5).

## FORBIDDEN — what cannot be inferred (each prohibition is a theorem)

1. Shrinking $I(O)$, $J_Q(O)$, $\Sigma$, or any $\Delta_{ab}$ by historical **frequencies** (DE-H2/H3); architecturally: no $M_\phi\to I_\theta$ edge (DR-M impossibility guard).
2. Any **unconditional** guarantee below the loss-typed floor $R_{\mathrm{set}}(J_Q(O),\mathcal A,L)$ — equivalently below the frozen radius in the scalar case (DR-F2/F3, DE-L5, DE-O4).
3. **Inner approximations** presented as certificates — of sets, of order sets, or of population classes (DR-F4(iii), DR-J outer semantics, DR-L2(iv)).
4. **Hidden measures**: any single-valued output at a tie without a declared $\tau$; any posterior without declared LIK; any "uniform prior" without a declared invariance class (DR-S3, DE-H4, DE-S5).
5. **Cross-population transfer without a declared transport class**; cross-fiber borrowing without a declared modulus on $C$ (DE-T3, DR-L2(ii)/CI-A5(iii)).
6. Off-coverage extrapolation without a declared member-level class (frozen F18); Tier-2/3 ranking outputs presented as identified (DE-R6); claims of "optimal" where only $\eta$-optimal is certified (DR-S5).

## OUTPUT — what a prediction represents

$$\big(\ \mathcal A^*_\eta\ \text{or its }\tau\text{-selection, or }\mathrm{abstain}\ ;\ \ \mathrm{Ledger}\ \big)$$

The action output **represents the criterion-optimal act under the declared information — never the true value, the true order, or the identified member.** The Ledger rows, mandatory and separately typed (DE-U1–U6 prove none is recoverable from the others):
1. identification row, verbatim from $I_\theta$: $\widehat J$ (or envelope summary), $R_{\mathrm{set}}$-floor, flags — **invariant in all declarations** (H1);
2. conditional rows with axiom tags: class-risk of the selected action, regret over $\widehat{\mathcal Q}$, ranking tier (1: identified / 2: decision-robust / 3: tie-broken) with the $p$-interval;
3. tolerance row: $\eta$ and the tightness diagnostics (DR-M2);
4. echo row: every declaration consumed, including $\tau$ — the answer carries its own audit trail (H6).

## FAILURE — when the model must abstain (or flag)

1. **Empty section** (support values unrealizable under the declared class): misspecification flag fires; no action (frozen unrealizability semantics).
2. **Off-coverage query with no declared class**: identification row reports $R_{\mathrm{set}}=\infty$; if the context demands a finite unconditional guarantee → abstain (F18 + DR-F2).
3. **Ranking at Tier 3 with no declared $\tau$**: abstain (or emit the action *set*) — a strict ordering here would be a hidden measure (DR-S3, DE-R3).
4. **Empty transported class** ($\widehat{\mathcal Q}=\emptyset$ under the declared transport): shift detected — the declarations are jointly inconsistent with history; surface, do not silently widen (DR-L; analogue of the misspecification flag one level up).
5. **Validity self-check failure** (envelope non-containment detected on held-out identified structure; coverage audit failure of $\widehat{\mathcal Q}$): the certificates are void; the only honest outputs are flags and the frozen minimax fallback (DR-M1 monotone degradation).
6. **Tie-boundary proximity at demanded single-valued continuous output**: the exact selector is discontinuous (DR-S4); emit $\mathcal A^*_\eta$ or flag localized error — a continuous single-valued output here is provably wrong somewhere.

---

**Compilation criterion.** A system is a valid realization iff: it factors as $D_\psi\circ(I_\theta\times M_\phi)$ with the typed channels (DR-M impossibility guard); each factor satisfies its validity predicate ($V_I,V_M,V_D$); and its emissions instantiate this contract's Output/Failure clauses. Its claims are then falsifiable by: the frozen protocols P1–P10, NP-1…6, T1–T5; plus the Phase-8 checks — floor consistency (no emitted unconditional guarantee under $R_{\mathrm{set}}$), outer-envelope audits (order-set containment on constructed families, e.g. the DR-J3 reversal witness), no-hidden-measure audits (symmetric-input/symmetric-output tests at ties, DE-S4/DR-S3), and the three-width accounting of DR-L4 (fiber-count, coverage, transport ablations must move exactly their own width term).
