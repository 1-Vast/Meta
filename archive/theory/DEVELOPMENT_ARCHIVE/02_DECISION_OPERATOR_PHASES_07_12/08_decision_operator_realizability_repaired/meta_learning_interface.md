# Meta-Learning Interface — REPAIRED (Part V)

> **Status:** Phase-8.1, 2026-08-03. Supersedes `../08_decision_operator_realizability/meta_learning_interface.md`. Repairs audit targets **T1** ($V_D$/DR-M1 floor typing), **T4** ($M_\phi$ domain typing), **T7** (loss-typed fallback). The three-operator separation, the spaces, and the minimality guard are carried; the operator signatures and composition theorem are restated as **DR-M-R / DR-M1-R / DR-M2-R**.

---

## 0. Spaces (carried, one addition)

As before ($\mathcal X,\mathcal O,\mathcal W,\mathcal K_m,\mathrm{Flags},\mathrm{Ctx},\mathrm{Ledger}$), plus:
- $\mathcal G$: the space of **decision specifications** $\gamma=(g,\,Q,\,\mathcal A,\,L,\,\text{criterion},\,\delta,\,\text{axiom tags})$ — the index the population operator was missing;
- $\mathfrak Q_\gamma$: convex weak-*-compact classes of laws on the $g$-pushforward space of $\gamma$;
- the Ledger's worst-case rows are re-typed: **(floor-bound row)** $R_{\mathrm{set}}(\widetilde J,\mathcal A,L)$ from certified witnesses (present when witnesses are computed, else marked absent); **(guarantee row)** $G_{\mathrm{cert}}(\widehat J,\mathcal A,L)$; the pair brackets the uncomputable true floor (DR-F4-R).

---

## 1. Identification operator $I_\theta$ (carried)

$$I_\theta:\ \mathcal W\times\mathcal O\times\mathcal X^m\ \longrightarrow\ \mathcal K_m\times\mathrm{Flags},\qquad V_I:\ \widehat J\supseteq J_Q(O)\ \text{under the declared closure class (outer envelope; order projection }\widehat\Sigma\supseteq\Sigma).$$
Optional certified extension (DR-F5-R): a finite witness list $\widetilde J\subseteq J_Q(O)$ of verified admissible members' query vectors — an inner *floor device*, never a feasibility report; emitting $\widetilde J$ in place of $\widehat J$ is a type error.

## 2. Population adaptation operator $M_\phi$ — retyped (T4)

**Indexed form (primary):**
$$M_\phi:\ \mathcal W\times\mathcal G\ \longrightarrow\ \mathfrak Q_\gamma\times(0,1],\qquad (w,\gamma)\ \mapsto\ (\widehat{\mathcal Q}_\gamma,\ 1-\delta).$$
The decision specification $\gamma$ — including $g$, $Q$, and the context — is an **explicit input**; the codomain is the class space *of that index*. $V_M$: under $\gamma$'s declared rungs ((IID)/(C-IID)/declared concentration; transport; per-task gate), $\Pr(\text{true conditioned law}\in\widehat{\mathcal Q}_\gamma)\ge1-\delta$, with **simultaneous** coverage over $\gamma$'s finite target family per DR-L3-R.

**Universally indexed form (equivalent alternative):**
$$M_\phi:\ \mathcal W\times\{\text{axiom tags}\}\ \longrightarrow\ \prod_{\gamma\in\mathcal G_0}\big(\mathfrak Q_\gamma\times(0,1]\big)\quad\text{for a declared index set }\mathcal G_0,$$
subject to **projective consistency**: if $g'=h\circ g$ (a coarser pushforward), then $h_*\widehat{\mathcal Q}_{\gamma}\subseteq\widehat{\mathcal Q}_{\gamma'}$ — the family may not contradict itself across specifications — and to a **joint** confidence accounting over $\mathcal G_0$ (one union-bound/uniform declaration covering the whole family, DR-L2-R(v)). Either form is admissible; an output class whose index appears in neither the input nor a declared $\mathcal G_0$ is untyped and forbidden.

Hard constraint carried: no $M_\phi\to I_\theta$ edge (DE-H2/H3).

## 3. Decision operator $D_\psi$ — retyped $V_D$ (T1, T7)

$$D_\psi:\ \big(\mathcal K_m\times\mathrm{Flags}\ [\times\ \text{witness lists}]\big)\times\big(\mathfrak Q_\gamma\times(0,1]\big)\times\mathrm{Ctx}\ \longrightarrow\ \big(2^{\mathcal A}\setminus\{\emptyset\}\ \cup\ \{\text{declared abstain}\in\mathcal A\}\big)\times\mathrm{Ledger}.$$

**$V_D$ (repaired):**
- **(guarantee validity)** every emitted unconditional worst-case number is computed as $\sup_{v\in\widehat J}L(\cdot,v)$ on the outer envelope — hence a valid upper bound on true worst-case risk (DR-F4-R(b)); the best such number is $G_{\mathrm{cert}}$;
- **(floor honesty)** the word "floor" attaches only to $R_{\mathrm{set}}(J_Q(O))$ symbolically or to certified witness values $R_{\mathrm{set}}(\widetilde J)$ numerically (DR-F4-R(c)); presenting any outer-envelope value as a floor is a type violation;
- **(honest selection)** $\mathcal A^*_\eta$ or declared-$\tau$ selection, with $\eta$ in the Ledger; single-valued without $\tau$ only under the DR-S2 uniqueness certificate; branch-switching caveat DR-S4-R;
- **(loss-typed fallback — T7)** every degraded/failure behavior is evaluated under the declared $(\mathcal A,L)$: off-coverage or vacuous-class configurations report $G_{\mathrm{cert}}(\widehat J,\mathcal A,L)$ — infinite only when $L$ is unbounded on the feasible set; finite ($\le1$, $\le c$) for $0$–$1$/ranking/abstention losses. Abstention is a declared action with a declared loss, selected when criterion-optimal or when demanded single-valuedness is untypable — never a universal-$\infty$ reflex;
- Phase-7 axioms H1–H6 carried.

---

## 4. Composition — repaired theorems

**Theorem DR-M1-R (honest composition). [proved]**
If $V_I$, $V_M$ (repaired), $V_D$ (repaired) hold, every emission of $\mathbb D_{\mathrm{real}}=D_\psi\circ(I_\theta\times M_\phi)$ is valid under its stated type: guarantee rows by outerness of $\widehat J$; floor-bound rows by witness membership + monotonicity; conditional rows by the simultaneous coverage of $\widehat{\mathcal Q}_\gamma$ under the declared rungs and by DR-S5's $\eta$-claims; separation rows by H1/DE-O3. Degradation is monotone and loss-typed: a vacuous $\widehat J$ or $\widehat{\mathcal Q}_\gamma$ collapses its factor to the conservative endpoint (report flags; DE-T4 minimax under the declared $(\mathcal A,L)$) without falsifying any emission. *The Phase-8 version of this theorem invoked the retracted floor reading of outer envelopes; no step above does.* $\square$

**Theorem DR-M2-R (tightness calculus, re-typed). [proved]**
The auditable end-to-end deficit is the **bracket width** $G_{\mathrm{cert}}(\widehat J)-R_{\mathrm{set}}(\widetilde J)$ plus the conditional-row widths (DR-L4-R's three terms transferred with the loss's BV/Lipschitz constant) plus $\eta$. For Lipschitz losses, $G_{\mathrm{cert}}\le R_{\mathrm{set}}(J)+\ell_v\,d_H(\widehat J,J)$; the bracket closes at rate of the envelope tightness plus witness richness. Discrete-loss boundary caveat carried. $\square$

**Minimality guard (carried).** Merging $M_\phi$ into $I_\theta$ violates DE-H2/H3; merging $M_\phi$ into $D_\psi$ without the explicit class interface violates DR-S3/H6; unary forms fail at the identification level (OP-10). Three typed objects remain the minimum.
