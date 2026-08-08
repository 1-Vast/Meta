# Meta-Learning Interface (Part V)

> **Status:** Phase-8, 2026-08-03. Phases 0–7 frozen and cited. New results carry **DR-M** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Mandate: define **exactly three** mathematical objects — their domains, codomains, and validity predicates. Subscripts $\theta,\phi,\psi$ mark the approximable degrees of freedom; they range over *unspecified* families. No neural implementation, parameterization, or training loss is chosen or implied.

---

## 0. Spaces

- $\mathcal X$: query domain (structureless, frozen). $\mathcal O=(\mathcal X\times\mathbb R)^{\le k}\times[0,\infty)$: current-task observation records $(S_t,\varepsilon)$.
- $\mathcal W$: archive-evidence space — finite relational records of historical traces $\{(D_i, y_i, c_i)\}_{i\le n}$ (designs, values, optional auxiliary labels).
- $\mathcal K_m$: nonempty closed subsets of $\mathbb R^m$; $\mathrm{Flags}$: the frozen partiality flags (unrealizable data / off-coverage / unbounded section) $\cup$ Phase-8 flags (tie-boundary, empty transport class).
- $\mathfrak Q_g$: convex weak-*-compact classes of laws on the decision-relevant pushforward space of $g$ (Part III), with confidence tags.
- $\mathrm{Ctx}$: decision contexts $(\mathcal A, L, \text{criterion}, \tau, \eta,\text{axiom declarations})$ — all declared, all echoed.
- $\mathrm{Ledger}$: the Phase-7 four-row object (`decision_uncertainty_ledger.md`) extended by the Phase-8 rows ($R_{\mathrm{set}}$-floor value; $\eta$-optimality tolerance).

---

## 1. Identification operator $I_\theta$

$$I_\theta:\ \mathcal W\times\mathcal O\times\mathcal X^m\ \longrightarrow\ \mathcal K_m\times\mathrm{Flags},\qquad (w,\,o,\,Q)\ \mapsto\ (\widehat J,\ \mathrm{flags}).$$

**Validity predicate $V_I$ [the correctness definition, not an aspiration]:** whenever the declared closure class is correct, $\widehat J\ \supseteq\ J_Q(O)$ (outer envelope; equivalently, for order decisions, $\widehat\Sigma\supseteq\Sigma(J_Q(O))$ after projection), with flags surfacing the three frozen partiality sources. **Tightness objective:** $d_H(\widehat J, J_Q(O))$ (quality, never validity). All frozen operator axioms (A1–A10 of Phase 3; Phase-5 realizability constraints) apply to $I_\theta$ unchanged; $\theta$ is the slot the archive's *feasibility channel* (DE-H1(i)) is allowed to shape.

## 2. Population adaptation operator $M_\phi$

$$M_\phi:\ \mathcal W\times\{\text{axiom declarations}\}\ \longrightarrow\ \mathfrak Q_g\times(0,1],\qquad (w,\,\text{tags})\ \mapsto\ (\widehat{\mathcal Q},\ 1-\delta).$$

**Validity predicate $V_M$:** under the declared rungs (EXCH or C-EXCH; TRANS radius; per-task gate), $\Pr\big(\text{true conditioned population law}\in\widehat{\mathcal Q}\big)\ge1-\delta$ — outer validity at the population level (DR-L2(iv)), with the DR-L3 interval class as the canonical construction. $\phi$ is the slot the archive's *frequency channel* (DE-H1(ii)) is allowed to shape. **Hard constraint:** $M_\phi$'s output never feeds back into $I_\theta$'s set (DE-H2/H3 — the typing rule that makes frequency-driven shrinkage a type error, not a runtime error).

## 3. Decision operator $D_\psi$

$$D_\psi:\ \big(\mathcal K_m\times\mathrm{Flags}\big)\times\big(\mathfrak Q_g\times(0,1]\big)\times\mathrm{Ctx}\ \longrightarrow\ \big(2^{\mathcal A}\setminus\{\emptyset\}\ \cup\ \{\mathrm{abstain}\}\big)\times\mathrm{Ledger}.$$

**Validity predicate $V_D$:** the Phase-7 honesty axioms H1–H6, plus the Phase-8 additions:
- **(floor consistency)** every emitted unconditional guarantee is $\ge R_{\mathrm{set}}(\widehat J,\mathcal A,L)$ — valid because $\widehat J$ is outer (DR-F4(i));
- **(honest selection)** the action output is $\mathcal A^*_\eta$ (or its declared-$\tau$ selection) with the tolerance $\eta$ in the Ledger (DR-S5); single-valued output without $\tau$ only under the DR-S2 uniqueness certificate;
- **(abstention)** the codomain's $\mathrm{abstain}$ element is mandatory in the cases enumerated in `THEORY_TO_MODEL_INTERFACE.md` §Failure.

$\psi$ is the slot where the *criterion computation* (conditioning, $\Gamma$-minimax evaluation, tie-break application) may be approximated.

---

## 4. Composition

$$\mathbb D_{\mathrm{real}}\;=\;D_\psi\circ\big(I_\theta\times M_\phi\big):\quad (\text{archive},\ \text{current observations},\ Q,\ \text{context})\ \longmapsto\ (\text{action set/abstain},\ \text{Ledger}).$$

**Theorem DR-M1 (honest composition). [proved]**
If $V_I$, $V_M$, $V_D$ hold, then every statement emitted by $\mathbb D_{\mathrm{real}}$ is valid: unconditional rows by DR-F4(i) (outer $\widehat J$ $\Rightarrow$ conservative floor and radius), conditional rows by $V_M$ (outer class $\Rightarrow$ $\Gamma$-minimax risk is a valid upper bound with confidence $1-\delta$) and DR-S5 ($\eta$-optimality claims), and the separation rows by H1/DE-O3 (selection is epistemically inert). Moreover the composition degrades **monotonically**: any single validity predicate weakened to its conservative limit (vacuous $\widehat J=$ everything admissible; vacuous $\widehat{\mathcal Q}=$ all laws) collapses that factor to the frozen fallback (identification: report only flags; decision: DE-T4(iii) minimax) without invalidating any other factor. $\square$

**Theorem DR-M2 (tightness error calculus). [proved]**
Approximation losses compose additively through declared moduli: floor inflation $\le\ell_v\,d_H(\widehat J,J)$ (DR-F4(ii), Lipschitz losses); conditional-risk inflation $\le$ the DR-L4 three-width sum transferred with the loss's BV/Lipschitz constant (DE-L3(iii)); selection tolerance $\eta\ \ge\ 2\times$ total risk-approximation error (DR-S5(i)). The tightness budget is therefore *auditable end-to-end*, and every term sits in the DE-U7 decomposition with its Phase-8 name. Discrete-loss caveat: at order boundaries the calculus is one-sided only (DR-F4 remark, DR-S4). $\square$

**Impossibility guard (why exactly three objects). [proved]**
Two mergers are forbidden by frozen theorems, so the factorization is minimal, not stylistic: merging $M_\phi$ into $I_\theta$ (one operator producing a "posterior-narrowed set") violates DE-H2/H3 — frequency would shrink the set; merging $M_\phi$ into $D_\psi$ without the explicit $\mathfrak Q_g$ interface hides the measure — violating DR-S3's no-hidden-measure dichotomy and H6 auditability. A unary form was already impossible at the identification level (OP-10). Three is the minimum that types the information flow honestly. $\square$
