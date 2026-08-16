# Honest Selection Operator — REPAIRED (Part II)

> **Status:** Phase-8.1, 2026-08-03. Supersedes `../08_decision_operator_realizability/honest_selection_operator.md`. The audit **passed** this file except for one over-broad claim; the repair narrows **DR-S4 → DR-S4-R** and re-anchors the risk objects to the repaired certificate types. DR-S1, DR-S2, DR-S3, DR-S5 are carried verbatim by citation (statements restated only where the floor re-typing touches them).

---

## 1. Carried results (audit: pass)

- **DR-S1 [carried].** $\mathcal A^*(J,\Delta,L)=\arg\min_a\rho(a)$ nonempty closed under compactness + l.s.c.; the maximal honest answer.
- **DR-S2 [carried].** Strict quasiconvexity $\Rightarrow$ singleton; squared-loss prediction under any $\Gamma$-minimax class never needs a tie-break.
- **DR-S3 [carried].** Equivariance obstruction: a $G$-equivariant single-valued selector exists iff $\mathcal A^*$ has a $G$-fixed point; free symmetry on a non-singleton $\mathcal A^*$ (sign-swapped rankings) admits none. The dichotomy is exhaustive and forced: **Option A** return $\mathcal A^*$; **Option B** return $\tau(\mathcal A^*)$ for an explicitly declared, data-independent, declared-structure-only tie-break $\tau$. Hidden measures (including any undeclared reference measure $\mu_0$) are forbidden — they are neither required (Option A exists) nor permitted (they are undeclared preferences).
- **DR-S5 [carried].** $\eta$-argmin honesty: $\sup_a|\hat\rho-\rho|\le\eta/2$ makes every $\hat\rho$-minimizer $\eta$-optimal for $\rho$; the valid claim is "$\eta$-optimal", never "optimal"; $\eta$ joins the ledger.

**Certificate re-anchoring (T1 propagation).** Where $\rho$ is the worst-case criterion computed on the outer envelope, $\rho(a)=\sup_{v\in\widehat J}L(a,v)$: minimizing it is exactly the $G_{\mathrm{cert}}$ policy of DR-F4-R(b), and the emitted per-action number $\rho(\hat a)$ is a **guarantee**, valid because $\widehat J$ is outer. No selection-layer statement may re-type it as a floor. Floor statements attached to a selection cite inner witnesses only (DR-F5-R).

---

## 2. DR-S4-R: discontinuity, correctly scoped (replaces DR-S4)

**Theorem DR-S4-R. [proved]**
Let $\{P_t\}_{t\in[0,1]}$ be a continuous path of decision problems (continuous variation of $J$ and/or $\Delta$ in the sense that $\rho_t(a)$ is jointly continuous). Suppose there are disjoint closed sets $\mathcal A_0,\mathcal A_1\subseteq\mathcal A$ at positive distance (separated components — e.g. the two discrete orderings) such that $\mathcal A^*(P_0)\subseteq\mathcal A_0$ and $\mathcal A^*(P_1)\subseteq\mathcal A_1$. Then **every** single-valued selector along the path has a jump of size at least $\operatorname{dist}(\mathcal A_0,\mathcal A_1)$ at some $t$ — a **branch-switching tie** is crossed, and continuity is impossible there.
**Scope (the audit's correction, adopted):** no jump is forced when the argmin merely becomes non-unique within a *connected* region — e.g. an absolute-loss median interval growing and shrinking admits continuous selections (midpoint of the argmin interval). The realizability warning applies exactly to separated-branch decisions: rankings, discrete structured actions, sign choices — not to convex prediction.
*Proof.* Suppose a selector $s$ with $s(t)\in\mathcal A^*(P_t)$ is continuous on $[0,1]$, in the branch-switching case $\bigcup_t\mathcal A^*(P_t)\subseteq\mathcal A_0\cup\mathcal A_1$. Then $T_0=\{t:s(t)\in\mathcal A_0\}$ and $T_1=\{t:s(t)\in\mathcal A_1\}$ are disjoint, nonempty ($0\in T_0$, $1\in T_1$), cover $[0,1]$, and are closed (preimages of closed sets under a continuous map). This contradicts the connectedness of $[0,1]$. Hence $s$ is discontinuous, and since its values lie in $\mathcal A_0\cup\mathcal A_1$, some discontinuity has jump size $\ge\operatorname{dist}(\mathcal A_0,\mathcal A_1)$. $\square$

**Consequence for approximators (carried, narrowed).** Continuous approximators of single-valued selectors incur irreducible localized error near branch-switching ties only; the honest realization there is the set-valued $\mathcal A^*_\eta$ (or declared abstention priced by its declared loss — DR-F1-R(iii)), not a continuous strict output. Convex-prediction selections carry no such obstruction.

---

## 3. Repaired summary

$$\boxed{\begin{array}{c}\text{Valid outputs: }\mathcal A^*\ \text{(or realizably }\mathcal A^*_\eta\text{) — or }\tau(\mathcal A^*)\ \text{with }\tau\ \text{explicitly declared. No hidden measures.}\\ \text{Single-valued canonical selection exists iff }\mathcal A^*\text{ is a singleton (certified, e.g. strict quasiconvexity).}\\ \text{Forced discontinuity is exactly the separated branch-switching case (DR-S4-R); emitted worst-case numbers are outer-typed guarantees, never floors.}\end{array}}$$
