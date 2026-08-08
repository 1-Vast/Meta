# Loss-Typed Robust Information Floor (Part I)

> **Status:** Phase-8 (decision-operator realizability), 2026-08-03. Phases 0–7 are **frozen and cited, not modified**. New results carry **DR-F** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Purpose: replace the scalar uncertainty quantities (radius $\tfrac12\omega(2\varepsilon)$; per-query Chebyshev radius) by a single decision-typed object that any future approximating system must respect as a floor.

---

## 1. The object

**Definition DR-F0.** For a nonempty identified joint pushforward $J\subseteq\mathbb R^m$ (canonically $J=J_Q(O)$), an action set $\mathcal A$, and a loss $L:\mathcal A\times\mathbb R^m\to[0,\infty]$:

$$R_{\mathrm{set}}(J,\mathcal A,L)\;=\;\inf_{a\in\mathcal A}\ \sup_{v\in J}\ L(a,v)\ \in[0,\infty],$$

and its randomized companion
$$R_{\mathrm{rand}}(J,\mathcal A,L)\;=\;\inf_{\xi\in\Delta(\mathcal A)}\ \sup_{v\in J}\ \int L(a,v)\,d\xi(a)\;\le\;R_{\mathrm{set}}.$$

$R_{\mathrm{set}}$ is the minimax value of the *identified-set game*: the best worst-case loss achievable by any action against an adversary confined to the admissible values. It is loss-typed by construction — no scalar geometry is presupposed.

**Proposition DR-F1 (the mandated specializations). [proved]**
(i) *Scalar absolute error* ($m=1$, $\mathcal A=\mathbb R$, $L=|a-v|$): $R_{\mathrm{set}}=\operatorname{rad}(J)=\tfrac12\operatorname{diam}\overline J$ for bounded $J$ (Chebyshev radius; $+\infty$ if unbounded) — the frozen conditional radius, and in the worst case over data exactly the frozen minimax identity $\tfrac12\omega_{x,D}(2\varepsilon)$ (Theorem 1). The floor is a strict generalization, not a replacement.
(ii) *Multi-query sup-loss:* $R_{\mathrm{set}}=\max_j\operatorname{rad}(\operatorname{proj}_jJ)$ — coordinatewise radii (frozen F1 Rem. 1.1 / DE-J2).
(iii) *Ranking losses:* $\mathcal A$ = orderings (with or without abstention), $L$ = $0$–$1$ or discordance counts. Both signs admissible (DE-R1(iii)) gives $R_{\mathrm{set}}=1$, $R_{\mathrm{rand}}=\tfrac12$ for a pair; the listwise case is computed from the admissible-order object (DR-J3, `joint_query_decision_objects.md`) — e.g. the two-member reversal witness there has $R_{\mathrm{set}}=2$ discordant pairs.
(iv) *Structured actions:* $\mathcal A$ arbitrary (permutations, subsets, monotone maps) — the definition never used linear or metric structure on $\mathcal A$. Existence of a minimizer needs only $\mathcal A$ compact and $L(\cdot,v)$ l.s.c. (DE-T4(i) machinery). $\square$

---

## 2. The floor theorem

**Theorem DR-F2 (no decision method beats the floor without additional information about the current member). [proved]**
Let $\Phi$ be **any** decision rule — deterministic or randomized, and depending on the observations $O$, the entire archive/history $H$, and any declared decision information $\Delta$. Then for every realization of $(O,H,\Delta)$:
$$\sup_{f\in I(O)}\ \mathbb E_{a\sim\Phi(O,H,\Delta)}\ L\big(a,\;e_Q(f)\big)\ \ \ge\ \ R_{\mathrm{rand}}\big(J_Q(O),\mathcal A,L\big),$$
and $\ge R_{\mathrm{set}}$ when $\Phi$ is deterministic.
*Proof.* Condition on $(O,H,\Delta)$: $\Phi$ is then a fixed (mixed) action $\xi$. Every $v\in J_Q(O)$ is realized by some $f\in I(O)$, and every such $f$ is consistent with $O$ — hence with everything the rule saw *about the current member*. So $\sup_{f\in I(O)}\mathbb E L(a,e_Q(f))=\sup_{v\in J_Q(O)}\int L(a,v)\,d\xi\ge R_{\mathrm{rand}}$ by definition of the infimum. History and declarations entered only through the choice of $\xi$, over which the infimum was already taken. $\square$

**Corollary DR-F3 (the exact meaning of "without additional information"). [proved]**
The floor moves only through its first argument. The two legitimate routes are:
(i) **new observational evidence about the current member** (larger $D$, smaller $\varepsilon$, auxiliary fiber label per CI-A) — shrinks $J_Q(O)$: identification, frozen;
(ii) **declared member-level structure** (stronger closure class) — shrinks $\mathcal F$, hence $J_Q(O)$: an axiom, frozen rules.
Population/frequency information moves neither (DE-H2/H3); it can lower only **conditional, tagged** risk statements. A system emitting an *unconditional* guarantee $<R_{\mathrm{rand}}(J_Q(O),\mathcal A,L)$ at any configuration contradicts DR-F2 — the loss-typed form of DE-L5(i)/DE-O4. $\square$

---

## 3. Approximation calculus — why the floor is realizable

A trainable system will not hold $J$ exactly; it holds an approximation $\widehat J$. The floor survives exactly under **outer semantics** (the Phase-5 envelope discipline, now loss-typed):

**Theorem DR-F4 (conservative validity + Lipschitz tightness). [proved]**
(i) *Monotonicity:* $J\subseteq J'\Rightarrow R_{\mathrm{set}}(J,\cdot,\cdot)\le R_{\mathrm{set}}(J',\cdot,\cdot)$ (sup over a larger set, then inf). Hence any **outer** approximation $\widehat J\supseteq J_Q(O)$ yields $R_{\mathrm{set}}(\widehat J)\ge R_{\mathrm{set}}(J_Q(O))$: the reported floor is valid (merely conservative), and every guarantee derived from it is honest.
(ii) *Tightness:* if moreover $d_H(\widehat J,J)\le h$ (Hausdorff) and $L$ is $\ell_v$-Lipschitz in $v$ uniformly in $a$, then
$$R_{\mathrm{set}}(J)\ \le\ R_{\mathrm{set}}(\widehat J)\ \le\ R_{\mathrm{set}}(J)+\ell_v\,h,$$
(each $\hat v\in\widehat J$ lies within $h$ of some $v\in J$, so per-action sups differ by $\le\ell_v h$). Approximation error degrades the floor linearly and one-sidedly.
(iii) *Falsity of inner semantics:* an under-approximation $\widehat J\subsetneq J$ can report a floor strictly below $R_{\mathrm{set}}(J)$ — a **false certificate**, the loss-typed recurrence of the frozen "no closure assumption, no valid radius". Inner approximations are forbidden for certificate purposes at every level of this program. $\square$

**Remark (discontinuous losses).** For $0$–$1$/discordance losses, (ii)'s Lipschitz hypothesis fails; the valid statement is (i) alone plus: the floor computed from the outer admissible-order object $\widehat\Sigma\supseteq\Sigma(J)$ (DR-J3) is conservative. Tightness then depends on order-boundary geometry, not on $h$ — flagged as the realizability price of discrete losses (see DR-S4 for the matching selection-side discontinuity).

---

## 4. Summary

$$\boxed{\ R_{\mathrm{set}}(J,\mathcal A,L)=\inf_{a}\sup_{v\in J}L(a,v)\ \text{ is the unique loss-typed floor: it specializes to every frozen radius, applies to ranking and}\atop\text{structured actions, is unbeatable by any rule using only current-member information (DR-F2), and survives approximation exactly under outer semantics (DR-F4).}\ }$$
