# Joint-Query Decision Theory (Part V)

> **Status:** Phase-7, 2026-08-03. Frozen corpus cited, not modified. New results carry **DE-J** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Object of study: the finite query set $Q=\{x_1,\dots,x_m\}$ and the joint identified pushforward
> $$J_Q(O)=\{(f(x_1),\dots,f(x_m)):f\in I(O)\}\subseteq\mathbb R^m,$$
> with marginal sections $J_{x_j}(O)=\operatorname{proj}_j J_Q(O)$ (the frozen per-query intervals/sets). Question: for which decisions do the marginals suffice, and when is the joint object (or a joint weighting) irreducibly required?

---

## 1. Two sufficiency anchors

**Proposition DE-J1 (pushforward sufficiency — restated). [proved]** Every decision whose loss depends on $f$ only through its $Q$-values is a function of $J_Q(O)$ (plus a weighting on it, if one is declared). $J_Q(O)$ is the maximal object this phase ever needs; the question is when it can be *compressed to its marginals*. (DE-S1.)

**Proposition DE-J2 (sup-loss anchor — frozen). [frozen citation]** For minimax under sup-loss $\max_j|a_j-v_j|$, the per-query midpoints are exactly optimal and the marginal sections suffice (F1 Rem. 1.1: per-query Chebyshev centers are jointly optimal under sup-loss). This is the unique multi-query case already settled by the frozen theory, and the mandate's observation — scalar marginal intervals suffice for scalar minimax absolute error — is its $m=1$ instance.

---

## 2. The five mandated decisions

Throughout, "weighting" means a law $\mu$ on $J_Q(O)$ (or an ambiguity set of such); "minimax" means the criterion consuming no weighting.

**DE-J3 (point prediction, squared loss). [proved]**
$L(a,v)=\sum_j(a_j-v_j)^2$ is additively separable.
(i) *Bayes under any weighting $\mu$:* $a_j^*=\int v_j\,d\mu$ — coordinatewise posterior means; only the **marginal laws** $\mu_j$ enter (Fubini). Marginals of the *weighting* suffice; note the weighting's marginals are constrained by, but not determined by, the set marginals.
(ii) *Minimax:* $\inf_a\sup_{v\in J_Q}\|a-v\|_2^2$ is the squared Euclidean Chebyshev radius of the **joint** set; marginal sections do not determine it in general, though they do bound it (the joint set sits in the product box).

**DE-J4 (point prediction, absolute loss). [proved]**
$L(a,v)=\sum_j|a_j-v_j|$, separable.
(i) *Bayes:* $a_j^*=$ any median of $\mu_j$ — marginal laws suffice.
(ii) *Minimax:* the joint set is load-bearing. Witness: $Q=\{x_1,x_2\}$, $J^{\mathrm{diag}}=\{(t,t):t\in[0,1]\}$ vs $J^{\mathrm{sq}}=[0,1]^2$ — **identical marginals** $[0,1]\times[0,1]$ (both realizable: $J^{\mathrm{diag}}$ by a one-parameter family with $f_t\equiv t$ near both queries; $J^{\mathrm{sq}}$ by a C1-type family). Over $J^{\mathrm{sq}}$: $\sup$ separates, unique optimum $a=(\tfrac12,\tfrac12)$. Over $J^{\mathrm{diag}}$: $\sup_t(|a_1-t|+|a_2-t|)$ is attained at $t\in\{0,1\}$ (convexity in $t$), giving $\max(a_1+a_2,\,2-a_1-a_2)$, minimized on the whole segment $\{a_1+a_2=1\}\cap[0,1]^2$. Same marginals, different optimal-action sets: the marginal compression already misreports *which actions are optimal*.

**DE-J5 (asymmetric / quantile loss). [proved]**
Pinball loss per coordinate, $\rho_\tau(u)=\tau u^+ +(1-\tau)u^-$.
(i) *Bayes:* $a_j^*=$ the $\tau$-quantile of $\mu_j$ — requires a **marginal law**; the identified interval carries no quantile, so no set-only rule reproduces Bayes-quantile behavior (this is DE-S4's impossibility in loss-specific form).
(ii) *Minimax on a marginal interval $[l,u]$ — closed form:* $\sup_{v}\rho_\tau(v-a)=\max\{\tau(u-a),\,(1-\tau)(a-l)\}$, minimized where the two branches meet:
$$\boxed{\ a^*=(1-\tau)\,l+\tau\,u\ }\qquad\text{value }\ \tau(1-\tau)(u-l).$$
So minimax *does* have a canonical asymmetric answer — but the asymmetry parameter $\tau$ is declared loss structure, consistent with DE-S4 (the loss itself is part of the declared context, not of the identified set).

**DE-J6 (pairwise ranking). [proved]**
The decision-relevant functional is $v_a-v_b$, whose identified set is $\Delta_{ab}(O)=\{f(x_a)-f(x_b):f\in I(O)\}$ — a **linear image of the joint set, invisible to the marginals**. Witness triple with identical marginals $[0,1]^2$:
- $J^{\mathrm{diag}}=\{(t,t)\}$: $\Delta_{ab}=\{0\}$ — the difference is *identified* (tie), minimax error $0$;
- $J^{\mathrm{anti}}=\{(t,1-t)\}$: $\Delta_{ab}=[-1,1]$ — sign fully ambiguous;
- $J^{\mathrm{sq}}=[0,1]^2$: $\Delta_{ab}=[-1,1]$ — ambiguous, but with different joint geometry than $J^{\mathrm{anti}}$.
Minimax absolute error for estimating $v_a-v_b$: $0$ vs $1$ vs $1$ — under identical marginals. This is the decision-layer face of frozen Thm 7.1/C9 (differences identifiable when values are not, and conversely). Ranking theory proper: `ranking_decision_theory.md`.

**DE-J7 (general joint loss). [proved]**
For arbitrary $L:\mathcal A\times\mathbb R^m\to\mathbb R$, the pair (joint set $J_Q(O)$, joint weighting on it) is **minimally sufficient**: sufficiency is DE-J1; minimality holds because (a) any object omitting joint-set information fails on DE-J4/J6's equal-marginal witnesses under some loss, and (b) any object omitting the weighting fails by DE-S3. No compression below the pair is uniformly valid over losses. $\square$

---

## 3. When do marginal sections suffice? The boundary theorem

**Theorem DE-J8. [proved, with scoped converse]**
Marginal information suffices in exactly the following regimes, and fails outside them in general:
(i) **[proved]** *Sup-loss minimax:* marginal sections suffice (DE-J2, frozen).
(ii) **[proved]** *Separable loss + product-form weighting:* Bayes risk splits coordinatewise; marginal laws suffice (Fubini; DE-J3(i)/J4(i) are instances — note any joint $\mu$ enters only via its marginals here, so "product-form" can be weakened to "marginals declared").
(iii) **[proved]** *Failure of (ii) without separability:* for $L=\ell(v_a-v_b,\cdot)$-type couplings (DE-J6), even a full product weighting on the marginal box misstates the risk, because the joint set/law of the difference is not determined by marginals — witnesses above.
(iv) **[proved]** *Failure of (i) beyond sup-loss:* separable minimax already needs the joint set (DE-J4(ii)); the frozen radius/diameter caveat (Thm 1 scope note: "for other joint losses the radius can exceed half the diameter") is the same phenomenon seen from the identification side.
Summary: **marginal sections are sufficient iff the criterion–loss pair factors through coordinates — separable Bayes with declared marginal laws, or sup-loss minimax. Every genuinely joint decision (ranking above all) consumes the joint pushforward.** $\square$

*Consequence for the interface.* The frozen operator's per-query output (center + radius per $x$) is decision-complete only for the regimes (i)–(ii). A decision extension serving ranking or non-separable losses must carry $J_Q(O)$ — or at least the identified sets of the decision-relevant linear functionals (e.g. $\Delta_{ab}$, which the frozen linear theory computes exactly via Thm 7.1's row-space test). This is an *interface widening*, not a change to $I$: $J_Q(O)$ and every $\Delta_{ab}$ are images of the same frozen admissible set.
