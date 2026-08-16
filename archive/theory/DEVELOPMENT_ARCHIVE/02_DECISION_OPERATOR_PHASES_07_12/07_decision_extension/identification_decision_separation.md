# Identification–Decision Separation (Part I)

> **Status:** Phase-7 (decision extension), 2026-08-03. Phases 0–6 are **frozen and cited, not modified**. The identified object is fixed: given observations $O$, the nonempty admissible set $I(O)\subseteq\mathcal F$ (Theorem 1 / CP-2; under Phase-6 relaxations, the union-fiber admissible set of CI-D — the results below hold for any nonempty $I(O)$, so the Phase-6 operator plugs in verbatim). New results carry **DE-S** numbers, tagged **[proved] / [conditional on stated assumptions] / [impossible] / [open]**. Nothing in this file changes $I$; everything in it is about what $I$ does and does not determine downstream.

---

## 1. The two maps

**Identification (frozen).**
$$I:\ O\ \longmapsto\ I(O)\;=\;\{f\in\mathcal F:\ \max_i|\tilde y_i-f(x_i)|\le\varepsilon\}\ \ne\ \emptyset,$$
and for a finite query set $Q=\{x_1,\dots,x_m\}$ its pushforward
$$J_Q(O)\;=\;\{(f(x_1),\dots,f(x_m)):f\in I(O)\}\subseteq\mathbb R^m .$$
$I$ is a function of the observations and the declared family/noise model alone. Its optimality theory (center + radius) is Theorem 1 and is not reopened.

**Decision (new).** A *decision context* is a triple $\kappa=(\mathcal A,L,\mathsf C)$: an action set $\mathcal A$, a loss $L:\mathcal A\times\mathbb R^m\to\mathbb R$ (through the query values; DE-S1 justifies this reduction), and a selection criterion $\mathsf C$. The decision map is
$$D:\ \big(\underbrace{J_Q(O)}_{\text{identified object}},\ \underbrace{\Delta}_{\text{declared decision information}}\big)\ \longmapsto\ a\in\mathcal A,$$
where $\Delta$ contains everything the criterion needs beyond the set (weightings, ambiguity classes, reference points, declared invariances — Part II determines the *minimal* such object). The separation demanded by this phase is the statement, proved below, that the $\Delta$-slot cannot in general be empty.

**Proposition DE-S1 (pushforward sufficiency). [proved]**
If the loss depends on $f$ only through $(f(x))_{x\in Q}$, then for every criterion whose input is the family of loss profiles, the decision problem is a function of $J_Q(O)$ alone; $I(O)$ carries no further decision-relevant information.
*Proof.* The profile of an action, $f\mapsto L(a,(f(x))_{x\in Q})$, factors through the evaluation map $e_Q:f\mapsto(f(x))_{x\in Q}$; hence the *set* of achievable loss values of every action, and any functional of it, is determined by $e_Q(I(O))=J_Q(O)$. $\square$
*Caveat (weighting-level failure):* a weighting on $I(O)$ is **not** determined by its pushforward to $J_Q(O)$ when $e_Q$ is non-injective, but its decision-relevant content is exactly the pushforward weighting on $J_Q(O)$ — so the reduction survives at the weighted level too, with "weighting on $J_Q(O)$" as the reduced object.

---

## 2. What identification *does* determine about decisions

Identification is not decision-mute. It determines exactly the **dominance structure**.

**Definition.** For actions $a,a'$ write $a\preceq_{\mathrm{dom}} a'$ ("$a$ weakly dominates") iff $L(a,v)\le L(a',v)$ for all $v\in J_Q(O)$; strict if additionally $<$ somewhere on $J_Q(O)$.

**Proposition DE-S2 (the identified set determines exactly admissibility). [proved]**
(i) $\preceq_{\mathrm{dom}}$ is a function of $J_Q(O)$ (immediate). (ii) *Scalar complete class in miniature:* let $m=1$, $J:=J_{x}(O)$ bounded, $\mathcal A=\mathbb R$, and $L(a,v)=\varphi(a-v)$ with $\varphi$ strictly decreasing on $(-\infty,0]$, strictly increasing on $[0,\infty)$, $\varphi(0)=0$ (covers absolute and squared loss). Then the undominated actions are exactly the closed interval $[\inf J,\sup J]$; and for squared loss the Bayes actions over all probability weightings supported on $J$ sweep exactly the closed convex hull $\overline{\operatorname{conv}}\,J=[\inf J,\sup J]$.
*Proof.* Any $a>\sup J$ is strictly dominated by $\sup J$ ($|a-v|>|\sup J-v|$ for all $v\in J$); symmetrically below. For $a,a'\in[\inf J,\sup J]$, $a\ne a'$: pick $v\in J$ with $v$ on $a$'s side (exists since $\inf,\sup$ are approached), then $L(a,v)<L(a',v)$ near the relevant end, and symmetrically — neither dominates. Bayes-squared actions are means $\int v\,d\pi(v)$, $\pi\in\Delta(J)$; means of probability measures on a bounded $J\subseteq\mathbb R$ sweep exactly $[\inf J,\sup J]$ (point masses give closure points; convexity gives everything between). $\square$

So: **identification delivers a partial order on actions (dominance) and hence the admissible action set — and nothing finer.** The rest of this file proves the "nothing finer".

---

## 3. Why $I$ alone cannot determine $D$

**Theorem DE-S3 (two-context underdetermination). [proved]**
There is no map $d:\{\text{nonempty subsets of }\mathbb R\}\to\mathbb R$ such that $d(J)$ is an optimal action in every admissible decision context with identified pushforward $J$.
*Proof.* Take $J=\{0,1\}$, realized inside the frozen theory (e.g. $\mathcal F$ of C1-type: two members agreeing on $D$ within $2\varepsilon$ with query values $0$ and $1$; both are in $I(O)$). Context $\kappa_1$: squared loss, Bayes criterion with weighting $(\tfrac13,\tfrac23)$ on $(0,1)$ — unique optimum $a=\tfrac23$. Context $\kappa_2$: same loss, weighting $(\tfrac23,\tfrac13)$ — unique optimum $a=\tfrac13$. Both weightings are supported on $I(O)$ and are consistent with every observation about the *current* member (weightings are not observables of the current member; no datum in $O$ discriminates them). A single value $d(\{0,1\})$ differs from at least one of $\tfrac23,\tfrac13$, hence is suboptimal in at least one admissible context. $\square$

**Theorem DE-S4 (equivariance impossibility: no canonical non-minimax rule). [proved]**
Call a selection rule $\sigma$ mapping nonempty bounded sets $J\subseteq\mathbb R$ to actions *intrinsic* if it depends on $J$ alone, and *affine-equivariant* if $\sigma(\alpha J+\beta)=\alpha\,\sigma(J)+\beta$ for all $\alpha\ne0,\beta\in\mathbb R$ (reflections included). Then every intrinsic affine-equivariant $\sigma$ satisfies $\sigma(J)=c$ whenever $J$ is symmetric about $c$ — i.e. it agrees with the Chebyshev center (the frozen minimax action for every loss $\varphi(|a-v|)$, $\varphi$ strictly increasing) on all symmetric sets.
*Proof.* If $J=2c-J$, apply equivariance with $(\alpha,\beta)=(-1,2c)$: $\sigma(J)=\sigma(2c-J)=2c-\sigma(J)$, so $\sigma(J)=c$. For symmetric $J\subseteq[c-r,c+r]$ with $c\pm r\in\overline J$, $\sup_{v\in J}\varphi(|a-v|)=\varphi(\max(|a-(c-r)|,|a-(c+r)|))$ is uniquely minimized at $a=c$. $\square$
*Consequence.* Any rule that at some symmetric identified set outputs a non-center action — every strictly asymmetric Bayes rule, every Bayes-quantile decision with $\tau\ne\tfrac12$, every asymmetric Hurwicz mixture — must either use information beyond the identified set or break affine equivariance, i.e. import a declared asymmetry of the value line. Either way it consumes a $\Delta$. Non-minimax decisions are never canonical functions of the identified set.
*Scope (honest).* Equivariance pins $\sigma$ only on sets with nontrivial affine symmetry; on an asymmetric set like $\{0,1,10\}$ (whose affine symmetry group is trivial) equivariance alone forces nothing. The theorem is an impossibility of *canonicity*, established where the symmetry has teeth — which suffices, since a canonical rule must be canonical everywhere.

**Proposition DE-S5 (no canonical weighting; the gauge obstruction). [proved, scoped]**
A "uniform prior over the identified set" is not a canonical object:
(i) on a continuum $J$, "uniform" is relative to a declared group — Lebesgue-uniform on an interval is affine-canonical but not invariant under the order-isomorphism group, and on $I(O)\subseteq\mathcal F$ itself no canonical reference measure exists;
(ii) the family's parametrization is unidentifiable (frozen CP §4.5: gauge and even the member set are not determined by finite windows), so "uniform over parameters" and "uniform over values" differ and neither is privileged by any observation;
(iii) *exception, stated honestly:* on a **finite** identified value set, counting measure is invariant under all permutations — a genuinely canonical weighting — but it is canonical only *given* the decision to weight values rather than members or parameters, and (ii) shows that choice is itself unidentified. Declaring the invariance class **is** decision information. $\square$

---

## 4. Same identified set, different optimal actions — worked constructions

All contexts below share the identified pushforward; every listed weighting is supported inside it; all actions differ. (Realizability of the sets inside the frozen theory: intervals arise from the Lipschitz class C2; finite sets from C1-type families.)

**(a) $J=[0,1]$, $\mathcal A=\mathbb R$.**

| Criterion (the $\Delta$ consumed) | Optimal action |
|---|---|
| minimax absolute error (no $\Delta$; frozen Thm 1) | $0.5$ |
| minimax pinball loss, $\tau=0.9$ (declared asymmetry; DE-J6) | $0.9$ |
| Bayes squared, weighting $=$ Beta-like mass toward $0$ (declared law) | $\approx$ its mean, e.g. $0.25$ |
| $\Gamma$-minimax squared over a TV-ball around uniform (declared class) | an interval-valued argmin near $0.5$, radius set by the ball |

**(b) $J=\{0,1,10\}$, $\mathcal A=\mathbb R$, all weightings uniform where used.**

| Criterion | Optimal action |
|---|---|
| minimax (absolute or squared) — Chebyshev center | $5$ |
| Bayes squared, uniform on the three values | $11/3\approx3.67$ |
| Bayes absolute, uniform on the three values (a median) | $1$ |
| Bayes squared, uniform on *parameters* of a family hitting $10$ twice (gauge choice, DE-S5(ii)) | $21/4=5.25$ |

One set; four defensible actions; the difference between any two is carried entirely by $\Delta$, not by $O$. This is the content of DE-S3 made concrete.

---

## 5. Summary

$$\boxed{\ I\ \text{determines the dominance order and the admissible action set — exactly; selecting}\atop\text{one undominated action requires a declared object }\Delta\ \text{that no observation about the current member supplies.}\ }$$

The frozen theory already contains the unique $\Delta$-free selection: the minimax/Chebyshev rule, which DE-S4 shows is the only canonical (intrinsic, equivariant) behavior. Every departure from it is a decision, not an inference. Part II (`minimal_decision_primitive.md`) determines the weakest admissible $\Delta$.
