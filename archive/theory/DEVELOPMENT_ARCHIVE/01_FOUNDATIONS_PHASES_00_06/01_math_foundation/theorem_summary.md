# Theorem Summary — Few-Shot Identifiability of Function-Family Members

> **Status:** frozen mathematical foundation, 2026-08-02. Source: `../00_raw_outputs/identifiability_treatise_raw.md` (adversarially verified; see provenance note there). Section references (§, Thm numbers) point into the raw treatise. This document states results; proofs live in the source.
>
> **Scope discipline:** no neural architecture, no application mapping. Any future model design must trace each of its components back to a numbered item here.

---

## 1. The central mathematical problem

Let $\mathcal F\subseteq\mathbb R^{\mathcal X}$ be a nonempty family of real-valued functions on an arbitrary set $\mathcal X$ (no topology, metric, or representation assumed). Data come in three strata:

- **Archive:** many members $f_\alpha$, each known only on its own finite design $D_\alpha$; covered set $U=\bigcup_\alpha D_\alpha$.
- **Sample:** a new member $f_\beta$ observed as $\tilde y\in\mathbb R^k$ with $\max_i|\tilde y_i-f_\beta(x_i)|\le\varepsilon$ on a design $D=\{x_1,\dots,x_k\}$, $k\le5$.
- **Query:** the value $f_\beta(x)$ is required at $x\notin D$.

**Question:** when, and exactly how well, can any estimator $\Phi:\mathbb R^k\to\mathbb R$ recover $f_\beta(x)$?

The entire theory reduces to two derived objects:

- the **trace set** $T_{D,x}=\{(f|_D,f(x)):f\in\mathcal F\}\subseteq\mathbb R^{k+1}$;
- the **trace modulus** $\omega_{x,D}(t)=\sup\{|f(x)-g(x)| : f,g\in\mathcal F,\ \max_i|f(x_i)-g(x_i)|\le t\}$,
  generalized to any functional $\varphi$ as $\omega_\varphi(t)$.

Every positive and negative result below is a statement about these two objects.

---

## 2. The strongest impossibility result

**Exact minimax identity (Thm 1).** With no assumptions on $\mathcal F$ beyond nonemptiness, for every design, query, and noise level, in $[0,\infty]$:

$$\inf_\Phi\ \sup_{f,\ \tilde y\ \text{consistent}}\ |\Phi(\tilde y)-f(x)| \;=\; \tfrac12\,\omega_{x,D}(2\varepsilon).$$

The bound survives randomized estimators (in expected worst-case error); the achieving rule is the midpoint (Chebyshev center) of the consistent-value section and depends on $\varepsilon$; the prefactor ½ is exact because the target is scalar (for sup-loss over many queries the per-query midpoint remains exactly optimal; for other joint losses the radius can exceed half the diameter).

**Strongest corollary (Cor. 1.4).** For $\mathcal F=\mathbb R^{\mathcal X}$ — or any family containing pairs that agree on $D$ and differ arbitrarily at $x$ — the minimax error is $\varepsilon$ at each design point and $+\infty$ at **every** other point. *The set of points where the data carry any information (finite minimax error) is exactly $D$.* Knowing every other member of the family completely does not help unless the family itself constrains the trace set.

**The impossibility recurs at family level (Thm 9.2).** Even for a $d$-dimensional linear family: at any query $x$ outside the archive-covered set $U$ (and outside $D$), the set of values consistent with *all* archive and sample data is either all of $\mathbb R$, or $\{0\}$ in the single degenerate case where the data force $f_\beta=0$ (an exception that vanishes when $\varepsilon>0$). Archive data about other points never transports information to an uncovered query; only member-level assumptions (continuity toward $\overline U$, analyticity, RKHS membership) can.

---

## 3. Minimal assumptions that restore identifiability

A strict hierarchy; each level is *equivalent* to a grade of inference (the strongest form of minimality — every sufficient assumption implies it):

| Level | Assumption | Exactly equivalent to | Certificate |
|---|---|---|---|
| **A0** | $\omega_{x,D}(0)=0$: the trace set $T_{D,x}$ is the **graph of a function** $\Psi$ over its projection $\mathcal F\vert_D$ | exact-data identifiability of $f_\beta(x)$ at $(D,x)$ | Prop. 1.5 (iff); $\Psi$ unique on $\mathcal F\vert_D$; needs no structure on $\mathcal X$ |
| **A1** | $\omega_{x,D}(0^+)=0$ | stable inference (error $\to0$ with $\varepsilon$); optimal modulus is exactly $\varepsilon\mapsto\tfrac12\omega(2\varepsilon)$ | strictly stronger than A0: tanh counterexample (§7 below) |
| **A2** | $\mathcal F$ inside a linear space of dimension $d\le k$ with a rank-$d$ design | **zero-error** recovery in the linear category | iff: exact identifiability of a linear family at a $k$-point design forces $\dim\le k$; NOT necessary for stable approximate inference (RKHS ball satisfies A1 at every design with infinite dimension) |

**Design-dependence is unavoidable in multivariate continuous settings.** On intervals, $d$-dimensional Haar spaces make *every* $d$-point design unisolvent. But (Mairhuber–Curtis–Sieklucki, Thm 2.2 + Prop. 2.3): on any space containing a triod — in particular open subsets of $\mathbb R^n$, $n\ge2$ — no $d\ge2$-dimensional space of continuous functions is unisolvent at all designs. So above dimension one, identifiability is a **joint property of family and design**, and the rank conditions of §4 must be checked at the actual $D$.

**For generic-design guarantees in nonlinear parametric families, the minimal regularity is quasianalyticity** (o-minimal definability + identity-theorem property for differences), not smoothness: Thm 8.1 requires it and the flat-$C^\infty$ counterexample (§7) shows smoothness alone yields no genericity statement of any kind.

---

## 4. Main theorem and proposition statements

- **Thm 1 (exact minimax).** Minimax error $=\tfrac12\omega_{x,D}(2\varepsilon)$; midpoint rule optimal; survives randomization.
- **Prop. 1.5 (identifiability boundary).** Exact identifiability $\iff\omega(0)=0\iff$ graph property; $\Psi$ unique on realizable traces; every zero-error rule equals $\Psi$ there.
- **Thm 2.2 (Mairhuber–Curtis–Sieklucki).** For compact Hausdorff $\mathcal X$, $|\mathcal X|\ge d$, $d\ge2$: a $d$-dim Haar subspace of $C(\mathcal X)$ exists iff $\mathcal X$ embeds as a closed subset of $S^1$, with $d$ odd if $\mathcal X\cong S^1$.
- **Prop. 2.3 (triod obstruction).** Direct swap-argument proof, no compactness: any triod-containing domain kills design-independent unisolvence for $d\ge2$ (continuous functions only; discontinuous pullback escape exists).
- **Thm 4.1 (summary dimension).** $d\le m_{\min}\le 2d+1$ for continuous member summaries (invariance of domain / Menger–Nöbeling); equality $m_{\min}=d$ iff the parameter space embeds in $\mathbb R^d$; measurable summaries collapse to dimension 1; infinite-dimensional families admit no finite continuous summary.
- **Thm 5.1 (rank theorem).** Linear functional $\lambda^\top c$ identifiable iff $\lambda\in\operatorname{row}(G)$; value $f_\beta(x)$ identifiable iff $\phi(x)\in\operatorname{row}(G)$; full member iff $\operatorname{rank}G=d$; constrained set $C$: iff $\ker G\cap(C-C)=\{0\}$.
- **Thm 5.2 (exact minimax, linear).** Minimax error $=\varepsilon\Lambda_*(x)$, $\Lambda_*(x)=\min\{\|w\|_*:G^\top w=\phi(x)\}$; attained by a linear rule (Smolyak 1965 / Marchuk–Osipenko 1975 lineage); Hahn–Banach duality proof.
- **Cor. 5.3 (sensitivity).** Optimal weights $\hat w=(G^+)^\top\phi(x)$; $\Lambda_2\le\|\phi(x)\|/\sigma_{\min}(G)$; weights are the sensitivity profile; weights sum to 1 when constants lie in $V$.
- **Thm 5.4 (nonlinear local rank).** Jacobian rank $d$ at $\theta_0$ ⇒ local identifiability with stability $\|J^+\|$; constant rank $r<d$ ⇒ $(d-r)$-dimensional local ambiguity manifolds.
- **Thm 7.1 (relative values, linear).** Difference identifiable iff $\phi(x_a)-\phi(x_b)\in\operatorname{row}(G)$; possible with neither value identifiable; joint value identifiability strictly stronger than difference identifiability.
- **Thm 7.2 (gauge taxonomy).** Additive family $g_0+\mathbb R\mathbf1$: all differences identifiable at $k=0$, values at $k=1$. Multiplicative $\{ag_0\}$: ratios at $k=0$, values at $k=1$. Affine $\{ag_0+b\}$, $\{1,g_0\}$ independent: affine invariants at $k=0$; at $k=1$ differences iff $g_0(x_a)=g_0(x_b)$; everything at $k=2$ with $g_0(x_1)\ne g_0(x_2)$.
- **Thm 8.1 (generic global identifiability).** Real-analytic separating $d$-parameter family on bounded connected domain: for $k\ge2d+1$, all designs outside a **closed** null subanalytic set give global injectivity of $\theta\mapsto(f_\theta(x_i))$.
- **Thm 8.2 (stable version).** Adding infinitesimal separation: generic $k\ge2d+1$ designs give a bi-Lipschitz evaluation map on compact convex subsets, hence $\|\hat\theta-\theta\|\le2\varepsilon\sqrt k/c(K)$.
- **Thm 8.3 (necessity).** Continuous exact recovery of a family containing a $d$-cell forces $k\ge d$. With $k\le5$: exact stable recovery only for $d\le5$; generic-global guarantee only for $d\le2$.
- **Thm 9.1 (archive identification).** Common unisolvent core $X_0$, rank-$d$ archive matrix: $V|_U$ identified iff the per-point rank conditions hold at every $x\in U$ (necessary and sufficient within the class of exactly-$d$-dimensional spaces, $d$ known); ambiguity otherwise is an affine subspace of dimension $d-\operatorname{rank}$, fully realized.
- **Thm 9.2 (archive impossibility dichotomy).** Off the covered set: consistent value set is $\mathbb R$ or (degenerately) $\{0\}$.
- **Derived procedure (§10).** Archive → trace sets on $U$; sample → consistent section; output → Chebyshev center. Forced by Thm 1 given its premises; linear/RKHS/Lipschitz closed forms (min-dual-norm weights, kernel interpolant, envelope average) drop out rather than being assumed; misspecified data handled by projection with a doubling bound $\tfrac12\omega(2\varepsilon+2\eta)$.

---

## 5. The exact identifiable object

**General family, exact data:** precisely the functionals $\varphi:\mathcal F\to\mathbb R$ that **factor through the trace map** $f\mapsto f|_D$ (constant on its fibers). Closed under arbitrary composition. With noise: $\varphi$ is stably identifiable iff $\omega_\varphi(0^+)=0$, with optimal error exactly $\tfrac12\omega_\varphi(2\varepsilon)$.

**Linear family $V$, basis $\{\phi_j\}$, evaluation matrix $G$:**

- identifiable *linear* functionals of $c_\beta$ = $\operatorname{row}(G)=\operatorname{Ann}(N_D)$, where $N_D=\{v\in V:v|_D=0\}$ — an $r$-dimensional space, $r=\operatorname{rank}G\le k$;
- the identified part of the member is the projection $P_{\operatorname{row}(G)}\,c_\beta$ ($r$ numbers);
- identifiable queries: exactly $\{x:\phi(x)\in\operatorname{row}(G)\}$, with value $\phi(x)^\top G^+\tilde y$;
- the **minimal data-computable summary** $z_\beta$ has dimension exactly $r=\operatorname{rank}G$; its existence is *equivalent to identifiability*, never an extra assumption ($z_\beta=\tilde y$ works whenever anything works).

**Family-level (common) information:** the trace sets $T_S$ for finite $S$ — and an archive can reveal them **only for $S$ inside the covered set $U$** (Thm 9.1 gives the exact conditions; Thm 9.2 the exact limit).

---

## 6. The exact non-identifiable components

- **General family:** the fiber partition of the trace map — for values at $x$: the section $\{t:(y,t)\in T_{D,x}\}$, of worst-case diameter $\omega_{x,D}(0)$ ($=2\rho$, twice the radius of information). Nonzero diameter is *the* obstruction; it can be finite (partial information, e.g. Lipschitz class: diameter $2\min_i\operatorname{dist}(x,x_i)$) or infinite (no information).
- **Linear family:** the coset $c_\beta+\ker G$ restricted to $N_D$ — dimension $d-r$ — which corrupts exactly the queries with $\phi(x)\notin\operatorname{row}(G)$. Residual gauge after $k$ observations = stabilizer of the observed trace within the family's symmetry group; identifiable functionals are its invariants (relative values survive precisely the additive gauge).
- **Family level:** everything at queries off the covered set $U$ (Thm 9.2 dichotomy) — non-identifiability recurs one level up, and no volume of archive data at other points repairs it.
- **Without continuity of the summary:** dimension itself is non-identifiable as a complexity measure (Borel collapse to dimension 1) — minimal-dimension statements are meaningful only in the continuous/stable category.

---

## 7. All important counterexamples

| # | Counterexample | What it defeats |
|---|---|---|
| C1 | $\mathcal F=\mathbb R^{\mathcal X}$ | any inference off $D$; information exactly on $D$ (Cor. 1.4) |
| C2 | 1-Lipschitz class; $S(y)$ = McShane–Whitney interval; $\omega(0)=2\min_i\operatorname{dist}(x,x_i)$ | identifiability without partial information being zero — informative $\ne$ identifying |
| C3 | $f_\theta(x_1)=\tanh\theta$, $f_\theta(x)=\theta$ | A0 ⇒ A1: identifiable but $\omega(t)=\infty$ for all $t>0$ — noise destroys everything |
| C4 | $\sin(\theta x)$, $\Theta=\mathbb R$, $d=1$: $k=1$ fails at every design; $k=2$ fails at every design via $\theta=\frac{\pi}{2b}-\frac{\pi}{a},\ \theta'=\frac{\pi}{2b}+\frac{\pi}{a}$; $k=3$ works at pairwise-irrational-ratio designs | sharpness of the $2d+1$ threshold; the linear-vs-nonlinear gap ($k=d$ vs $k=2d+1$) |
| C5 | $f_\theta(x)=\theta\,e^{-1/x}\mathbf 1_{x>0}$ — $C^\infty$, o-minimally definable, separating | genericity under mere smoothness/definability: bad designs have nonempty interior for every $k$; quasianalyticity is the true dividing line |
| C6 | tripod parameter space | $m_{\min}=d$ for summaries: here $m_{\min}=2>d=1$; only the sandwich $d\le m_{\min}\le2d+1$ holds |
| C7 | $n$-cube of tent functions inside the Lipschitz class | any finite-dimensional continuous summary of an infinite-dimensional family |
| C8 | Borel bijection $\mathbb R^d\to\mathbb R$ | dimension lower bounds without continuity |
| C9 | $H(x_1)=(0,0),H(x_2)=(1,0),H(x_a)=(0,1),H(x_b)=(1,1)$ in $V=\operatorname{span}\{1,h_1,h_2\}$ | "values needed for differences": the difference is identifiable while neither value is |
| C10 | $S^1$ with even $d$ (cyclic-shift sign argument) | the naive Mairhuber–Curtis biconditional; parity exclusion required |
| C11 | polynomials pulled back through a discontinuous bijection $[0,1]^2\to[0,1]$ | reading the triod obstruction as set-theoretic; it is topological (continuity load-bearing) |
| C12 | Lipschitz class, $D_1=\{-1\}$, $D_2=\{+1\}$, query $0$ | "equal radius ⇒ equal information": $\rho$ comparison is one-directional |
| C13 | forced-zero archive case ($y=0$, full-rank design) | the blanket "uncovered queries are completely unconstrained"; exact statement is the $\mathbb R$-or-$\{0\}$ dichotomy |
| C14 | two-branch family $\mathcal F_A\cup\mathcal F_B$ on disjoint supports | "adaptive design never helps": it helps unboundedly for nonconvex families (while gaining nothing for convex balanced ones) |
| C15 | equispaced vs Chebyshev nodes: $\Lambda\sim2^k/(ek\log k)$ vs $\sim\frac2\pi\log k$ | "any unisolvent design is as good as any other": designs of equal size differ exponentially |

---

## 8. All assumptions required by the theory

**Global (every result).**
- $\mathcal F$ nonempty; values in $\mathbb R$ (scalar target is load-bearing for the exact constant ½; sup-loss multi-query version holds, other joint losses only up to radius/diameter gap).
- Noise model: per-coordinate error, **closed** ball, known level $\varepsilon$; deterministic estimators (randomization covered by the lower bound only).
- Exact-data statements ($\varepsilon=0$) marked as such; the optimal rule depends on $\varepsilon$.

**Theorem 1 / Prop. 1.5:** nothing further. (This is the point.)

**Haar / Mairhuber–Curtis block:** the zero-count ⇔ nonsingularity equivalence needs only a $d$-dimensional space on a set with $\ge d$ points; Thm 2.2 needs $\mathcal X$ compact Hausdorff, $|\mathcal X|\ge d$, $d\ge2$, with the circle-parity exclusion; Prop. 2.3 needs continuity of the members and an embedded triod; exponential-family example needs distinct **real** exponents.

**Linear theory (Thms 5.1–5.3, 7.1):** $\{\phi_j\}$ linearly independent as functions; for the exact minimax formula: $\mathcal F$ = all of $V$ (unconstrained coefficients — for constrained $C$ the formula is an upper bound and the criterion becomes $\ker G\cap(C-C)=\{0\}$), $V$ known exactly (archive error adds $\delta\cdot C(x,G)\|c_\beta\|$ via subspace-angle perturbation), identifiability $\phi(x)\in\operatorname{row}(G)$ (otherwise both sides $+\infty$), linear scalar target, convex balanced noise ball.

**Nonlinear local theory (Thm 5.4):** $C^1$ parametrization on open $\Theta\subseteq\mathbb R^d$; constant-rank hypothesis for the negative direction.

**RKHS block:** positive-definite kernel; $K_{DD}$ nonsingular (else pseudoinverse/variational definitions); feasibility $\|s_y\|_{\mathcal H}\le1$ for the Golomb–Weinberger interval; the GP identity is an identity of formulas, not models (paths a.s. outside $\mathcal H$ when infinite-dimensional — Driscoll).

**Generic identifiability (Thms 8.1–8.2):** $\Theta\subset\mathbb R^d$, $X\subset\mathbb R^m$ open **bounded**, $X$ **connected** (componentwise separation if not), $F$ real-analytic **on a neighborhood of $\overline\Theta\times\overline X$** (boundedness and closure-analyticity are load-bearing for subanalytic projection theory); parametrization **separating**; for Thm 8.2 additionally **infinitesimally separating**, and constants degrade toward $\partial\Theta$ (compact convex $K$). Minimal known relaxation: o-minimal definability + empty-interior agreement sets (quasianalytic Denjoy–Carleman classes qualify).

**Necessity bounds (Thm 8.3, Thm 4.1):** the family carries a topology at least as fine as pointwise convergence and contains a $d$-cell; summary maps continuous (else C8 collapses dimension); cardinality of the query-quotient $\le\mathfrak c$ for any finite-dimensional summary to exist at all.

**Archive theorems (9.1–9.2):** archive observed exactly ($\varepsilon_{\text{archive}}=0$; noisy archive needs the subspace-angle addendum); model class = subspaces of dimension **exactly** $d$ with $d$ known (rank certifies $\ge d$ only, never $\le d$); common core $X_0\subseteq\bigcap_j D_j$ with $|X_0|=d$ and rank-$d$ archive matrix (sufficient route — overlap chaining can substitute); per-point rank conditions for the iff; Thm 9.2 stated for $x\notin U\cup D$.

**Budget consequence of $k\le5$:** exact stable recovery requires family dimension $d\le5$; generic-global identifiability for nonlinear analytic families is guaranteed only for $d\le2$; for $d\in\{3,4,5\}$ nonlinear families, identifiability must be established per-design (no generic guarantee at $k\le5$), and for linear families rank-$d$ designs exist for every $d\le5$ on any set.

---

*End of frozen foundation. Next stages (model development, application mapping) must cite items here by number and may not weaken assumptions silently.*
