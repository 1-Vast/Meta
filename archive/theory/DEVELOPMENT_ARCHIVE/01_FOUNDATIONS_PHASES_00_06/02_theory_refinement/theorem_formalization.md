# Formal Statement of the Identifiability Theory

> **Status:** theory-refinement layer, 2026-08-02. Source of record: `../00_raw_outputs/identifiability_treatise_raw.md`; frozen summary: `../01_math_foundation/theorem_summary.md`. This file restates the theory axiomatically: standing assumptions, definitions, formal statements with type tags (**[iff]**, **[sufficient]**, **[necessary]**, **[exact value]**, **[lower bound]**, **[dichotomy]**, **[counterexample]**), proof sketches, and the counterexample registry. Counterexample labels C1–C15 refer to §7 of the frozen summary and are reproduced in Part VIII.

---

## Part 0 — Standing assumptions

- **SA1.** $\mathcal X$ is a set; no topology, metric, measure, or representation is assumed unless a statement introduces it explicitly. $\mathcal F\subseteq\mathbb R^{\mathcal X}$ is **nonempty**.
- **SA2 (observation model).** A design is a finite $D=\{x_1,\dots,x_k\}\subset\mathcal X$ of distinct points, $k\le5$. Data: $\tilde y\in\mathbb R^k$ with $\|\tilde y-f_\beta|_D\|_\infty\le\varepsilon$ for the true member $f_\beta\in\mathcal F$; the noise ball is **closed**, the level $\varepsilon\ge0$ is **known**. $\varepsilon=0$ is "exact data".
- **SA3 (estimators).** An estimator is an arbitrary map $\Phi:\mathbb R^k\to\mathbb R$ (no measurability required). Randomized estimators are covered by lower bounds only (F1, Rem.).
- **SA4 (target).** The recovery target is a scalar: $f_\beta(x)$ at a query $x\notin D$, or more generally a functional $\varphi(f_\beta)$. Scalar-valuedness is load-bearing for the exact constant in F1 (Rem. F1.2).

---

## Part I — Definitions

- **D1 (trace map).** $R_D:\mathcal F\to\mathbb R^k$, $R_D f=f|_D=(f(x_1),\dots,f(x_k))$.
- **D2 (trace set / window).** $T_{D,x}=\{(f|_D,f(x)):f\in\mathcal F\}\subseteq\mathbb R^{k+1}$; more generally $T_S=\{f|_S:f\in\mathcal F\}$ for finite $S\subset\mathcal X$.
- **D3 (section).** $S_\varepsilon(\tilde y)=\{t\in\mathbb R:\exists(u,t)\in T_{D,x},\ \|u-\tilde y\|_\infty\le\varepsilon\}$. Write $S(\tilde y)$ when $\varepsilon$ is fixed; $S_0$ is the exact-data section map.
- **D4 (realizable data).** $\tilde y$ is realizable iff $S_\varepsilon(\tilde y)\ne\emptyset$ iff some member is $\varepsilon$-consistent with $\tilde y$.
- **D5 (trace modulus).** $\omega_{x,D}(t)=\sup\{|f(x)-g(x)|:f,g\in\mathcal F,\ \|f|_D-g|_D\|_\infty\le t\}$; for a functional $\varphi:\mathcal F\to\mathbb R$, $\omega_\varphi(t)=\sup\{|\varphi(f)-\varphi(g)|:\|f|_D-g|_D\|_\infty\le t\}$. Both are nondecreasing, $[0,\infty]$-valued, $\omega(0)\ge0$ by SA1.
- **D6 (identifiability).** $f_\beta(x)$ is *exactly identifiable at $(D,x)$* iff $f|_D=g|_D\Rightarrow f(x)=g(x)$ on $\mathcal F$ (equivalently $\omega_{x,D}(0)=0$). A functional $\varphi$ is identifiable iff it is constant on the fibers of $R_D$. *Stable* identifiability: worst-case error $\to0$ as $\varepsilon\to0$.
- **D7 (radius of information).** $\rho(x;D,\varepsilon)=\tfrac12\,\omega_{x,D}(2\varepsilon)$.
- **D8 (linear data).** $V=\operatorname{span}\{\phi_1,\dots,\phi_d\}$ with $\{\phi_j\}$ linearly independent as functions; $G\in\mathbb R^{k\times d}$, $G_{ij}=\phi_j(x_i)$; feature vector $\phi(x)=(\phi_1(x),\dots,\phi_d(x))^\top$; null set $N_D=\{v\in V:v|_D=0\}$; coefficient set $C\subseteq\mathbb R^d$ when $\mathcal F=\{f_c:c\in C\}$.
- **D9 (generalized Lebesgue function).** For a norm $\|\cdot\|$ on $\mathbb R^k$ with dual $\|\cdot\|_*$: $\Lambda_*(x)=\min\{\|w\|_*:G^\top w=\phi(x)\}$ ($=+\infty$ if infeasible).
- **D10 (Haar space).** $d$-dimensional $V$ of functions on a set with $\ge d$ points such that no nonzero member has $d$ zeros; equivalently (pure linear algebra) every $d$-point collocation matrix is nonsingular.
- **D11 (power function).** For a positive-definite kernel $K$ with RKHS $\mathcal H$ and nonsingular $K_{DD}$: $P_D(x)^2=K(x,x)-k_x^\top K_{DD}^{-1}k_x$.
- **D12 (parametric regularity).** A parametrization $\theta\mapsto f_\theta$, $\theta\in\Theta$, is *separating* iff injective into $\mathbb R^{\mathcal X}$ (componentwise on disconnected domains); *infinitesimally separating* iff no $\theta$ and $u\ne0$ give $\partial_uf_\theta\equiv0$.
- **D13 (archive).** A collection $\{(g_j,D_j)\}_{j\le n}$ of members with their designs, values known exactly on $D_j$; covered set $U=\bigcup_jD_j$.
- **D14 (equivalence of problems).** Two families induce *the same inference problem at $(D,x)$* iff for every $\varepsilon\ge0$ they have identical section maps $S_\varepsilon$ (equivalently, by F-O9 in `operator_formulation.md`, identical trace sets $T_{D,x}$).

---

## Part II — The fundamental theorem and the identifiability boundary

**F1 (Theorem — exact minimax). [exact value]**
*Hypotheses:* SA1–SA4 only.
*Statement:* $\displaystyle\inf_\Phi\ \sup_{f\in\mathcal F}\ \sup_{\|\tilde y-f|_D\|_\infty\le\varepsilon}|\Phi(\tilde y)-f(x)|\;=\;\tfrac12\,\omega_{x,D}(2\varepsilon)$ in $[0,\infty]$.
*Proof sketch:* Lower — a pair at trace distance $\le2\varepsilon$ shares the coordinatewise-midpoint data; any $\Phi$ errs $\ge\tfrac12|f(x)-g(x)|$ on one. Upper — for realizable $\tilde y$, $\operatorname{diam}S(\tilde y)\le\omega(2\varepsilon)$ (triangle inequality through $\tilde y$), with $\sup_{\tilde y}\operatorname{diam}S(\tilde y)=\omega(2\varepsilon)$; the midpoint rule $\Phi=\tfrac12(\inf S+\sup S)$ (arbitrary off the realizable set) attains $\tfrac12\omega(2\varepsilon)$ when finite; if infinite, the lower bound closes the identity. $\square$
*Rem. F1.1:* the noise argument $2\varepsilon$ is sharp: a shared $\varepsilon$-consistent datum exists for a pair iff their traces are within $2\varepsilon$.
*Rem. F1.2:* the constant $\tfrac12$ uses target $=\mathbb R$ (radius $=$ half-diameter on the line). Multi-query: under sup-loss the per-query midpoint rule is again exactly optimal with value $\sup_x\tfrac12\omega_{x,D}(2\varepsilon)$; under other joint losses only radius $\in[\tfrac12\operatorname{diam},\operatorname{diam}]$ holds.
*Rem. F1.3:* the lower bound survives randomization: $\max(\mathbb E|Z-a|,\mathbb E|Z-b|)\ge\tfrac12|a-b|$.
*Rem. F1.4:* the optimal rule depends on $\varepsilon$; the optimal error modulus is $\varepsilon\mapsto\tfrac12\omega(2\varepsilon)$.

**F2 (Proposition — identifiability boundary). [iff]**
Exact-data identifiability at $(D,x)$ $\iff$ $\omega_{x,D}(0)=0$ $\iff$ $T_{D,x}$ is the graph of a function $\Psi$ over $\mathcal F|_D$. Then $\Psi$ is unique on $\mathcal F|_D$ and every zero-error rule equals $\Psi$ there; off $\mathcal F|_D$ no optimality principle constrains a rule.
*Proof sketch:* all three are restatements of "$f|_D=g|_D\Rightarrow f(x)=g(x)$"; uniqueness because $\Psi$'s values are forced pointwise. $\square$

**F3 (Proposition — stability boundary). [iff]**
Stable identifiability $\iff\omega_{x,D}(0^+)=0$; Lipschitz-$L$ inference possible $\iff\omega(t)\le2Lt$ for all $t$. Strictness against F2: **C3** ($\tanh$ family: $\omega(0)=0$, $\omega(t)=\infty\ \forall t>0$).

**F4 (Corollary — absolute impossibility). [lower bound / $\infty$]**
If $\mathcal F$ contains, for every $M$, a pair agreeing on $D$ and differing by $\ge M$ at $x$ (e.g. $\mathcal F=\mathbb R^{\mathcal X}$, **C1**), every estimator — randomized included — has infinite worst-case error at $x$; for $\mathcal F=\mathbb R^{\mathcal X}$ the minimax error is $\varepsilon$ on $D$ and $+\infty$ off $D$: finite error exactly on $D$.

**F5 (Proposition — partial information). [exact value]**
For the 1-Lipschitz class on a metric space with exact data: $S_0(y)=[\max_i(y_i-\operatorname{dist}(x,x_i)),\ \min_i(y_i+\operatorname{dist}(x,x_i))]$ (feasibility: $|y_i-y_j|\le\operatorname{dist}(x_i,x_j)$), $\omega_{x,D}(0)=2\min_i\operatorname{dist}(x,x_i)$, minimax $=\min_i\operatorname{dist}(x,x_i)$. **[counterexample role]** identifiable $\ne$ informative: $\omega(0)>0$ yet finite (**C2**).

---

## Part III — Linear theory

**F6 (Theorem — rank). [iff, three forms]**
*Hypotheses:* D8, exact data $y=Gc_\beta$.
(i) Unconstrained ($C=\mathbb R^d$): $\lambda^\top c_\beta$ determined by $y$ $\iff\lambda\in\operatorname{row}(G)$, and then equals $\lambda^\top G^+y$ for every consistent $c$. (ii) $f_\beta(x)$ identifiable $\iff\phi(x)\in\operatorname{row}(G)$; the whole $c_\beta$ $\iff\operatorname{rank}G=d$ (forces $k\ge d$). (iii) Constrained: the *member* ($c_\beta$, equivalently all queries simultaneously) is identifiable on $C\iff\ker G\cap(C-C)=\{0\}$; per-query variant: $f_\beta(x)$ identifiable on $C\iff\phi(x)^\top$ vanishes on $\ker G\cap(C-C)$; a general functional is identifiable iff constant on each fiber $\{c\in C:Gc=y\}$.
*Proof sketch:* (i) $\lambda\perp\ker G$ is exactly constancy on fibers $c+\ker G$; the pseudoinverse formula is fiber-independent. (iii) two consistent coefficients differ by an element of $\ker G\cap(C-C)$. $\square$
*Sharpness:* identifiability below full rank is real for constrained/nonlinear $C$ (F11).

**F7 (Theorem — exact minimax, linear). [exact value; sufficiency of a linear rule]**
*Hypotheses:* D8–D9; $\mathcal F=V$ **unconstrained**; $V$ known exactly; closed noise ball in an arbitrary norm $\|\cdot\|$; $\phi(x)\in\operatorname{row}(G)$ (else both sides $+\infty$); scalar linear target.
*Statement:* minimax error $=\varepsilon\Lambda_*(x)$, attained by the linear rule $\tilde y\mapsto\hat w^\top\tilde y$ at any minimizer $\hat w$ of D9 (its worst-case error is exactly $\varepsilon\|\hat w\|_*$).
*Proof sketch:* upper — $|\hat w^\top(\tilde y-Gc)|\le\|\hat w\|_*\|\tilde y-Gc\|$. Lower — two-point argument (midpoint valid for any norm by homogeneity) plus the Hahn–Banach identity: the norm of $v\mapsto w_0^\top v$ on $\operatorname{col}(G)$ equals $\min\{\|w\|_*:G^\top w=G^\top w_0\}$; extension attaining the min exists in finite dimension. $\square$
*Attribution:* Smolyak (1965, exact information); Marchuk–Osipenko (1975, noisy values); fails for jittered nodes (Kacewicz–Plaskota 2003).
*Caution:* for constrained $C$ the formula is an upper bound only; canonicity of the linear rule is **not** claimed (see `operator_formulation.md`, the midrange example).

**F8 (Corollary — sensitivity). [exact identities]**
$\hat w=(G^+)^\top\phi(x)$ (any rank, given feasibility); $\Lambda_2(x)=\|(G^+)^\top\phi(x)\|_2\le\|\phi(x)\|_2/\sigma_{\min}(G)$; $\partial(\hat w^\top\tilde y)/\partial\tilde y_i=\hat w_i$ (the weights are the sensitivity profile); if $\mathbf1\in V$, $\sum_iw_i=1$ for every feasible $w$.

**F9 (Proposition — dimension is the price of exactness). [necessary]**
A linear family exactly identifiable at some $k$-point design for all queries has $\dim\le k$ (the restriction map is injective linear into $\mathbb R^k$). *Not necessary for stability:* the RKHS unit ball satisfies F3 at every design with $\dim=\infty$; its radius is the Golomb–Weinberger interval half-length $P_D(x)\sqrt{1-\|s_y\|_{\mathcal H}^2}$ (D11).

**F10 (Theorem — design-independence). [sufficient / impossibility]**
(i) **[sufficient]** Haar (D10) $\Rightarrow$ every $d$-point design identifies every query with exact data; on intervals: polynomials of degree $<d$, real-exponent exponential spans.
(ii) **[impossibility]** (Mairhuber–Curtis–Sieklucki) compact Hausdorff $\mathcal X$, $|\mathcal X|\ge d$, $d\ge2$: a $d$-dim Haar subspace of $C(\mathcal X)$ exists iff $\mathcal X$ embeds as a closed subset of $S^1$, with $d$ odd when $\mathcal X\cong S^1$ (**C10**). Direct triod version (no compactness): on any triod-containing domain, every $d\ge2$-dim space of continuous functions has some **size-$d$** design at which evaluation is singular (swap argument).
(iii) *Scope:* three restrictions keep (ii) honest. The obstruction concerns designs of size **exactly $d$**: with $k>d$ points, design-independent identifiability of a $d$-dimensional continuous family on a triod is possible — take $V=\operatorname{span}\{1,f\}$ with $f$ continuous on the triod assuming each value at most twice (e.g. arclength $s$ on one arm, $2s$ on the second, $-s$ on the third); then among any three distinct points two have different $f$-values, so every $3$-point design has rank-$2$ evaluation and identifies all queries. The obstruction is topological, not set-theoretic (**C11**: discontinuous pullback restores design-independence on the square); and $d=1$ is exempt. At $k=d$, the Haar property is *equivalent* to every-$d$-point-design identifiability (D10), so (i) is an iff at that design size.
(iv) Node quality: with $k$ nodes on $[-1,1]$, $\Lambda\sim\frac2\pi\log k$ (Chebyshev) vs $\Lambda\sim2^k/(ek\log k)$ (equispaced) (**C15**).

---

## Part IV — Parametric (nonlinear) families

**F11 (Theorem — local rank). [sufficient local / necessary local]**
$C^1$ family on open $\Theta\subseteq\mathbb R^d$, Jacobian $J(\theta)=(\partial_\theta f_\theta(x_i))\in\mathbb R^{k\times d}$. If $\operatorname{rank}J(\theta_0)=d$: local identifiability near $\theta_0$ with $\|\delta\theta\|\lesssim\|J(\theta_0)^+\|\|\delta y\|$. If $\operatorname{rank}J\equiv r<d$ near $\theta_0$: fibers are $(d-r)$-manifolds (constant-rank theorem) — locally non-identifiable, ambiguity tangent to $\ker J$.

**F12 (Theorem — generic global identifiability). [sufficient, generic]**
*Hypotheses:* $\Theta\subset\mathbb R^d$, $X\subset\mathbb R^m$ open **bounded**, $X$ **connected**; $F$ real-analytic on a neighborhood of $\overline\Theta\times\overline X$; separating (D12).
*Statement:* for $k\ge2d+1$ there is a **closed** subanalytic $N\subset X^k$, $\dim N\le2d+k(m-1)<km$ (Lebesgue-null), such that every design in the open dense full-measure set $X^k\setminus N$ makes $\theta\mapsto(f_\theta(x_i))_i$ injective on all of $\Theta$.
*Proof sketch:* agreement sets of distinct parameters are proper analytic subsets of connected $X$ (identity theorem), $\dim\le m-1$; the incidence set over parameter pairs is relatively compact semianalytic with fiber dimension $\le k(m-1)$, total $\le2d+k(m-1)$; bounded-subanalytic projections do not raise dimension; the o-minimal frontier theorem closes the bad set. $\square$
*Load-bearing:* boundedness; connectedness; closure-analyticity. Minimal relaxation: o-minimal definability $+$ empty-interior agreement sets (quasianalytic Denjoy–Carleman qualifies).

**F13 (Theorem — stable generic version). [sufficient, generic]**
Add infinitesimal separation (D12): for $k\ge2d+1$, generic designs give an evaluation map that is injective and an immersion; on each compact convex $K\subset\Theta$, $\|\Phi(\theta)-\Phi(\theta')\|\ge c(K)\|\theta-\theta'\|$, hence $\|\hat\theta-\theta\|\le2\varepsilon\sqrt k/c(K)$ for $\varepsilon$-consistent estimates. (Immersion alone is generic at $k\ge2d$; the combined guarantee carries $2d+1$.)

**F14 (Proposition — necessity, sharpness, dividing line).**
(i) **[necessary]** exact recovery — by *any* map, continuity of the rule not needed — of a family containing a continuously and injectively parametrized $d$-cell forces $k\ge d$: exact recovery makes the evaluation map injective on the cell, and it is continuous by hypothesis on the parametrization, so a continuous injection $[0,1]^d\to\mathbb R^k$ exists, impossible for $k<d$ by invariance of domain.
(ii) **[sharpness]** (**C4**) $f_\theta(x)=\sin(\theta x)$, $\Theta=\mathbb R$, $d=1$: $k=1$ fails at every design; $k=2$ fails at every design via the explicit pair $\theta=\frac{\pi}{2b}-\frac\pi a,\ \theta'=\frac{\pi}{2b}+\frac\pi a$; $k=3=2d+1$ succeeds at designs with pairwise irrational ratios (pigeonhole over translation/reflection cases). The $2d+1$ threshold is attained; linear families need only $k=d$.
(iii) **[counterexample]** (**C5**) $f_\theta(x)=\theta e^{-1/x}\mathbf1_{x>0}$: $C^\infty$, o-minimally definable, separating — bad designs have nonempty interior for every $k$. Quasianalyticity, not smoothness, is the dividing line.
(iv) Budget: $k\le5\Rightarrow$ exact stable recovery only for $d\le5$; F12's guarantee only for $d\le2$.

---

## Part V — Functional identifiability and relative values

**F15 (Theorem — relative values, linear). [iff; witness]**
$f_\beta(x_a)-f_\beta(x_b)$ identifiable $\iff\phi(x_a)-\phi(x_b)\in\operatorname{row}(G)$ — possible with neither value identifiable (**C9**). For $V=\operatorname{span}\{1,h_1,h_2\}$, $k=2$: the criterion is $H(x_a)-H(x_b)\in\operatorname{span}\{H(x_1)-H(x_2)\}$ (zero multiple allowed; degenerate $H(x_1)=H(x_2)$ handled separately). Joint identifiability of both values is strictly stronger than of the difference; given the difference identifiable, $\phi(x_a)\in\operatorname{row}(G)\iff\phi(x_b)\in\operatorname{row}(G)$. Quantitatively $\omega_{\delta_a-\delta_b}$ can vanish while $\omega_{\delta_a}=\omega_{\delta_b}=\infty$.

**F16 (Proposition — gauge taxonomy). [exact observation counts]**
Convention: a partially defined functional is identifiable iff constant on (fiber $\cap$ domain).
(a) $\mathcal F=g_0+\mathbb R\mathbf1$: all differences at $k=0$ (at archive-covered pairs), all values at $k=1$.
(b) $\mathcal F=\{ag_0\}$: ratios at $k=0$ where defined (definedness at the true member itself unidentifiable at $k=0$); values at $k=1$ if $g_0(x_1)\ne0$.
(c) $\mathcal F=\{ag_0+b\}$, $\{1,g_0\}$ independent: affine invariants at $k=0$; at $k=1$ differences iff $g_0(x_a)=g_0(x_b)$ (value $0$); everything at $k=2$ if $g_0(x_1)\ne g_0(x_2)$. If $\{1,g_0\}$ dependent, the family degenerates to constants and $k=0$ identifies all differences as $0$.
*Principle:* residual ambiguity $=$ stabilizer of the observed trace in the family's symmetry; identifiable functionals $=$ its invariants.

---

## Part VI — The archive level

**F17 (Theorem — archive identification). [iff within the fixed-$d$ class]**
*Hypotheses:* D13, archive exact; model class $=$ subspaces of dimension **exactly** $d$, $d$ known; common core $X_0\subseteq\bigcap_jD_j$, $|X_0|=d$; archive matrix $A=(g_j|_{X_0})$ of rank $d$ (data-checkable).
*Statement:* restriction $V\to\mathbb R^{X_0}$ is bijective and the $g_j$ span $V$; for $x\in U$, the representing row $a_x$ (with $g(x)=a_x\cdot g|_{X_0}$ on $V$) is determined **iff** the submatrix of $A$ over members observed at $x$ has rank $d$; otherwise the ambiguity is an affine subspace of dimension $d-\operatorname{rank}$, every point realized by a genuine candidate family. All conditions holding on $U$ $\Rightarrow$ $V|_U$ fully identified; new-member inference on $U$ proceeds by F6 with rows $a_{x_i}$. Common core is one sufficient route (overlap chaining can substitute).

**F18 (Theorem — archive impossibility). [dichotomy]**
For $x\notin U\cup D$: the set of values $f_\beta(x)$ consistent with archive and sample, over all exactly-$d$-dimensional candidate families, is **either $\mathbb R$ or $\{0\}$** — the latter precisely when the data force $f_\beta=0$ (**C13**; the exception vanishes for $\varepsilon>0$).
*Proof sketch:* modify a basis of $V$ at the single point $x$ arbitrarily: restrictions to $U$ unchanged, archive and sample matched, $f_\beta(x)=c_\beta\cdot t$ sweeps $\mathbb R$ when $c_\beta\ne0$. $\square$
*Meaning:* F4 recurs at family level; only member-level assumptions transport information off $U$.

---

## Part VII — Summary dimension

**F19 (Theorem — dimension of member summaries). [necessary / sufficient / counterexamples]**
For $\mathcal F$ injectively, bicontinuously parametrized by $\Theta$ (pointwise topology):
(i) **[necessary]** $\Theta\supseteq$ open subset of $\mathbb R^d$ $\Rightarrow$ any continuous injective $z:\mathcal F\to\mathbb R^m$ has $m\ge d$ (invariance of domain).
(ii) **[sufficient]** $\Theta$ compact metric of covering dimension $d$ $\Rightarrow$ some continuous injective $z$ into $\mathbb R^{2d+1}$ (Menger–Nöbeling). Sandwich: $d\le m_{\min}\le2d+1$.
(iii) **[counterexample]** tripod: $m_{\min}=2>d=1$ (**C6**); equality $m_{\min}=d$ iff $\Theta$ embeds in $\mathbb R^d$.
(iv) **[counterexample]** without continuity: Borel collapse to dimension 1 (**C8**).
(v) **[impossibility]** infinite-dimensional families: no finite continuous summary (tent-cube embedding in the Lipschitz class, **C7**).
(vi) **[obstruction]** if the query-quotient of $\mathcal F$ has cardinality $>\mathfrak c$: no finite-dimensional summary of any kind (pigeonhole).

**F20 (Proposition — minimal data-computable summary). [exact value]**
Linear case: the minimal dimension of a continuous summary computable from the data and sufficient (with family info) for all identifiable queries is exactly $r=\operatorname{rank}G$. *Upper:* coordinates of $y$ in a basis of $\operatorname{col}(G)$. *Lower:* injectivity on the $r$-dimensional space $\operatorname{col}(G)$ plus invariance of domain. Existence of *any* data-computable summary is **equivalent to identifiability** — never an extra assumption ($z=\tilde y$ works whenever anything does).

---

## Part VIII — Counterexample registry (formal roles)

| Label | Object | Defeats / sharpens |
|---|---|---|
| C1 | $\mathcal F=\mathbb R^{\mathcal X}$ | F4: information exactly on $D$ |
| C2 | 1-Lipschitz class | F5: partial information without identifiability |
| C3 | $\tanh$ family | F2 $\not\Rightarrow$ F3: identifiable, unstable |
| C4 | $\sin(\theta x)$ | F14(ii): $2d+1$ attained; $k=d$ insufficient nonlinearly |
| C5 | flat $C^\infty$ family | F12 needs quasianalyticity, not smoothness |
| C6 | tripod parameter space | F19(iii): $m_{\min}>d$ possible |
| C7 | tent-cube in Lipschitz class | F19(v): no finite continuous summary |
| C8 | Borel bijection $\mathbb R^d\to\mathbb R$ | F19(iv): continuity load-bearing |
| C9 | secant witness | F15: difference identifiable, values not |
| C10 | $S^1$, even $d$ | F10(ii): parity exclusion required |
| C11 | discontinuous pullback | F10(iii): obstruction is topological |
| C12 | mirrored one-point designs | design comparison is one-directional in $\rho$ |
| C13 | forced-zero archive case | F18: dichotomy, not blanket unconstrainedness |
| C14 | two-branch family | adaptivity can help unboundedly (nonconvex); never helps (convex balanced) |
| C15 | equispaced vs Chebyshev nodes | equal-size designs differ exponentially |

---

*Every statement above carries its type tag; nothing is asserted beyond its stated hypotheses. The operator-theoretic layer (canonical mapping, its uniqueness and structure) is developed in `operator_formulation.md`; the conceptual reading in `core_principle.md`.*
