# Identifiability of a family member from $k\le 5$ evaluations

> **Provenance.** Frozen research record, 2026-08-02. Produced as a pure-mathematics derivation (no assumed representation of the family, no algorithm posited in advance). Every technically risky claim block was adversarially refereed by independent verification agents: the duality/minimax block and the RKHS/power-function block were **confirmed** outright; the trace-modulus theorem, Haar/Mairhuber–Curtis block, rank/relative-values block, dimension-count theorem, and archive block were confirmed **with qualifications**, all of which are incorporated below. The sharpness arithmetic for the $\sin(\theta x)$ example and the summary-dimension facts (Borel isomorphism, cube embedding, invariance of domain) were re-verified by hand; the minimal-dimension claim was corrected per an independent completeness audit (sandwich $d\le m_{\min}\le 2d+1$ with tripod counterexample). No claim below stands as originally drafted if a referee found a defect; this file is the corrected, post-verification text.

---

## 0. Formalization

$\mathcal X$ is a set with no assumed structure; $\mathcal F\subseteq\mathbb R^{\mathcal X}$ is a nonempty family. Three strata of data:

- **Archive** $\mathcal A$: for many $\alpha$, the restriction $f_\alpha|_{D_\alpha}$ on a finite $D_\alpha\subset\mathcal X$. Its union of design points is the *covered set* $U=\bigcup_\alpha D_\alpha$.
- **Sample**: for the new member, $\tilde y\in\mathbb R^k$ with $\max_i|\tilde y_i-f_\beta(x_i)|\le\varepsilon$ on the design $D=\{x_1,\dots,x_k\}$, $k\le 5$ (exact data is $\varepsilon=0$).
- **Query**: a point $x\notin D$; an *estimator* is any map $\Phi:\mathbb R^k\to\mathbb R$ (deterministic; randomization is addressed in Rem. 1.3).

The only object linking sample to query is the **trace set**

$$T_{D,x}\;=\;\{(f|_D,\,f(x)) : f\in\mathcal F\}\ \subseteq\ \mathbb R^{k+1},$$

and, quantitatively, the **trace modulus**

$$\omega_{x,D}(t)\;=\;\sup\Big\{\,|f(x)-g(x)| \;:\; f,g\in\mathcal F,\ \max_i|f(x_i)-g(x_i)|\le t\Big\},\qquad t\ge 0 .$$

More generally, for any functional $\varphi:\mathcal F\to\mathbb R$ set $\omega_\varphi(t)=\sup\{|\varphi(f)-\varphi(g)|:\max_i|f(x_i)-g(x_i)|\le t\}$; the case $\varphi=\delta_x$ recovers $\omega_{x,D}$. Everything below is an assertion about these objects.

---

## 1. The strongest non-identifiability statement (Question 1)

**Theorem 1 (exact minimax; no assumptions on $\mathcal F$).**
Fix $D$, $x$, $\varepsilon\ge0$. With values in $[0,\infty]$,

$$\inf_{\Phi}\ \sup_{f\in\mathcal F}\ \sup_{\tilde y\,:\,\|\tilde y-f|_D\|_\infty\le\varepsilon}\ \big|\Phi(\tilde y)-f(x)\big| \;=\; \tfrac12\,\omega_{x,D}(2\varepsilon).$$

*Proof.* **Lower bound.** If $f,g$ have $\max_i|f(x_i)-g(x_i)|\le2\varepsilon$, the coordinatewise midpoint data $\tilde y_i=\tfrac12(f(x_i)+g(x_i))$ is $\varepsilon$-consistent with both, so any $\Phi$ errs at least $\tfrac12|f(x)-g(x)|$ on one of them; take the supremum over pairs. **Upper bound.** For realizable $\tilde y$ let $S(\tilde y)=\{h(x):h\in\mathcal F,\ \max_i|h(x_i)-\tilde y_i|\le\varepsilon\}$; any two members consistent with the same $\tilde y$ have trace distance $\le2\varepsilon$ by the triangle inequality, so $\operatorname{diam}S(\tilde y)\le\omega_{x,D}(2\varepsilon)$, with equality attained in the supremum over realizable data (midpoint construction again). When $\omega_{x,D}(2\varepsilon)<\infty$, the rule $\Phi(\tilde y)=\tfrac12(\inf S+\sup S)$ on realizable data (arbitrary elsewhere — realizability is exactly $S\neq\emptyset$) has worst-case error exactly $\tfrac12\sup_{\tilde y}\operatorname{diam}S(\tilde y)$. When $\omega_{x,D}(2\varepsilon)=\infty$ the lower bound alone gives the equality. $\square$

Three remarks, each load-bearing:

**Remark 1.1 (why the constant is exact).** The prefactor $\tfrac12$ uses that the target is $\mathbb R$: in the real line, Chebyshev radius $=$ half-diameter and the midpoint is the Chebyshev center. For *joint* recovery at several queries under a joint norm, the radius of the consistent set can exceed half its diameter. However, for the sup-loss over a query set the per-query midpoint rule is again exactly optimal: $\;\inf_\Phi\sup(\text{sup-loss})=\sup_x\tfrac12\omega_{x,D}(2\varepsilon)$ (the $\ge$ direction because a joint rule induces per-query rules; the $\le$ direction by running the midpoint rule at each query).

**Remark 1.2 (the optimal rule knows $\varepsilon$).** The optimal error *modulus* is $\varepsilon\mapsto\tfrac12\omega_{x,D}(2\varepsilon)$, achieved by an $\varepsilon$-dependent rule; no single rule is claimed optimal for all $\varepsilon$ simultaneously.

**Remark 1.3 (randomization does not help).** For any random output $Z$, $\max\big(\mathbb E|Z-a|,\mathbb E|Z-b|\big)\ge\tfrac12|a-b|$, so the lower bound survives for randomized estimators in expected worst-case error.

**Corollary 1.4 (strongest non-identifiability).** For $\mathcal F=\mathbb R^{\mathcal X}$ (or any family containing, for each $M$, a pair agreeing on $D$ and differing by $\ge M$ at $x$): the minimax error is $\varepsilon$ at each $x_i\in D$ and $+\infty$ at every $x\notin D$. *The set of points at which the data carry any information at all (finite minimax error) is exactly $D$.* Membership of $f_\beta$ in a family constrains nothing unless the family itself does: complete knowledge of every other member is worthless when the trace set is a full cylinder over $\mathbb R$.

**Proposition 1.5 (tautological characterization — the exact boundary).** With exact data, $f(x)$ is determined by $f|_D$ for every $f\in\mathcal F$ **iff** $\omega_{x,D}(0)=0$ **iff** $T_{D,x}$ is the graph of a function $\Psi$ over its projection $\mathcal F|_D$. Then $\Psi$ is unique on $\mathcal F|_D$ and *every* zero-error rule agrees with $\Psi$ there. (Off $\mathcal F|_D$, i.e., on unrealizable data, no optimality principle constrains the rule; see §10.4.)

**Worked instance (partial information without identifiability).** $\mathcal F=\{f:\mathcal X\to\mathbb R,\ 1\text{-Lipschitz}\}$ on a metric space. Data $y$ is feasible iff $|y_i-y_j|\le \operatorname{dist}(x_i,x_j)$ (McShane's condition), and then

$$S(y)=\Big[\max_i\big(y_i-\operatorname{dist}(x,x_i)\big),\ \min_i\big(y_i+\operatorname{dist}(x,x_i)\big)\Big],$$

the endpoints attained by the McShane/Whitney envelope extensions, which lie in $\mathcal F$. Consequently $\omega_{x,D}(0)=2\min_i\operatorname{dist}(x,x_i)$ (extremal pair $\pm\min_i\operatorname{dist}(\cdot,x_i)$), and the exact-data minimax error is $\min_i\operatorname{dist}(x,x_i)$: never zero off $D$, but finite — a family can be informative without being identifying. (This is the unbounded Lipschitz class; adding a sup-norm bound $M$ changes $\omega$ by an $M$-dependent truncation.)

---

## 2. Minimal assumptions for inference (Question 2)

The assumptions form a strict hierarchy; each level is *exactly* equivalent to a grade of inference, which is the strongest possible form of minimality.

**A0 (identifiability): $\omega_{x,D}(0)=0$.** Necessary and sufficient for exact-data inference at $(D,x)$ (Prop. 1.5). Minimality is not rhetorical: any assumption whatsoever that enables inference *implies* A0, because A0 *is* the statement that inference is possible. Its content: identifiability is a property of the finite-dimensional set $T_{D,x}$ alone — no topology, metric, or representation of $\mathcal F$ enters.

**A1 (stability): $\omega_{x,D}(0^+)=0$.** By Theorem 1, worst-case error $\to0$ as $\varepsilon\to0$ **iff** $\omega(t)\to0$ as $t\downarrow0$, and the optimal modulus is exactly $\tfrac12\omega(2\varepsilon)$; Lipschitz inference with constant $L$ is possible iff $\omega(t)\le 2Lt$. A1 is strictly stronger than A0:

**Counterexample 2.1 (identifiable, not stable).** $\mathcal X=\{x_1,x\}$, $\mathcal F=\{f_\theta:\theta\in\mathbb R\}$ with $f_\theta(x_1)=\tanh\theta$, $f_\theta(x)=\theta$. The trace map is injective, so $\omega(0)=0$ and exact data determine $f_\beta(x)$; but $\omega(t)=\infty$ for every $t>0$, so *any* noise destroys everything. Identifiability without stability is genuinely possible.

**A2 (exact finite-sample recovery, linear category).** *Proposition:* a linear family $\mathcal F=V$ that satisfies A0 at some $k$-point design for all queries must have $\dim V\le k$, because the restriction map $V\to\mathbb R^D$ is then injective linear. Conversely $\dim V=d\le k$ plus a rank-$d$ design gives zero-error recovery (§5). So in the linear category, **finite dimension is exactly the price of zero-error recovery** — but it is *not* necessary for stable approximate inference: the unit ball of an infinite-dimensional RKHS satisfies A1 at every design (its modulus is computed by the power function, §6) while $\omega(0)>0$ at generic queries. The correct reading: A2 buys exactness; A1 alone buys convergence.

**Design-independence is a topological luxury.** On an interval, $d$-dimensional *Haar spaces* (every nonzero element has $\le d-1$ zeros) make **every** set of $d$ distinct points unisolvent — for a $d$-dimensional space $V$ of functions on *any* set with $\ge d$ points, "every nonzero $v$ has $\le d-1$ zeros" $\iff$ "every $d$-point collocation matrix is nonsingular" (pure linear algebra: a null vector of the matrix is a coefficient vector of a member vanishing at all $d$ points). Polynomials of degree $<d$ and $\mathrm{span}\{e^{s_1t},\dots,e^{s_dt}\}$ (distinct real $s_j$) qualify. But:

**Theorem 2.2 (Mairhuber–Curtis–Sieklucki).** For compact Hausdorff $\mathcal X$ with $|\mathcal X|\ge d$ and $d\ge2$: $C(\mathcal X)$ contains a $d$-dimensional Haar subspace **iff** $\mathcal X$ is homeomorphic to a closed subset of the circle $S^1$ — and if $\mathcal X\cong S^1$ itself, additionally $d$ must be odd (slide the $d$ nodes cyclically: the collocation determinant is continuous, never zero, and returns multiplied by the cyclic sign $(-1)^{d-1}$, forcing $d$ odd; trigonometric polynomials realize every odd $d$).

**Proposition 2.3 (triod obstruction, direct proof — no compactness).** If $\mathcal X$ contains a triod (three arcs meeting at a point — e.g., any open subset of $\mathbb R^n$, $n\ge2$, or any manifold of dimension $\ge2$), then for every $d\ge2$-dimensional $V\subset C(\mathcal X)$ some $d$ distinct points make evaluation singular. *Proof:* fix $x_3,\dots,x_d$ distinct outside two arms; continuously exchange $x_1,x_2$ by routing one through the third arm, keeping all points distinct; the collocation determinant is continuous along the path and ends with its sign flipped (row transposition), so it vanishes at some intermediate distinct configuration. $\square$

So in multivariate continuous settings, identifiability is unavoidably a **joint** property of family and design: the rank conditions of §5 must be checked at the actual $D$. Two scoping caveats keep this honest: (i) the obstruction is topological, not set-theoretic — pulling polynomials back through a discontinuous bijection $[0,1]^2\to[0,1]$ produces a $d$-dimensional family of (discontinuous) functions on the square for which *every* $d$-point design is unisolvent; on a bare set, design-independent unisolvence exists for every $d$ (inject $\mathcal X\hookrightarrow\mathbb R$, use Vandermonde); (ii) $d=1$ (constants) is exempt on any $\mathcal X$.

---

## 3. The information ledger (Question 3)

**General family (no structure).** Let $R_D:\mathcal F\to\mathbb R^k$ be the trace map.

- **Common to the family:** the trace sets $T_{S}$ for finite $S\subset\mathcal X$ — equivalently the family viewed through all finite windows. This is the *only* family-level information that exists, and (§10) the only part an archive can reveal is $T_S$ for $S$ inside the covered set $U$.
- **Specific to a member:** which fiber of $R_D$ it occupies, i.e., its equivalence class under "equal on $D$" — refined, as $D$ grows within $U$, toward the member itself.
- **Identifiable from $k$ evaluations (exact data):** *exactly* the functionals $\varphi:\mathcal F\to\mathbb R$ that factor through $R_D$ (constant on fibers). This class is closed under arbitrary composition — any function of identifiable quantities is identifiable. With noise, the graded version: $\varphi$ is stably identifiable iff $\omega_\varphi(0^+)=0$, and its optimal error is $\tfrac12\omega_\varphi(2\varepsilon)$ (Theorem 1 verbatim with $\varphi$ in place of $\delta_x$).
- **Unidentifiable:** the fiber partition itself — for values: the sections $\{t:(y,t)\in T_{D,x}\}$, of diameter $\omega_{x,D}(0)$ in the worst case.

**Linear family $V=\operatorname{span}\{\phi_1,\dots,\phi_d\}$** (basis; $G\in\mathbb R^{k\times d}$, $G_{ij}=\phi_j(x_i)$; $N_D=\{v\in V: v|_D=0\}$). The exact sequence $0\to N_D\to V\to V|_D\to 0$ dualizes to: the identifiable *linear* functionals of the coefficient vector are **exactly** $\operatorname{row}(G)=\operatorname{Ann}(N_D)$ — an $r$-dimensional space, $r=\operatorname{rank}G\le k$ ("exactly" scoped to linear functionals at $\varepsilon=0$; nonlinear functions of them are identifiable too, by the general ledger). Member-specific content: $c_\beta\in\mathbb R^d$; identified part: $P_{\operatorname{row}(G)}c_\beta$ ($r$ numbers); unidentifiable remainder: the coset $c_\beta+\ker G$, dimension $d-r$, which pollutes exactly the queries with $\phi(x)\notin\operatorname{row}(G)$.

---

## 4. The finite-dimensional object $z_\beta$ (Question 4)

Existence must be split by what "$z_\beta$" is allowed to be. Nothing below assumes existence.

**(a) Cardinality obstruction (the honest first check).** If the quotient of $\mathcal F$ by "equal at all required queries" has cardinality $>\mathfrak c$, **no** finite-dimensional $z_\beta$ exists at all — not even a non-measurable one — by pigeonhole into $\mathbb R^m$. Since $\mathcal X$ is arbitrary, this branch is live (e.g. $\mathcal F=\mathbb R^{\mathcal X}$ with $|\mathcal X|>\mathfrak c$; already $|\mathbb R^{\mathbb R}|=2^{\mathfrak c}$).

**(b) Finitely many required queries $x^{(1)},\dots,x^{(q)}$:** $z_\beta=(f_\beta(x^{(j)}))_j\in\mathbb R^q$ trivially suffices for *any* family. Bare existence is therefore never the question; the question is whether $z_\beta$ is **computable from the $k$ observations** — and that holds **iff** the problem is identifiable (A0 at each query), in which case $z_\beta=\tilde y\in\mathbb R^k$ itself works. Existence of a data-computable $z_\beta$ is *equivalent to*, not an extra assumption beyond, identifiability.

**(c) Sharp minimal dimension of the data-computable summary (linear case).** $m_{\min}=r=\operatorname{rank}G$. *Upper:* coordinates of $y$ in a basis of $\operatorname{col}(G)$. *Lower:* exact data ranges over the $r$-dimensional subspace $\operatorname{col}(G)$; a continuous $z=\zeta(y)$ that (with family info) determines all identifiable values must be injective on $\operatorname{col}(G)$, so $m\ge r$ by invariance of domain.

**(d) All queries, continuous category.** Here the answer is a sandwich, not an equality:

**Theorem 4.1.** Let $\mathcal F$ be injectively and bicontinuously parametrized by a space $\Theta$ (pointwise-convergence topology on $\mathcal F$).
(i) If $\Theta$ contains an open subset of $\mathbb R^d$, any continuous injective $z:\mathcal F\to\mathbb R^m$ forces $m\ge d$ (invariance of domain: compose a chart with $\mathbb R^m\hookrightarrow\mathbb R^d$; the image would be open yet lie in an $m$-flat).
(ii) If $\Theta$ is compact metric of covering dimension $d$, some continuous injective $z$ into $\mathbb R^{2d+1}$ exists (Menger–Nöbeling embedding). So $d\le m_{\min}\le 2d+1$.
(iii) $m_{\min}=d$ is **not** guaranteed: for $\Theta=$ a tripod (dimension 1), no continuous injection into $\mathbb R$ exists (the image would be an interval; deleting the branch point's image leaves $\le2$ components, but the tripod minus its branch point has 3), so $m_{\min}=2>d=1$. Equality holds exactly when $\Theta$ embeds in $\mathbb R^d$ — e.g., open parameter domains and linear families.
(iv) **Continuity is load-bearing:** any two uncountable Polish spaces are Borel isomorphic, so a *measurable* exact summary of dimension 1 always exists when the parameter space is standard Borel. "Minimal dimension" is meaningful only in the continuous/stable category.
(v) For genuinely infinite-dimensional families no finite $m$ works: the 1-Lipschitz class on $[0,1]$ contains homeomorphic copies of $[0,\delta]^n$ for every $n$ ($n$ disjoint tents of height $t_i\le\delta$ and slope $\le1$; the heights are read off at the peaks, and the sup-distance equals $\max_i|t_i-t_i'|$), so a continuous injection into $\mathbb R^m$ would embed an $n$-cube, $n>m$ — impossible.

**Conclusion for Question 4:** the finite-dimensional object exists *and is computable from the data* precisely when the problem is identifiable; its minimal dimension is $\operatorname{rank}G$ (linear case) and is sandwiched by $[\dim,\,2\dim+1]$ in general; asking for it without identifiability, or for an infinite-dimensional family at all queries, is provably vacuous.

---

## 5. Rank, uniqueness, stability, sensitivity (Question 5)

Standing setting: $V$, $G$, $\phi(x)$ as above; exact statements verified.

**Theorem 5.1 (rank).** For $\mathcal F=V$: $\lambda^\top c$ is determined by $y=Gc$ (uniformly in $c$) **iff** $\lambda\in\operatorname{row}(G)$, and then equals $\lambda^\top G^+y$ for every consistent $c$. In particular $f_\beta(x)$ is identifiable iff $\phi(x)\in\operatorname{row}(G)$; all of $c_\beta$ iff $\operatorname{rank}G=d$ (forcing $k\ge d$). For a constrained coefficient set $C\subsetneq\mathbb R^d$: identifiability on $C$ **iff** $\ker G\cap(C-C)=\{0\}$ — so nonlinear or constrained families can be identifiable below full rank.

**Theorem 5.2 (exact minimax and optimal weights; hypotheses: $\mathcal F=V$ unconstrained, known $V$, closed noise ball in an arbitrary norm $\|\cdot\|$ on $\mathbb R^k$, scalar linear target).** Assume $\phi(x)\in\operatorname{row}(G)$ and define the **generalized Lebesgue function**

$$\Lambda_*(x)\;=\;\min\{\|w\|_*\;:\;G^\top w=\phi(x)\}\quad(\text{dual norm; minimum attained}).$$

Then the minimax error equals $\varepsilon\,\Lambda_*(x)$ **exactly**, attained by the *linear* rule $\tilde y\mapsto \hat w^\top\tilde y$ at any minimizer $\hat w$ (whose worst-case error is exactly $\varepsilon\|\hat w\|_*$). *Proof core:* upper bound by duality of the pairing; lower bound by the two-point argument plus the Hahn–Banach identity $\sup\{w_0^\top v:v\in\operatorname{col}(G),\|v\|\le1\}=\min\{\|w\|_*:G^\top w=G^\top w_0\}$ (norm of a functional on a subspace $=$ min over extensions; the extension realizing the min exists in finite dimension without choice). If $\phi(x)\notin\operatorname{row}(G)$, both sides are $+\infty$. Linearity costing nothing is Smolyak's theorem (1965) for exact information; the noisy-values case is Marchuk–Osipenko (1975); survey: Micchelli–Rivlin (1977). (The result fails for jittered *nodes* as opposed to jittered values — a nearby trap, Kacewicz–Plaskota 2003.)

**Corollary 5.3 (Euclidean case, sensitivity).** $\hat w=(G^+)^\top\phi(x)$ for any rank (given feasibility), and $\Lambda_2(x)=\|(G^+)^\top\phi(x)\|_2\le\|\phi(x)\|_2/\sigma_{\min}(G)$. The weight vector **is** the sensitivity profile: $\partial\hat f(x)/\partial\tilde y_i=\hat w_i$. If $\mathbf 1\in V$ then $\sum_i w_i=1$ for *every* feasible $w$ (pair $G^\top w=\phi(x)$ with the coefficient vector of $\mathbf 1$).

**Theorem 5.4 (nonlinear local rank).** For a $C^1$ family $\theta\mapsto f_\theta$, $\theta\in\Theta\subseteq\mathbb R^d$ open, with Jacobian $J(\theta)=\big(\partial_\theta f_\theta(x_i)\big)\in\mathbb R^{k\times d}$: if $\operatorname{rank}J(\theta_0)=d$, the evaluation map is locally injective near $\theta_0$ with local stability $\|\delta\theta\|\lesssim\|J(\theta_0)^+\|\,\|\delta y\|$; if $\operatorname{rank}J\equiv r<d$ on a neighborhood, the constant-rank theorem makes the fibers $(d-r)$-manifolds — locally non-identifiable, with the tangent space of the ambiguity equal to $\ker J$.

**Uniqueness.** Unique recovery of the member $=$ triviality of the ambiguity set: $\ker G\cap(C-C)=\{0\}$ (linear/constrained), injectivity of the evaluation map (nonlinear; §8 for when this is generic).

**Stability under archive error.** If $V$ is known only up to a subspace estimated from the archive with principal-angle error $\delta$ (Davis–Kahan/Wedin regime), the recovery error acquires an additive model term: $|\hat f(x)-f_\beta(x)|\le\varepsilon\Lambda_*(x)+\delta\cdot C(x,G)\|c_\beta\|$, with $C(x,G)$ controlled by $\sigma_{\min}(G)$ and $\|\phi(x)\|$ — degrading gracefully in $\delta$ but *blowing up as the design approaches rank deficiency*, which is the correct qualitative warning: ill-conditioned designs amplify both noise and model error.

---

## 6. When two designs carry different information, and the quantity that measures it (Questions 6–7)

**Derivation, not decree.** All information about $f_\beta(x)$ resides in the section $S(\tilde y)$ of the trace set (§1); a design's quality at $x$ is the size of that section in the worst case; Theorem 1 says the resulting scalar is *exactly* the optimal error. So the derived quantity is

$$\rho(x;D,\varepsilon)\;=\;\tfrac12\,\omega_{x,D}(2\varepsilon)\qquad(\text{the radius of information at }x).$$

It is monotone under design refinement ($D\subseteq D'$ shrinks every section), needs no metric or topology to define, and Theorem 1 is its certificate of exactness. Its closed forms under the three structural regimes:

| family | $\rho(x;D,\cdot)$ | regime |
|---|---|---|
| 1-Lipschitz class | $\min_i\operatorname{dist}(x,x_i)$ $(+\,\varepsilon)$ | geometric coverage |
| linear $V$, $\dim d$ | $\varepsilon\,\Lambda_*(x)$ | conditioning / Lebesgue function |
| RKHS unit ball, exact data | $P_D(x)\sqrt{1-\|s_y\|_{\mathcal H}^2}$ | power function |

where $P_D(x)^2=K(x,x)-k_x^\top K_{DD}^{-1}k_x$ satisfies the three verified identities: distance from $K(\cdot,x)$ to $\operatorname{span}\{K(\cdot,x_i)\}$ $=$ worst-case interpolation error over the unit ball $=$ Gaussian conditional standard deviation (an identity of formulas across genuinely different models — GP paths a.s. lie outside $\mathcal H$ when $\dim\mathcal H=\infty$, Driscoll), and the exact-data consistent set is precisely the Golomb–Weinberger interval $s_y(x)\pm P_D(x)\sqrt{1-\|s_y\|^2}$.

**Two designs differ at $x$ exactly in one direction:** $\rho(x;D_1,\varepsilon)\ne\rho(x;D_2,\varepsilon)$ implies they carry different information, and the difference can be enormous — for polynomial interpolation with $k$ nodes on $[-1,1]$, $\Lambda\sim\frac2\pi\log k$ at Chebyshev nodes versus $\Lambda\sim 2^{k}/(e\,k\log k)$ at equispaced nodes (Turetskii 1940 / Schönhage 1961; Luttmann–Rivlin for Chebyshev). The converse is **false**: $\rho$ is a scalar compression, and designs with equal radii need not be informationally equivalent (Lipschitz class on $\mathbb R$, $D_1=\{-1\}$, $D_2=\{+1\}$, query $0$: equal $\rho$, different section maps). The complete comparison is the Blackwell-type preorder "every $D_1$-section is contained in a $D_2$-section"; $\rho$ is its canonical scalar quotient.

**Adaptivity.** For convex, balanced $\mathcal F$ and a linear-functional target, adaptive selection of the $x_i$ gains nothing over the best nonadaptive design (Bakhvalov; Gal–Micchelli). For nonconvex families it can gain unboundedly: with $\mathcal F=\mathcal F_A\cup\mathcal F_B$ supported on disjoint regions, one adaptive evaluation identifies the branch and spends the remaining $k-1$ points inside it, while any nonadaptive design must split its budget. With $k\le5$ this dichotomy is decision-relevant.

---

## 7. Relative values (Question 8)

The general principle is the ledger of §3: a functional is identifiable iff constant on data-fibers; for an affine functional, iff its linear part annihilates the ambiguity set. Applied to $\ell=\delta_{x_a}-\delta_{x_b}$:

**Theorem 7.1 (linear case).** $f_\beta(x_a)-f_\beta(x_b)$ is identifiable iff $\phi(x_a)-\phi(x_b)\in\operatorname{row}(G)$ — which can hold when neither $\phi(x_a)$ nor $\phi(x_b)$ lies in $\operatorname{row}(G)$. Witness ($V=\operatorname{span}\{1,h_1,h_2\}$, $k=2$, $H=(h_1,h_2)$): $H(x_1)=(0,0)$, $H(x_2)=(1,0)$, $H(x_a)=(0,1)$, $H(x_b)=(1,1)$ — the difference is identifiable, neither value is. The exact geometric criterion: $H(x_a)-H(x_b)\in\operatorname{span}\{H(x_1)-H(x_2)\}$ (scalar multiple, zero allowed; if $H(x_1)=H(x_2)$ only the trivial difference is identifiable). The correct strictness statement: *joint identifiability of both values is strictly stronger than identifiability of their difference; and given the difference is identifiable, $\phi(x_a)\in\operatorname{row}(G)\iff\phi(x_b)\in\operatorname{row}(G)$.* Quantitatively, $\omega_{\delta_a-\delta_b}$ can vanish while $\omega_{\delta_a}=\omega_{\delta_b}=\infty$.

**Theorem 7.2 (gauge taxonomy — when only a transformation of a common $g_0$ is determined).** With the convention that a partially defined functional is identifiable iff constant on (fiber $\cap$ domain):
- $\mathcal F=g_0+\mathbb R\mathbf 1$ (additive member effect): **all pairwise differences are identifiable with $k=0$** — they are family-level information — while no absolute value is; $k=1$ identifies everything. (Honesty clause from §9: at $k=0$ this holds at archive-covered pairs, since $g_0$ itself is known only there.)
- $\mathcal F=\{a\,g_0\}$: all ratios identifiable at $k=0$ where defined (whether the ratio is defined at the true member — $a\ne0$ — is itself unidentifiable at $k=0$); one evaluation with $g_0(x_1)\ne0$ identifies everything.
- $\mathcal F=\{a\,g_0+b\}$ with $\{1,g_0\}$ independent: at $k=0$ exactly the affine invariants (e.g. ratios of differences) are identifiable; at $k=1$ a one-parameter gauge mixing scale and offset survives, and $f(x_a)-f(x_b)$ is identifiable **iff** $g_0(x_a)=g_0(x_b)$ (value $0$); at $k=2$ with $g_0(x_1)\ne g_0(x_2)$, everything. (If $\{1,g_0\}$ is dependent the family collapses to constants and $k=0$ already identifies all differences as $0$ — parameter non-identifiability without member non-identifiability.)

The gauge picture: the residual ambiguity after $k$ observations is the stabilizer of the observed trace inside the family's symmetry; identifiable functionals are its invariants. Relative values survive precisely the additive gauge.

---

## 8. The weakest property that flips the problem, quantified (Question 9)

**The weakest property is A0** — $T_{D,x}$ is a graph over its projection — by Prop. 1.5 it is implied by every assumption that enables inference. It needs no structure on $\mathcal X$ whatsoever. Everything else in this treatise is the *computation* of A0/A1 under added structure. The decisive quantitative theorem in the smooth regime:

**Theorem 8.1 (generic global identifiability).** Let $\Theta\subset\mathbb R^d$, $X\subset\mathbb R^m$ be open and bounded, $X$ connected, $F$ real-analytic on a neighborhood of $\overline\Theta\times\overline X$, and the parametrization separating ($\theta\ne\theta'\Rightarrow f_\theta\ne f_{\theta'}$; on each component, if $X$ is disconnected). Then for $k\ge 2d+1$ there is a **closed** subanalytic set $N\subset X^k$ of dimension $\le 2d+k(m-1)<km$ (hence Lebesgue-null) such that every design in the open, dense, full-measure set $X^k\setminus N$ makes $\theta\mapsto(f_\theta(x_i))_i$ **injective on all of $\Theta$**.
*Proof skeleton:* for $\theta\ne\theta'$ the agreement set is a proper analytic subset of connected $X$ (identity theorem), of dimension $\le m-1$ and measure zero; the incidence set over pairs is relatively compact semianalytic with fibers of dimension $\le k(m-1)$, so its dimension is $\le 2d+k(m-1)$; projections of bounded subanalytic sets do not raise dimension, and the o-minimal frontier theorem upgrades the bad set to its closure with the same bound. Boundedness is load-bearing (unbounded subanalytic sets lose projection stability); connectedness is not removable.

**Theorem 8.2 (stable version).** If additionally the family is infinitesimally separating (no $\theta$, $u\ne0$ with $\partial_uf_\theta\equiv0$), then for $k\ge 2d+1$ generic designs make the evaluation map simultaneously injective and an immersion; on every compact convex $K\subset\Theta$ there is $c(K)>0$ with $\|\Phi(\theta)-\Phi(\theta')\|\ge c\|\theta-\theta'\|$, so $\varepsilon$-consistent estimates satisfy $\|\hat\theta-\theta\|\le 2\varepsilon\sqrt k/c(K)$ (constants degrade as $K\uparrow\Theta$). Immersion alone is generic already at $k\ge2d$, but the combined guarantee carries the $2d+1$ threshold.

**Theorem 8.3 (necessity).** If the evaluation map is injective (exact recovery) and continuous on a family containing a $d$-cell, then $k\ge d$ (invariance of domain). With $k\le5$: exact stable recovery confines the family to $d\le5$; the generic-global guarantee of Thm 8.1 is available for $d\le2$.

**Sharpness (hand-verified arithmetic).** $f_\theta(x)=\sin(\theta x)$, $\Theta=\mathbb R$, $d=1$, separating ($\frac{d}{dx}\sin(\theta x)|_0=\theta$):
- $k=1$ fails for **every** design ($\theta\mapsto\sin(\theta x_1)$ is periodic, or constant at $x_1=0$).
- $k=2$ fails for **every** design $\{a,b\}$, $ab \ne 0$: the pair $\theta=\frac{\pi}{2b}-\frac{\pi}{a}$, $\theta'=\frac{\pi}{2b}+\frac{\pi}{a}$ collides — $\theta'a=\theta a+2\pi$, and $\theta' b=\frac\pi2+\frac{\pi b}{a}$, $\theta b=\frac\pi2-\frac{\pi b}{a}$ with $\sin(\frac\pi2\pm u)=\cos u$.
- $k=3=2d+1$ succeeds for all designs with pairwise irrational ratios: a collision forces each $x_i$ into a "translation" case ($(\theta'-\theta)x_i\in2\pi\mathbb Z$) or a "reflection" case ($(\theta'+\theta)x_i\in\pi+2\pi\mathbb Z$); by pigeonhole two indices share a case, and either case forces a rational ratio (the reflection case cannot have $\theta+\theta'=0$ since $0\notin\pi+2\pi\mathbb Z$).

So $2d+1$ is attained: the gap between linear families (where $k=d$ suffices — unisolvent designs always exist by the independence lemma, and generically for analytic $\phi_j$ on connected domains since $\det G$ is analytic and $\not\equiv0$) and nonlinear ones is real. (For unbounded $\Theta$ the good set here is full-measure but not open — rational ratios are dense — consistent with the boundedness hypothesis of Thm 8.1.)

**The dividing line is quasianalyticity, not smoothness.** $f_\theta(x)=\theta\,e^{-1/x}\mathbf 1_{x>0}$ on $X=(-1,1)$ is $C^\infty$, definable in $\mathbb R_{\exp}$, separating — yet every agreement set contains $(-1,0]$, so bad designs have *positive measure and nonempty interior for every $k$*: no genericity statement of any kind holds. The exact minimal hypotheses for Thm 8.1: (a) definability in some o-minimal expansion of $\mathbb R$ (for the dimension theory) **plus** (b) the identity-theorem property for differences (empty-interior agreement sets). Quasianalytic Denjoy–Carleman families satisfy both (Rolin–Speissegger–Wilkie), so the theorem extends to them; mere $C^\infty$ satisfies neither usefully (by Whitney, agreement sets of smooth families are arbitrary closed sets).

---

## 9. The archive: what the family-level information is, and its limits

**Theorem 9.1 (identification of $V$ from an incomplete archive — exact iff, within the class of $d$-dimensional spaces, $d$ known).** Members $g_1,\dots,g_n$ of an unknown $d$-dimensional $V$, member $j$ observed on $D_j$; suppose a common core $X_0\subseteq\bigcap_jD_j$, $|X_0|=d$, with the $n\times d$ matrix $A=(g_j|_{X_0})$ of rank $d$ (data-checkable). Then the restriction $V\to\mathbb R^{X_0}$ is bijective ($X_0$ is unisolvent — surjectivity from the rows of $A$, bijectivity from $\dim V=d$) and the $g_j$ span $V$. For $x\in U$: the representing row $a_x$ with $g(x)=a_x\cdot g|_{X_0}$ for all $g\in V$ is determined **iff** the submatrix of $A$ over members observed at $x$ has rank $d$; otherwise the ambiguity is an affine subspace of dimension $d-\operatorname{rank}$, *every point of which is realized by a genuine candidate family* — so the rank condition is necessary, not just sufficient. When it holds at every $x\in U$, $V|_U$ is fully identified, and the new member's inference at queries in $U$ proceeds by Theorem 5.1 with rows $a_{x_i}$. (The common core is one sufficient route; chaining through pairwise design overlaps can substitute for it.)

**Theorem 9.2 (impossibility off the covered set — the recursion, exact dichotomy).** For a query $x\notin U\cup D$: over all $d$-dimensional families consistent with archive and sample, the set of possible values $f_\beta(x)$ is **either all of $\mathbb R$, or $\{0\}$** — the latter precisely in the degenerate case where the data force $f_\beta=0$ (e.g. $y=0$ under a full-rank design; with $\varepsilon>0$ this exception disappears and the answer is always all of $\mathbb R$). *Construction:* alter a basis of $V$ at the single point $x$ arbitrarily; restrictions to $U$ are unchanged, so the archive and sample are matched while $f_\beta(x)=c_\beta\cdot t$ sweeps $\mathbb R$ whenever $c_\beta\ne0$.
**Meaning:** the original non-identifiability theorem (Cor. 1.4) recurs one level up. Family-level information at a query exists only where the archive covers it, or where an assumption on the *members* (continuity toward $\overline U$, analyticity, RKHS membership — the last yielding interval bounds, not values) transports it. No amount of archive data about *other* points substitutes.

---

## 10. The procedure, derived rather than posited (Question 10)

The derivation is now forced, and its premises are exactly the theorems above — listed, because "forced" is honest only relative to them:

1. **Archive $\to$ trace sets** on the covered set (premise: Thm 9.1's rank conditions; scope certificate: Thm 9.2 — do not promise queries off $U$ without a member-level assumption; stability certificate: the Davis–Kahan remark of §5).
2. **Sample $\to$ consistent section** $S(\tilde y)$ of the trace set (premise: the noise model and its norm; under the model, realizability is automatic).
3. **Output $\to$ the midpoint (Chebyshev center) of $S(\tilde y)$** — forced by Theorem 1: any other output has strictly larger worst-case error at some consistent data.
4. **Off the consistent tube** (misspecified data at distance $\eta$ from realizability): project onto the realizable set and run the rule at inflated tolerance; the error certificate degrades from $\tfrac12\omega(2\varepsilon)$ to $\tfrac12\omega(2\varepsilon+2\eta)$ plus the family's approximation error at the query — the doubling is the price of a convention where optimality is silent.

The classical formulas now *drop out as closed forms of step 3*, none of them assumed: for linear $V$, the midpoint is the min-dual-norm weighted estimate $\hat w^\top\tilde y$ of Theorem 5.2 (exactly optimal, linear for free); for the RKHS ball, it is the kernel interpolant $s_y(x)$, center of the Golomb–Weinberger interval; for the Lipschitz class, the average of the two envelope extensions. The algorithmic zoo is the shadow of one theorem applied to three trace geometries.

---

## Certificate ledger

| Construction | Certificate delivered |
|---|---|
| Trace modulus $\omega$ | exact minimax identity (Thm 1); lower bound survives randomization |
| Non-identifiability | Cor. 1.4: information exactly on $D$; $\infty$ elsewhere |
| A0 graph property | necessary *and* sufficient (Prop. 1.5); uniqueness of $\Psi$ on $\mathcal F|_D$ |
| A0 vs A1 | separating counterexample ($\tanh$ family, Ex. 2.1) |
| A2 finite dimension | iff-price of zero-error recovery (linear); RKHS ball shows not necessary for stability |
| Design-independence | Haar iff-lemma (any set); Mairhuber–Curtis–Sieklucki with parity exclusion; triod swap proof; bare-set escape |
| Rank theorem | Thm 5.1 (iff, both constrained and unconstrained) |
| Stability | Thm 5.2: exact minimax $=\varepsilon\Lambda_*(x)$, linear-optimal (Hahn–Banach); $\sigma_{\min}$ bound; weights-sum-to-one |
| Sensitivity | weight vector $=$ derivative of estimate; archive-error term via subspace angle |
| $z_\beta$ existence | cardinality obstruction; equivalence with identifiability; $m_{\min}=\operatorname{rank}G$; sandwich $d\le m_{\min}\le2d+1$ with tripod counterexample; Borel collapse without continuity; Lipschitz cube-embedding impossibility |
| Design comparison | $\rho$ exactness (Thm 1); one-directional implication with equal-$\rho$ counterexample; node asymptotics $\frac2\pi\log k$ vs $2^k/(ek\log k)$; adaptivity dichotomy |
| Relative values | annihilator iff (Thm 7.1) with explicit witness and degenerate cases; gauge taxonomy with exact observation counts |
| Weakest property | minimality of A0 (tautologically necessary); Thm 8.1 with closed-null exceptional set; necessity $k\ge d$; sharpness $2d+1$ ($\sin$ family, both failures explicit); quasianalyticity dividing line (flat $C^\infty$ counterexample) |
| Archive step | Thm 9.1 exact iff within fixed-$d$ class; Thm 9.2 dichotomy ($\mathbb R$ or forced $\{0\}$) — the recursion |
| Procedure | forced by Thm 1 given listed premises; misspecification doubling bound |

**The single-sentence resolution:** with no assumptions the data determine nothing off the design (and the minimax theorem says exactly how much "nothing" costs); the weakest property that changes this is that the family's trace set over design-plus-query be the graph of a function, the optimal error is identically half the trace modulus at twice the noise level, and every classical inference formula — interpolation weights, kriging, envelope averaging — is the closed form of one derived object: the Chebyshev center of the trace-set section, computed in a trace geometry the archive can only ever reveal on the points it covers.
