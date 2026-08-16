# The Core Principle

> **Status:** theory-refinement layer, 2026-08-02. Source of record: `../00_raw_outputs/identifiability_treatise_raw.md`; formal statements: `theorem_formalization.md` (F-numbers); operator layer: `operator_formulation.md` (OP-numbers). All new propositions in this file (CP-1 … CP-3) were adversarially refereed; qualifications are incorporated.

Setting as in F-Part 0: nonempty $\mathcal F\subseteq\mathbb R^{\mathcal X}$ on an arbitrary set; design $D=\{x_1,\dots,x_k\}$, $k\le5$; data $\tilde y$ with $\|\tilde y-f_\beta|_D\|_\infty\le\varepsilon$; query $x\notin D$. Trace set $T=T_{D,x}$; section $S_\varepsilon(\tilde y)$; modulus $\omega_{x,D}$.

---

## 1. The minimal mathematical object that must be recovered

There are two objects, one per level, and both admit exact minimality statements.

### 1.1 Family level: the trace set — nothing more, nothing less

**CP-1 (Sufficiency — refereed and strengthened).** For fixed $(D,x,\varepsilon)$, not only the minimax risk and the optimal rules but **every estimator's entire worst-case risk profile** is a functional of $T$ alone: the risk of any $\Phi$ equals $\sup_{(u,t)\in T}\sup_{\|\tilde y-u\|_\infty\le\varepsilon}|\Phi(\tilde y)-t|$, and the section map is $S_\varepsilon(\tilde y)=\{t:\exists(u,t)\in T,\ \|u-\tilde y\|_\infty\le\varepsilon\}$ — literally a function of $(T,\tilde y,\varepsilon)$.

**CP-2 (Minimality).** Conversely, $T$ is recovered from the **exact-data** section map as its graph: $T=\{(y,t):t\in S_0(y)\}$. Hence two families induce the same inference problem at $(D,x)$ *for all $\varepsilon\ge0$, including $\varepsilon=0$,* **iff** they have the same trace set; no strictly coarser invariant of the family determines the problem.
*The $\varepsilon=0$ clause is load-bearing:* families realizing $T=\mathbb Q\times\{0\}$ and $T'=\mathbb R\times\{0\}$ ($k=1$) have identical sections, realizable sets, risks, and optimal rules for **every** $\varepsilon>0$, yet $T\ne T'$. At any fixed positive noise level the trace set is identified only up to such dense-in-section modifications.

Under identifiability (F2), the trace set collapses to its graph function: the minimal family-level object *is the function* $\Psi_{D,x}:\mathcal F|_D\to\mathbb R$, unique on its domain.

### 1.2 Member level: the fiber, and its exact minimal summary

The data can never point past the fiber: what is recovered about $f_\beta$ is its equivalence class under $f\sim g\iff f|_D=g|_D$ (exact data), refined by nothing. A functional of the member is computable from exact data **iff it is constant on every fiber over realizable data** — iff it factors through the trace map (F6, general ledger).

**CP-3 (minimal member summary — query-relative; refereed with correction).** Fix a query set $Q$ and let the linear structure of F-D8 hold. The data determine $c_\beta$ exactly up to the coset $c_\beta+\ker G$. The minimal continuous data-computable summary sufficient for the identifiable queries **in $Q$** is
$$z_\beta = P_W\,c_\beta,\qquad W=\operatorname{span}\{\phi(x):x\in Q,\ \phi(x)\in\operatorname{row}(G)\},\qquad \dim z_\beta=\dim W,$$
and no continuous summary of dimension $<\dim W$ suffices (the decoding constraint forces a continuous injection of a $\dim W$-dimensional linear image of the data; invariance of domain). Only for the class of **all** identifiable linear functionals does this equal $P_{\operatorname{row}(G)}c_\beta$ with dimension exactly $r=\operatorname{rank}G$ — automatic when $Q\supseteq D$, since the rows of $G$ span $\operatorname{row}(G)$.
*Strictness of the query-relativity:* $V=\operatorname{span}\{1,t\}$ on $\mathbb R$, $D=\{0,1\}$ (so $r=2$), $Q=\{2\}$: the minimal object is the single scalar $\Psi_{D,2}(y)=2y_2-y_1$, of dimension $1<r$.

**The principle in one line:** what must be recovered is *(the trace set; the position of the data within it)* — the first from the family, the second from the observations, and nothing else exists to be recovered.

---

## 2. What the observed family provides

**The family provides the constraint, and only the constraint.** Its entire contribution is the geometry of $T$ — which pairs $(u,t)$ are jointly possible. Three exact statements:

1. **It provides nothing about the fiber.** By CP-1, the family enters the problem only through $T$; by the realization proposition (OP-1: *every* nonempty subset of $\mathbb R^{k+1}$ is the trace set of some family, using only that $x_1,\dots,x_k,x$ are distinct), no further property of the family — cardinality, symmetry, parametrization — can shrink a section beyond $T$'s geometry. The decomposition *constraint (family) $\times$ selection (data)* is exact, with no cross-term (CP-1 + OP-1, refereed).

2. **It is learnable only through finite windows on the covered set.** What an archive can reveal is at most the window sets $T_S=\{f|_S:f\in\mathcal F\}$ for finite $S\subseteq U$ (F17 gives the exact conditions in the linear class; F18 the exact limit off $U$). A useful closure fact (refereed): the class of windows absorbs queries — $T_{D\cup\{x\},x'}$ is a coordinate permutation of the joint trace over $(D,\{x,x'\})$ — so "equal trace sets on all finite windows of $U$" is equivalent to "equal finite restriction sets on $U$", and then **every** finite-data inference problem over $U$ (multi-query, multi-window, every estimator's risk) is identical between the two families.

3. **The member set itself is *not* provided — only its windows.** Calibration example (refereed): on $\mathcal X=\mathbb N$, the family of finitely-supported $0$–$1$ sequences and the family of *all* $0$–$1$ sequences have identical finite windows, hence are indistinguishable from any archive plus any finite sample — yet their member sets differ. Even at the family level, only $T$-level structure is ever identifiable; the indexing/parametrization (the "gauge") never is.

---

## 3. What the finite observations provide

**The observations provide the selection.** Their entire contribution:

1. **Fiber selection.** The data replace the $k=0$ prior section $T_x=\{f(x):f\in\mathcal F\}$ by the posterior section $S_\varepsilon(\tilde y)$ — a subset of $T_x$ cut out by the $k$ consistency constraints. The information supplied is exactly this reduction; its worst-case residual size is the radius $\rho=\tfrac12\omega_{x,D}(2\varepsilon)$ (F1), and per-functional the residual is $\tfrac12\omega_\varphi(2\varepsilon)$.

2. **Which functionals become determined.** Exactly those constant on the residual fibers: at $\varepsilon=0$, the functionals factoring through the trace map; in the linear case, the functionals in $\operatorname{row}(G)$ — $r=\operatorname{rank}G\le k$ linear dimensions of the member, never more than $k$ regardless of the family's size.

3. **Monotone refinement, non-monotone estimates.** Each added observation shrinks every section (OP-8), so the *guarantee* improves monotonically; the *estimates* themselves need not move monotonically (midpoints of nested intervals can oscillate: $[0,10]\supset[9,10]\supset[9,9.2]$ has midpoints $5,\ 9.5,\ 9.1$).

4. **Design quality is real and quantified.** Two designs of equal size can supply exponentially different information about the same query ($\rho$-comparison, F10(iv)); with $k\le5$ the placement of the five points, not their number, is frequently the binding constraint.

---

## 4. What is fundamentally impossible to recover

The complete list, each item with its certificate and its exact scope (refereed):

1. **Off-design values, absent constraints.** For any family containing pairs that agree on $D$ and differ arbitrarily at $x$ (in particular $\mathcal F=\mathbb R^{\mathcal X}$): infinite minimax error at every $x\notin D$, for every estimator, randomized included. Information exists exactly on $D$. *(Certificate: two-point argument; F4.)*

2. **The position within a section.** Given $\tilde y$, no rule distinguishes members of the same $\varepsilon$-consistent fiber; the irreducible conditional error is the section half-diameter, attained. Linear form: the coset $c_\beta+\ker G$, invisible precisely at queries with $\phi(x)\notin\operatorname{row}(G)$. *(Certificate: definition of section + F6.)*

3. **Family structure off the covered set — scoped.** Within the linear exactly-$d$-dimensional candidate class of F17/F18: at $x\notin U\cup D$ the consistent value set is all of $\mathbb R$, or $\{0\}$ in the forced-zero degenerate case. *(Certificate: F18. The dichotomy is a linear-class statement; for general families the correct general statement is item 1 applied one level up — no window of the archive constrains an uncovered point without member-level assumptions.)*

4. **Sub-resolution distinctions.** No rule distinguishes members whose traces differ by $\le2\varepsilon$; the noise floor $\tfrac12\omega(2\varepsilon)$ is irreducible. *(Certificate: F1 lower bound.)*

5. **The gauge.** Two parametrized families with equal finite windows on the covered region induce literally identical inference problems for every finite dataset (by the window-closure fact of §2.2) — identical achievable risks, identical optimal outputs. The parametrization, the indexing, and even the member set (§2.3) are unidentifiable; only $T$-level structure ever is. *(Certificate: CP-1/CP-2 + window closure; no topology enters.)*

---

## 5. The two-level principle

The archive problem is **structurally parallel** to the base problem — the same decomposition into *constraint* and *selection* recurs with the family in the role of the member: the archive supplies finite windows (selections) of an unknown constraint object, and the gauge-invariant objects at the upper level are again trace sets, not coordinates (a basis of $V$ is exactly the kind of indexing that §4.5 proves unidentifiable). The parallel is organizational, not a substitution instance of F1: the exact identifiability content at the upper level is supplied solely by the archive theorems (F17 for the positive direction, F18 for the limit). Its consequence stands: **the theory is closed under its own recursion** — the impossibility that governs members governs families, one level up, with the covered set playing the role of the design.

---

*Summary of the principle.* The problem factors exactly into a constraint and a selection. The constraint is the trace set: minimal, sufficient, gauge-free, learnable only through covered windows. The selection is the data: it buys at most $k$ dimensions of the member, exactly the functionals constant on its residual fiber. Everything impossible is impossible for one of five listed reasons, each certified. There is nothing else in the problem.
