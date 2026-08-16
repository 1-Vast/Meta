# RESEARCH HANDOFF — Few-Shot Identifiability of Function-Family Members

> **Frozen mathematical specification. 2026-08-02.**
> Purpose: transfer the mathematical foundation to a future researcher who will design computational systems. This document contains **no architectures, no application domains, no engineering choices** — only established mathematics and its open edges.
>
> **Corpus map (sources of record, in dependency order):**
> - `00_raw_outputs/identifiability_treatise_raw.md` — the full derivation with proofs (§0–§10, certificate ledger).
> - `01_math_foundation/theorem_summary.md` — frozen foundation (F-statement summaries, counterexamples C1–C15, assumption register).
> - `02_theory_refinement/` — `core_principle.md` (CP-numbers), `theorem_formalization.md` (formal statements F1–F20), `operator_formulation.md` (OP-numbers).
> - `03_meta_learning_principle/` — `meta_learning_abstraction.md` (MP-numbers, predictions P1–P10), `operator_specification.md` (interface, axioms A1–A10), `failure_modes.md` (certified catalogue).
>
> **Provenance:** every substantive claim was adversarially refereed by independent verification passes (three rounds across the phases); refuted or overclaimed drafts were corrected before freezing. Numbered references below (F, C, CP, OP, MP, A, P) point into the corpus.
>
> **Epistemic status:** the corpus contains **no conjectures**. Every statement below is proven under its stated hypotheses, with one flagged exception (the archive-noise quantitative constant, marked OPEN — §9). All guarantees are **worst-case (minimax)**; the theory is distribution-free and contains no average-case statements.

---

## 1. Research objective

**The abstract problem.** Let $\mathcal F\subseteq\mathbb R^{\mathcal X}$ be a nonempty family of related real-valued functions on a set $\mathcal X$ carrying no assumed structure. Three strata of data:

- **Archive:** finitely many previous members $f_\alpha$, each observed exactly on its own finite design $D_\alpha$; covered set $U=\bigcup_\alpha D_\alpha$.
- **New-member sample:** for a new member $f_\beta$, values $\tilde y\in\mathbb R^k$ with $\max_i|\tilde y_i-f_\beta(x_i)|\le\varepsilon$ on a design $D=\{x_1,\dots,x_k\}$, $k\le5$, noise level $\varepsilon$ known (closed ball, $\varepsilon=0$ = exact data).
- **Query:** the value $f_\beta(x)$ is required at points $x\notin D$ — *query-dependent inference*: which queries are answerable, how accurately, and with what sensitivity all depend on $x$.

**The demand placed on the theory** (met by the corpus): weakest assumptions, fewest free quantities, and a certificate — necessary condition, sufficient condition, rank theorem, stability bound, lower bound, or counterexample — for every construction. No algorithm was posited in advance; the inference procedure was *derived* from an optimality principle.

---

## 2. Mathematical impossibility results

**The master identity (F1).** With no assumptions beyond nonemptiness of $\mathcal F$, define the trace modulus
$$\omega_{x,D}(t)=\sup\{|f(x)-g(x)|:f,g\in\mathcal F,\ \max_i|f(x_i)-g(x_i)|\le t\}.$$
Then the minimax error of *any* estimator (randomization included) is **exactly** $\tfrac12\,\omega_{x,D}(2\varepsilon)$ (scalar target; multi-query evaluation is exact only under sup-loss — F1 Rem. 1.1/1.2). Everything below is an instance.

**What cannot be recovered.**
1. **Off-design values, absent constraints (F4/C1).** For families containing pairs agreeing on $D$ with arbitrarily different query values (e.g. $\mathcal F=\mathbb R^{\mathcal X}$), the minimax error is $\varepsilon$ at each design point and $+\infty$ everywhere else. *Finite information exists exactly on $D$.* Complete knowledge of every other member is worthless unless the family constrains the joint trace set.
2. **Position within the residual ambiguity (F1/F6).** Given the data, members of the same consistent fiber are indistinguishable; irreducible error = section half-diameter. Linear form: the coset $c_\beta+\ker G$, invisible exactly at queries with $\phi(x)\notin\operatorname{row}(G)$.
3. **Family structure off the covered set (F18/C13).** Within the linear exactly-$d$-dimensional class: at $x\notin U\cup D$ the consistent value set is all of $\mathbb R$ — or $\{0\}$ in one degenerate forced-zero case that vanishes for $\varepsilon>0$. The impossibility *recurs one level up*: no volume of archive at other locations transports information to an uncovered point without a declared member-level assumption.
4. **Sub-resolution distinctions (F1).** Members with traces within $2\varepsilon$ are indistinguishable; the noise floor $\tfrac12\omega(2\varepsilon)$ is irreducible.
5. **The gauge (CP §4.5, CP §2.3).** Parametrization, indexing, and even the member *set* of the family are unidentifiable (two families with identical finite windows can have different member sets); only window-level structure is ever identifiable.

**Why finite observations are insufficient.** $k$ observations cut the prior by codimension at most $k$: in the continuous category no scheme extracts more than $k$ independent dimensions of member identity (F20/CP-3; without continuity, dimension claims collapse — Borel bijection $\mathbb R^d\to\mathbb R$, C8). If the required-query quotient of $\mathcal F$ has cardinality $>\mathfrak c$, no finite-dimensional summary of any kind exists (pigeonhole).

**The exact source of ambiguity.** All of it is the geometry of one object: the **trace set** $T_{D,x}=\{(f|_D,f(x)):f\in\mathcal F\}\subseteq\mathbb R^{k+1}$ and its data-section $S_\varepsilon(\tilde y)$. The ambiguity at $(D,x,\tilde y)$ *is* the section; its worst-case diameter is $\omega(2\varepsilon)$; nothing else contributes (CP-1/CP-2: minimal sufficiency).

**Strongest structural constraints discovered.**
- $k\ge d$ is necessary for exact recovery of families containing a continuously and injectively parametrized $d$-cell (invariance of domain; F14(i)) — so $k\le5$ confines exact stable recovery to task dimension $d\le5$.
- Nonlinear families can defeat every design of size $d$: for $\sin(\theta x)$, $k=1$ and $k=2$ fail at **every** design (explicit collision pair $\theta=\frac{\pi}{2b}-\frac{\pi}{a}$, $\theta'=\frac{\pi}{2b}+\frac{\pi}{a}$), while $k=3=2d+1$ succeeds at designs with pairwise-irrational ratios (a generic set) — the $2d+1$ bound is attained (F14(ii)/C4). Generic-design global identifiability for analytic families is guaranteed at $k\ge2d+1$ (F12) — hence guaranteed only for $d\le2$ at this budget; sharpness for general $d$ is open (§9, item 3).
- Design-independence is topologically forbidden in multivariate domains at design size $d$: on any triod-containing domain, every $d\ge2$-dimensional space of *continuous* functions has singular size-$d$ designs (Mairhuber–Curtis–Sieklucki + swap argument, F10; with $k>d$ points, or discontinuous families, escapes exist — scoping is exact).
- Smoothness is not enough for generic identifiability: non-quasianalytic families have positive-measure failing design sets for every $k$ (C5); the dividing line is o-minimal definability + the identity-theorem property.
- Identifiability does not imply stability: families exist with $\omega(0)=0$ but $\omega(t)=\infty$ for all $t>0$ (tanh family, C3) — exact at $\varepsilon=0$, worthless under any noise.

---

## 3. Minimal assumptions for identifiability

**Definitions** (F-Part I): trace map $R_D f=f|_D$; trace set $T_{D,x}$; window $T_S=\{f|_S\}$; section $S_\varepsilon(\tilde y)$; modulus $\omega_{x,D}$ (and per-functional $\omega_\varphi$); radius $\rho=\tfrac12\omega(2\varepsilon)$; evaluation matrix $G_{ij}=\phi_j(x_i)$ for linear families; null set $N_D=\{v:v|_D=0\}$; generalized Lebesgue function $\Lambda_*(x)=\min\{\|w\|_*:G^\top w=\phi(x)\}$; power function $P_D(x)$.

**The assumption hierarchy — each level exactly equivalent to a grade of inference (the strongest form of minimality: every sufficient assumption implies it):**

| Level | Assumption | Equivalent to | Status |
|---|---|---|---|
| A0 | $\omega_{x,D}(0)=0$ — the trace set is the **graph of a function** $\Psi$ over $\mathcal F|_D$ | exact-data identifiability at $(D,x)$ | proven iff (F2) |
| A1 | $\omega_{x,D}(0^+)=0$ | stable inference; optimal modulus exactly $\varepsilon\mapsto\tfrac12\omega(2\varepsilon)$ | proven iff (F3); strictly stronger than A0 (C3) |
| A2 | linear family of dimension $d\le k$ + rank-$d$ design | **zero-error** recovery (linear category) | proven iff-price (F9); *not* necessary for stable approximate inference (RKHS ball) |

**Necessary conditions (all proven):** A0 for exact inference (tautologically minimal); $\omega(0^+)=0$ for stability; $\phi(x)\in\operatorname{row}(G)$ per query, $\operatorname{rank}G=d$ for the member (F6); $k\ge d$ ($d$-cell scope); archive: $n\ge d$ tasks and rank-$d$ conditions at every covered point — necessary *and* sufficient within the exactly-$d$ class with $d$ known (F17); coverage of the query by the archive or a declared member-level assumption (F18).

**Sufficient conditions (all proven):** Haar spaces on intervals make every $d$-point design identifying (equivalent to design-independence at size $d$); rank conditions above; for analytic separating $d$-parameter families on bounded connected domains, generic designs of size $2d+1$ identify globally, with a closed null exceptional set (F12), and bi-Lipschitz stability with infinitesimal separation (F13); RKHS balls and Lipschitz classes give stable *partial* information at every design (exact **minimax** radii: the power-function interval; min-distance — conditional, data-dependent radii can be strictly smaller, MP-6).

**Proven vs conjectural.** Everything in this section is proven under stated hypotheses. The corpus contains no conjectures. One quantitative item is OPEN (archive-noise constant, §9). Constrained families can be identifiable below full rank ($\ker G\cap(C-C)=\{0\}$, F6(iii)) — the necessity statements are scoped accordingly.

---

## 4. The identifiable mathematical object

**"What is the smallest mathematical object that must be recovered?"** Two objects, one per level; both minimality statements are proven.

**Meta level — the window system.**
- *Notation:* $\mathbb T=(T_S)_{S\ \text{finite}\subset U}$, $T_S=\{f|_S:f\in\mathcal F\}\subseteq\mathbb R^S$, with exact projective consistency $\operatorname{proj}_S T_{S'}=T_S$ for $S\subseteq S'$.
- *Meaning:* the family seen through all finite windows of the covered region — the constraint linking observations to queries. Under A0 a joint window collapses to the graph of the transfer function $\Psi_{D,x}$.
- *Information source:* the archive, and only on covered windows; identified exactly under F17's rank conditions (within the exactly-$d$ linear class, $d$ known); bounded off $U$ by the F18 dichotomy.
- *Degrees of freedom:* in the linear regime, those of a $d$-dimensional subspace restricted to $U$; in general, the sets themselves. Constraints: projective consistency is necessary but **not sufficient** for the system to come from any family (MP-1: countable coverage or compact windows suffice; Waterhouse's empty inverse limit is a consistent, unrealizable window system). At $\varepsilon>0$, windows are data-determined only up to closure (CP-2) — the canonical object is the *closed* system.
- *Identifiability conditions:* minimal sufficiency is exact — the inference problem at $(D,x)$ is a function of $T_{D,x}$ alone, and $T_{D,x}$ is recovered from exact-data sections; the $\varepsilon=0$ clause is load-bearing (CP-2).

**Task level — the fiber quotient.**
- *Notation:* the class of $f_\beta$ under agreement on $D$; linear closed form: $P_{\operatorname{row}(G)}c_\beta=G^+y$.
- *Meaning:* everything the support determines about the member — the functionals constant on its residual fiber.
- *Information source:* the $k$ observations, and nothing else.
- *Degrees of freedom:* at most $r=\operatorname{rank}G\le k$ continuous dimensions; per query set $Q$, exactly $\dim\operatorname{span}\{\phi(x):x\in Q\ \text{identifiable}\}\le r$ (CP-3; equality when those features span $\operatorname{row}(G)$).
- *Identifiability conditions:* existence of a data-computable summary is **equivalent to identifiability** — never an extra assumption (F20); continuous-category dimension bounds are the sandwich $d\le m_{\min}\le2d+1$ with the tripod counterexample showing $m_{\min}>d$ occurs (F19).

---

## 5. The abstract adaptation principle

**The principle (derived, not posited):** *meta-learning is identification of the closed window system on the covered region; adaptation is sectioning it at the support values; prediction is the unique conditionally-minimax selection — the section's Chebyshev center — always accompanied by its radius.*

**Learned from previous functions:** the constraint $\mathbb T$, and only that — gauge-free, distribution-free, learnable exactly on covered windows (F17), void beyond them (F18).

**Inferred from limited observations:** the selection — the section cut, worth at most $k$ continuous dimensions of member identity. Where prior sections are bounded and the support data realizable, adaptation decomposes exactly (MP-2): prediction $=$ zero-shot baseline $\operatorname{cen}(T_x)$ $+$ correction $\Delta$, with $|\Delta|\le\tfrac12\operatorname{diam}T_x$, $\Delta\equiv0$ at $k=0$, and worst-case guarantee gain $\tfrac12(\operatorname{diam}T_x-\omega(2\varepsilon))\ge0$ — zero exactly at saturation. Where prior sections are unbounded (including the unconstrained linear regime), no canonical baseline exists and any anchor is a gauge choice (MP-3); the anchored form acquires exact minimax content in particular for Euclidean centered-ball coefficient families with exact realizable data (projection of the anchor onto the task's constraint set), *not* precisely for them (sphere counterexample).

**How query-dependent change arises — as a theorem, not a choice:** support and query meet inside one geometric object, the joint window $T_{D\cup\{x\}}$. Consequently: the validity region $\{x:\phi(x)\in\operatorname{row}(G)\}$ depends on the support configuration; the read-out of the shared state is query-dependent (sensitivity profile $w(x)=(G^+)^\top\phi(x)$, summing to $1$ on the validity region when constants lie in the family); the per-query minimal state varies with the query; the certificate (conditional radius) varies with the query; and in regimes with no finite continuous state at all (Lipschitz class), adaptation information is carried irreducibly by the query-coupled section.

**Two certified qualitative properties of optimal adaptation:** guarantees tighten monotonically with evidence while estimates move non-monotonically (OP-8); and the optimal operator is genuinely discontinuous in the support values at section-topology transitions, so continuous surrogates carry irreducible localized error (MP-4).

---

## 6. Computational interface

$$\mathbb A:\ \underbrace{\mathbb T}_{\text{previously learned object}}\times\underbrace{S_t=\{(x_i,\tilde y_i)\}_{i\le k}}_{\text{limited observations (a finite set)}}\times\underbrace{x}_{\text{query}}\times\underbrace{\varepsilon}_{\text{noise level}}\ \dashrightarrow\ \underbrace{\big(\operatorname{cen}S_\varepsilon(\tilde y),\ \tfrac12\operatorname{diam}S_\varepsilon(\tilde y)\big)}_{\text{estimated value + certificate}}$$

- **No unary form exists:** $A(S_t,x)$ without the family argument has infinite worst-case error (OP-10); the family argument is minimally-sufficiently the trace set (CP-2; minimality holds with the $\varepsilon=0$ clause — at fixed $\varepsilon>0$ the trace set is determined only up to dense-in-section modifications).
- **The target operator is not linear.** Even for linear families the canonical operator is nonlinear once $k\ge3$ (the two-point midrange being linear is a $k=2$ coincidence); linear *optimal selections* exist in the linear regime, but they are not the canonical, pointwise-dominant operator (OP-4a/OP-4b).
- **The output is a pair.** The point prediction alone loses the certificate; the change form is $\Delta=\mathbb A-\text{baseline}$ where MP-2 applies.
- **Partiality is structural** (three sources, to be surfaced, never silently extrapolated): unrealizable data (doubles as misspecification detector); off-coverage configurations; unbounded sections at covered-but-unidentifiable queries (radius $+\infty$ is the required output).
- **$\varepsilon$ is an argument:** no single rule is optimal across noise levels (F1 Rem. 1.4).

**What a future model must approximate:** the **closed window system on the covered region** — operationally, *conservative outer envelope pairs* $(\inf S_\varepsilon,\sup S_\varepsilon)$ per encountered configuration, from which center and radius follow. One-sided semantics is mandatory: outer enclosures keep certificates valid (merely conservative); any under-approximation is a false certificate. A finite archive gives only *inner* window approximations; a declared model-class closure assumption (e.g. exactly-$d$ linear) is what converts inner to outer — *no closure assumption, no valid radius*. Necessary properties of any admissible approximation: axioms A1–A10 of `operator_specification.md` (permutation invariance; affine equivariance including reflections; gauge invariance; $\varepsilon$- and support-monotonicity of the guarantee with the center exempt; projective consistency; reproduction on exact identifiable data; 1-Lipschitz selection; conditional sensitivity-sum; partiality surfacing). Error calculus: section error transfers to prediction error with constant 1 (tight); window error transfers through the modulus, $\tfrac12\omega(2\varepsilon+2h)+h$, possibly infinite (tanh); archive-noise constant OPEN.

---

## 7. Failure conditions

Full certified catalogue: `03_meta_learning_principle/failure_modes.md` (certified modes in groups A–F plus degenerate-case conventions in group G, each with certificate and observable signature). Headlines:

- **Adaptation impossible:** no shared structure (F4); query outside the validity region (section $=\mathbb R$); rank failure at the actual support configuration (size-$d$ design degeneracy in multivariate domains, F10); query off the covered set (F18); support size below the $d$-cell dimension (F14(i)); nonlinear global ambiguity at every configuration for $d\le k<2d+1$ (C4); non-quasianalytic flat directions with positive-measure failing configurations (C5); saturation $\omega(2\varepsilon)=\operatorname{diam}T_x$ (MP-2).
- **Assumptions break:** archive diversity below $d$ (rank failure — realized ambiguity, F17); unknown/misspecified $d$ (rank certifies $\ge d$ only); noisy archive (constant OPEN); misspecified $\varepsilon$ (two distinct failure directions); identifiable-but-unstable families under any noise (C3); consistent-but-unrealizable learned window systems (MP-1(iii)).
- **Degenerate solutions:** forced-zero archive case (C13); residual gauge at small $k$ — only invariants identifiable, absolute outputs fabricated (F16); singleton family (adaptation trivially zero); query in support (floor $\varepsilon$, not zero — F4); Borel-collapse pseudo-compressions that fail continuity audits (C8); scalar design scores mis-ranking supports (C12); monotone-update constraints contradicting optimal adaptation (OP-8).

---

## 8. Design constraints for future models — checklist

A future computational model **must**:

1. **Preserve the identifiable object.** Approximate the closed window system / conservative envelope pairs — not a point-prediction map alone (which forfeits the certificate and under-determines the object, OP-9).
2. **Not replace it with unconstrained embeddings.** Any internal representation of a task is subject to the proven dimension laws: at most $k$ continuous dimensions from $k$ observations; sandwich $[d,2d+1]$ with topological obstructions (tripod); *no* finite continuous representation for Lipschitz-type families — the model must then operate on sections/intervals, not vectors. An embedding that ignores these is provably lossy or fabricating; and the gauge (which coordinates) is unidentifiable, so no meaning may be attached to embedding coordinates beyond window-level invariants.
3. **Distinguish calibration from genuine adaptation.** The exact separation exists where prior sections are bounded and the data realizable (MP-2): zero-shot-attainable accuracy (baseline $\operatorname{cen}T_x$; plus gauge-invariants identifiable at $k=0$, F16) is family-level calibration, not adaptation; genuine adaptation is the section-narrowing correction $\Delta$ (with $|\Delta|\le\tfrac12\operatorname{diam}T_x$), whose worst-case **guarantee gain** is exactly $\tfrac12(\operatorname{diam}T_x-\omega(2\varepsilon))$. A model's claimed few-shot gains must be measured against the baseline, and vanish at saturation.
4. **Pass the mathematical falsification tests** (P1–P10, `meta_learning_abstraction.md` §6.9 — all adversarial protocols with computable constants): radius floor (P1); $k$-dimension capacity ceiling (P2); sensitivity sum (P3); query-dependent degradation curves (P4); support-placement ratios via exact small-$k$ constants (P5); off-coverage collapse under the F18 modification (P6); the $2d+1$ support-size cliff on the $\sin$ family (P7); noise-fragility divergence on tanh-type families (P8); discontinuity spikes (P9); non-monotone optimal updates (P10).
5. **Ship the radius with every prediction**, take $\varepsilon$ as input, surface the three partiality sources, and declare the model-class closure assumption that legitimizes its certificates (§6).
6. **Respect the scope:** all guarantees are worst-case; benchmark averages neither confirm nor refute the theory (failure-mode F3); multi-query evaluation is exact only under sup-loss.

---

## 9. Open questions

**Unresolved mathematical questions (the corpus's honest edges):**
1. **Archive-noise constant (OPEN, flagged in the corpus).** A refereeable quantitative bound for the composed operator under inexact archives: fix the subspace metric (e.g. largest principal angle), the coefficient norm, the smallness regime, and derive the explicit constant $C(x,G)$ with remainder — or prove a matching lower bound.
2. **Exact class for anchored-projection optimality.** Centered balls suffice; spheres also work (referee's counterexample). Characterize exactly the families whose sections are symmetric about the projected-anchor value.
3. **Sharpness of $2d+1$ for general $d$.** The $\sin$ family attains it at $d=1$. Exhibit analytic families attaining it for every $d$, or improve the generic bound below $2d+1$ for structured subclasses.
4. **Archive identification beyond the common core.** Deterministic overlap-chaining conditions replacing the common-core hypothesis in F17, with an exact combinatorial characterization.
5. **Realization of window systems between countable and compact.** MP-1 gives two sufficient regimes and a counterexample; the exact boundary (which index sets / which window classes) is open.
6. **Structure of the section-containment preorder.** Order-theoretic properties of the Blackwell-type support comparison (chains, lattice structure, existence of maximal $k$-point supports).
7. **Joint-loss constants.** Exact minimax constants for multi-query recovery under norms other than sup-loss (the radius/diameter gap is bounded but not pinned).
8. **$\varepsilon$-adaptive rules.** Whether a single rule can be simultaneously near-optimal across all noise levels, and with what factor (the corpus proves only per-$\varepsilon$ optimality).
9. **The F18 dichotomy beyond the linear exactly-$d$ class** — the correct general statement off the covered set for structured nonlinear classes.

**Assumptions requiring empirical validation when the framework is instantiated** (validation targets, not mathematical gaps): boundedness of prior sections (baseline existence); the effective task dimension $d$ and the adequacy of the exactly-$d$ closure class; the bounded-per-coordinate noise model versus stochastic alternatives; archive coverage of the intended query region; quasianalyticity-like regularity of the actual task family.

**Possible approximation strategies (mathematical directions only — no architectures):**
- approximate the *envelope pair* under the outer-enclosure constraint, with tightness as the quality objective and validity as a hard constraint;
- obtain outer windows from inner archive evidence via declared finite-dimensional closure classes (F17's mechanism), with the class dimension audited by rank certificates ($\ge d$ only);
- represent the meta-object as a discretized closed window system on covered configurations with projective consistency enforced as a constraint (MP-1's realizability caveats observed);
- realize the two-factor composition $A\circ\widehat T$ explicitly, propagating the OPEN flag of the archive-noise constant;
- exploit the regime closed forms as exactly solvable substructures (min-dual-norm weights; kernel interpolant; envelope averages) wherever the corresponding hypotheses are declared and checked.

---

**End of frozen specification.** A future system is correct relative to this document iff it approximates the object of §4, through the interface of §6, within the constraints of §8 — and its claims are refutable exactly by the protocols cited there. Nothing in this document says how to build it; everything in it says what "it" must be.
