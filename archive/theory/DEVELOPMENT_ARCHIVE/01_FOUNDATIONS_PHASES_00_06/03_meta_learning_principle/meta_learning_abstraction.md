# The Abstract Meta-Learning Principle

> **Status:** Phase-3 derivation, 2026-08-02. Sources of record: `../00_raw_outputs/identifiability_treatise_raw.md`, `../01_math_foundation/theorem_summary.md` (F-numbers, C-counterexamples), `../02_theory_refinement/` (CP-, OP-numbers). New results in this phase carry MP-numbers; every MP was adversarially refereed and the corrections are incorporated. Companion files: `operator_specification.md` (the interface), `failure_modes.md` (certified failure catalogue).

## 0. Scope and conventions

- **Task $=$ family member.** A task $t$ is a member $f_t\in\mathcal F\subseteq\mathbb R^{\mathcal X}$. A new task provides a support set $S_t=\{(x_i,y_i)\}_{i\le k}$, $k\le5$, with $\max_i|y_i-f_t(x_i)|\le\varepsilon$, $\varepsilon$ known. A query $x$ asks for $f_t(x)$. The *meta-level data* is the archive: finitely many previous tasks, each observed on its own finite support; covered set $U$.
- **Distribution-free, worst-case.** The theory contains no probabilistic assumption: the archive is a finite set of tasks, not a sample from a task distribution, and **every certificate below is minimax** (sup over consistent tasks and noise). Average-case performance under a benign task distribution is outside the theory; the falsifiable predictions of §6 are therefore stated as adversarial protocols, not benchmark expectations.
- **Scalar target.** The exact constants hold for one scalar query at a time; the certified multi-query bridge is sup-loss, where the per-query rule remains exactly optimal with value $\sup_x\tfrac12\omega_{x,D}(2\varepsilon)$ (F1 Rem. 1.1). Other joint losses carry only radius $\in[\tfrac12\mathrm{diam},\mathrm{diam}]$.
- Notation: window $T_S=\{f|_S:f\in\mathcal F\}$; joint window $T_{D,x}$; section $S_\varepsilon(\tilde y)$; modulus $\omega_{x,D}$; center $\operatorname{cen}(S)=\tfrac12(\inf S+\sup S)$.

---

## 1. From function identification to task adaptation

The four objects, derived — none of them is assumed to be a parameter, embedding, or latent vector.

### 1.1 Shared across tasks: the window system

The object shared across tasks is the **projectively consistent system of windows**
$$\mathbb T\;=\;\big(T_S\big)_{S\ \text{finite}\subset U},\qquad \operatorname{proj}_S T_{S'}=T_S\ \ (S\subseteq S'),$$
the family seen through all finite windows of the covered region. The consistency identity is exact (image of image; refereed). This — not a parameter vector, not an embedding — is what previous tasks can contribute: by CP-1/CP-2 the inference problem at any $(D,x)$ is a function of $T_{D,x}$ alone and of nothing finer, and by CP §2 the archive can reveal windows only on $U$.

**MP-1 (Realization theorem for window systems — new; refereed: confirmed).** Consistency is necessary but *not sufficient* for a system to be a family:
(i) if the covered region is **countable**, every projectively consistent system of nonempty windows is realized by a family (countable directed posets have cofinal chains; a surjective inverse sequence of nonempty sets has nonempty limit with surjective limit projections; the thread family has exactly the given windows);
(ii) if every window is **compact**, the same holds over arbitrary index sets (inverse limits of nonempty compacta);
(iii) in general, **no**: Waterhouse's empty inverse limit (Proc. AMS 36 (1972) 618) is verbatim a consistent window system — $T_S=\{$injections $S\to\mathbb N\}$ over the finite subsets of an uncountable $\mathcal X$ — realized by no family whatsoever. (Historical precedents: Higman–Stone 1954; Aronszajn; Henkin.)
*Caveat (CP-2):* at any fixed $\varepsilon>0$ windows are determined by data only up to dense modifications; realization statements should be read for closed systems.

### 1.2 Specific to one task: the fiber position

The task-specific object is the **position of $f_t$ within the constraint**: its equivalence class under agreement on observed points — never a coordinate vector per se. Coordinates for it exist exactly in the finite-dimensional continuously-embeddable regime (F19: sandwich $d\le m_{\min}\le2d+1$; tripod counterexample C6 shows $m_{\min}>d$ happens; the gauge — which coordinate system — is *never* identifiable, CP §4.5).

### 1.3 Inferable from a small support set: the quotient, at most $k$ deep

From $S_t$ one can infer **exactly the functionals of the task constant on its residual fiber** (the consistent set of members). In the linear regime this is the projection $P_{\operatorname{row}(G)}c_t=G^+y$ — computable from the support (refereed identity), living in a space of dimension $r=\operatorname{rank}G\le k$. In the *continuous category* no scheme extracts more than $k$ independent dimensions of task identity from $k$ observations (F20/CP-3; without continuity, dimension claims are void — Borel collapse C8). Minimality is query-relative: for a query set $Q$ the minimal continuous state has dimension $\dim\operatorname{span}\{\phi(x):x\in Q\ \text{identifiable}\}\le r$, with equality when the identifiable queries' features span $\operatorname{row}(G)$ (e.g. when reconstruction of the support itself is among the queries) — refereed correction.

### 1.4 Impossible to infer: the five certified voids

(i) anything off the support for unconstrained families (F4 — information lives exactly on $D$); (ii) the position within a residual section — irreducible error $=$ section half-diameter (F1); (iii) family structure at locations never covered by any previous task — within the linear exactly-$d$ class, the consistent value set off $U$ is $\mathbb R$ or the forced-zero $\{0\}$ (F18); for general families the correct statement is (i) applied one level up; (iv) sub-resolution distinctions — traces within $2\varepsilon$ are indistinguishable, floor $\tfrac12\omega(2\varepsilon)$ (F1); (v) the gauge — parametrization, indexing, and even the member *set* of the family (CP §2.3: identical finite windows, different member sets) are unidentifiable; only $\mathbb T$-level structure is.

---

## 2. The abstract meta-learning operator

### 2.1 The implied object

No unary operator $A(S_t,x)$ exists: against $\mathcal F=\mathbb R^{\mathcal X}$ any support-only rule has infinite worst-case error off the support (OP-10). The theory implies a **partial binary operator with a certificate output**:
$$\boxed{\ \mathbb A\big(\mathbb T;\;S_t,\;x,\;\varepsilon\big)\;=\;\Big(\operatorname{cen}\big(S_\varepsilon(\tilde y)\big),\ \tfrac12\operatorname{diam}S_\varepsilon(\tilde y)\Big)\ }$$
— section the joint window $T_{D\cup\{x\}}$ at the support values, output the center **and the radius**. It is minimax-optimal for every family realizing $\mathbb T$ (F1) and is the *unique* conditionally-minimax rule on realizable data when $\omega(2\varepsilon)<\infty$ (OP-3). Its three partiality sources are structural, not defects: unrealizable support data; queries outside covered windows; unbounded sections at covered-but-unidentifiable queries (center undefined, radius $+\infty$ is the honest output).

### 2.2 "How a function should change when evidence arrives"

The mandate's distinction — not a universal function, but a law of change — is realized exactly, with a boundedness hypothesis the general theory forces:

**MP-2 (Baseline–correction decomposition; refereed: confirmed under its hypotheses).** Let the prior value set $T_x=\{f(x):f\in\mathcal F\}$ (the $k=0$ section) be nonempty and **bounded**, and the support data realizable. Then with $B(x)=\operatorname{cen}(T_x)$ (the best zero-shot prediction: the $k=0$ minimax value is exactly $\tfrac12\operatorname{diam}T_x$, since $\omega_{x,\emptyset}(t)=\operatorname{diam}T_x$):
$$\mathbb A=\;B(x)\;+\;\Delta(\mathbb T;S_t,x,\varepsilon),\qquad |\Delta|\le\tfrac12\operatorname{diam}T_x,\qquad \Delta\equiv0\ \text{at}\ k=0,$$
and the **worst-case guarantee gain** of adaptation is
$$\tfrac12\big(\operatorname{diam}T_x-\omega_{x,D}(2\varepsilon)\big)\;\ge\;0,$$
zero precisely at *saturation* $\omega(2\varepsilon)=\operatorname{diam}T_x$ (noise so large the support constraints never bind).
*Scope:* when $T_x$ is unbounded — which includes the unconstrained linear regime, where $T_x=\mathbb R$ — no canonical baseline exists; see MP-3.

**MP-3 (Anchors are gauge choices; refereed hypotheses).** A canonical baseline is definable iff the prior sections are bounded ($B=\operatorname{cen}T_x$). Otherwise any anchored decomposition $\mathbb A=g_0(x)+\Delta$ requires choosing $g_0$, and that choice is a gauge fixing — certified unidentifiable (§1.4(v)). The anchored form acquires minimax content in a specific regime: for coefficient families that are Euclidean balls centered at $c_0$, with exact realizable data, the canonical prediction equals the value of the **projection of the anchor onto the task's constraint set** ($c_0+G^+(y-Gc_0)=P_{\{c:Gc=y\}}(c_0)$; Golomb–Weinberger center). This holds *in particular* for centered balls — **not** "precisely": the sphere $\{\|c-c_0\|=R\}$ has sections likewise symmetric about the projected anchor (referee's counterexample), so ball-ness is sufficient, not necessary.

Two certified phenomena about the law of change:

- **Guarantee-monotone, estimate-non-monotone (OP-8).** Additional support shrinks every section, so the certificate tightens monotonically — but the estimates need not move monotonically (nested sections $[0,10]\supset[9,10]\supset[9,9.2]$ have centers $5,\,9.5,\,9.1$). Adaptation is not a pointwise contraction toward truth.
- **MP-4 (Discontinuity obstruction — new).** $\mathbb A(\mathbb T;\cdot,x,\varepsilon)$ is genuinely discontinuous in the support values: for the two-member family $f\equiv0$, $g$ with $g(x_1)=1,g(x)=c$, at noise $2\varepsilon\ge1$ the center jumps from $0$ to $c/2$ as $\tilde y_1$ crosses the boundary where $g$ enters the consistent set. Consequently **any continuous approximator of $\mathbb A$ has irreducible sup-error $\ge$ half the jump** near section-topology transitions. This is a structural property of optimal adaptation, not an artifact.

---

## 3. Meta-level versus task-level information

**Minimum from the meta side.** For a single encounter $(D,x)$: exactly $T_{D,x}$ — minimally sufficient, with the $\varepsilon=0$ clause (CP-2). For all encounters within the covered region: the closed window system $\mathbb T$ on $U$ (§1.1), obtainable from the archive exactly under F17's rank conditions and bounded off $U$ by F18.

**Minimum from the task side.** The support values modulo section-equality — the data enter only through $S_\varepsilon(\tilde y)$. Linear regime, exact data: $G^+y\in\operatorname{row}(G)$, dimension $r\le k$; per query set $Q$, dimension $\dim\operatorname{span}\{\phi(x):x\in Q\ \text{identifiable}\}\le r$ (CP-3, with the refereed span condition for equality).

**The menu verdict** (mandate: do not pre-select). What the theory implies is, primarily, an item **not on the menu**:

| Candidate object | Status under the theory |
|---|---|
| **(other) set-valued enclosure correspondence + canonical selection** | **The primary implied object** — assumption-free stratum (OP §7.1–7.2): the correspondence $(\mathbb T,S_t,x)\mapsto S_\varepsilon(\tilde y)$ with its unique conditionally-minimax selection and radius. Everything below derives from it. |
| constraint set | Yes — the per-encounter face of the primary object: adaptation *is* constraint intersection (§4). |
| conditional operator | Yes — the selection is conditional estimation with set-theoretic (not probabilistic) conditioning. |
| sufficient statistic | Exists iff identifiable (F20 — equivalence, never an extra assumption); minimal dimensions known and query-relative (CP-3). |
| task coordinate system | Only in the finite-dimensional continuously-embeddable regime; dimensions sandwiched $[d,2d+1]$, tripod obstruction (C6); the gauge never identifiable. |
| task transformation | Only for gauge families ($g_0+$ group orbit; F16), with exact observation counts for killing the residual gauge. |
| update direction | Derived object $\Delta$ (MP-2) where a baseline exists. *Refereed correction:* $\Delta$ is **not** confined to identifiable directions — in partial-information regimes (Lipschitz, RKHS ball) the center moves at queries where nothing is identifiable; components of $\Delta$ along ambiguous directions exist but carry no certificate beyond the radius. |

---

## 4. Few-shot adaptation geometry

**Identifiable directions of change.** Linear regime: the quotient $\operatorname{row}(G)\cong V/N_D$ — at most $k$ dimensions. **Ambiguous directions:** $N_D=\{v\in V:v|_D=0\}$, the functions invisible on the support, dimension $d-r$; generally, the residual section spread $\omega$. The ambiguous component is not "unmoved" by adaptation (see MP-2 correction above) — it is *uncertified*.

**How the support constrains change.** Each observation is one section constraint; $k$ observations cut the prior set by codimension at most $k$ (linear: exactly $r$). The cut is geometry-, not count-, dominated: designs of equal size differ in yield by exact, computable factors (Lebesgue-type constants; at the theory's own budget $k\le5$ the polynomial-node gap is a modest computable constant, *not* the asymptotic exponential — refereed scope). The complete comparison of two supports is the **Blackwell-type section-containment preorder** (every $D_1$-section inside some $D_2$-section), of which the radius $\rho$ is the canonical scalar quotient — and strictly coarser (C12: equal radii, inequivalent supports).

**How the query interacts with the inferred change (MP-6, refereed).** The transported object is the *joint* window $T_{D\cup\{x\}}$ — support and query enter one geometric object, and this is why adaptation is query-dependent as a theorem, not a design choice:
1. the **validity region** $\{x:\phi(x)\in\operatorname{row}(G)\}$ depends on the support configuration; for $r<d$ it is generally proper;
2. the **read-out is query-dependent**: value at $x$ = pairing of the shared state $G^+y$ with $\phi(x)$; the sensitivity of the prediction to each support value is the $x$-dependent weight vector $w(x)=(G^+)^\top\phi(x)$ (weights summing to $1$ on the validity region when constants lie in the family — refereed);
3. the **per-query minimal state** varies with the query (CP-3 strictness: $\dim1<r=2$ example);
4. in regimes with **no finite continuous state at all** (Lipschitz class, F19(v)), adaptation information is carried irreducibly by the query-coupled section — there is nothing query-independent to transport;
5. the **certificate is query-dependent**: conditional radius at $x$, bounded by regime quantities ($\varepsilon\Lambda_*(x)$; $P_D(x)$-interval; for Lipschitz: conditional radius $\le\min_i\operatorname{dist}(x,x_i)$ with equality in the worst case over data — *refereed correction: the min-distance is the minimax radius, the conditional radius can be strictly smaller, even zero*).
6. **degenerate query positions:** for $x=x_i$ (query in support) the floor is $\varepsilon$, not $0$, for rich families (F4) — even memorization has a worst-case cost under noise.

**Saturation.** Where MP-2 applies, the adaptation gain $\tfrac12(\operatorname{diam}T_x-\omega(2\varepsilon))$ quantifies exactly when few-shot evidence helps at all; at $\omega(2\varepsilon)=\operatorname{diam}T_x$ the support is provably worthless at $x$ in the worst case.

---

## 5. Relation to existing abstract learning concepts

Comparisons only now, each with its demonstrated correspondence and its scope; no equivalence is claimed beyond what is proved.

- **Initialization-based adaptation.** Demonstrated: in the Euclidean linear regime with exact realizable data, the canonical prediction for a centered-ball family equals the value of the *projection of the meta-anchor onto the task's constraint set* (MP-3). Scope: sufficiency only (sphere counterexample kills "precisely"); for unbounded prior sections the anchor itself is an unidentifiable gauge choice. No general equivalence.
- **Metric-based adaptation.** Demonstrated: in the 1-Lipschitz regime the canonical operator *is* a distance-based rule — the average of the two McShane–Whitney envelopes, with minimax radius $\min_i\operatorname{dist}(x,x_i)$ (exact data; $\varepsilon>0$ via inflated envelopes). Scope: that regime.
- **Posterior-style inference.** Demonstrated: at $\varepsilon=0$ on the RKHS unit ball, the canonical center is the kernel interpolant — formally the Gaussian conditional mean (identity of formulas, not models: Driscoll; and under the theory's $\ell_\infty$ noise the canonical center is *not* the noisy-GP posterior mean). Scope: exact-data formula identity only.
- **Conditional function estimation.** The closest general match, by definition rather than coincidence: the canonical operator is conditional estimation with **set-theoretic conditioning** (the section), valid with no structure at all.
- **Operator learning.** The learnable unknown is not the sample-to-value operator but the **window system**; the operator is derived from it (OP-9). Learning the point-prediction operator alone forfeits the certificate (see `operator_specification.md`).
- **Low-dimensional task representation.** Exists iff the finite-dimensional continuous regime holds; dimensions are the sandwich $[d,2d+1]$ (not "$=d$": tripod), and data-computable states are query-relative (CP-3). For the Lipschitz class no finite continuous representation exists at all — yet the principle still operates, via sections.

---

## 6. Theoretical requirements

1. **Mathematical definition.** The principle: *meta-learning is identification of a projectively consistent window system $\mathbb T$ on the covered region; adaptation is sectioning $\mathbb T$ at the support values; prediction is the unique conditionally-minimax selection, always accompanied by its radius.* Formally the operator of §2.1.
2. **Required assumptions.** Base stratum: none beyond nonempty $\mathcal F$, known $\varepsilon$, closed $\ell_\infty$ noise, scalar target (sup-loss bridge for multi-query). Regime strata: as tabulated in `theorem_formalization.md` Part 0/F7/F12/F17. Meta stratum: archive exactness, class dimension $d$ known, rank conditions (F17); $n\ge d$ tasks necessary.
3. **Identifiable object.** Meta: closed window system on $U$ (exactly; F17/F18). Task: the fiber quotient — linear form $G^+y$, dimension $r\le k$.
4. **Adaptation object.** The section-narrowing; where prior sections are bounded, the correction $\Delta$ of MP-2 with gain identity and saturation threshold; anchored variants only as gauge choices (MP-3).
5. **Query prediction rule.** $\operatorname{cen}(S_\varepsilon(\tilde y))$ with radius; closed forms: $\phi(x)^\top G^+\tilde y$ (linear), kernel interpolant (RKHS, $\varepsilon=0$), envelope average (Lipschitz).
6. **Necessary conditions.** $\omega_{x,D}(0)=0$ for exact adaptation at $(D,x)$; $\omega(0^+)=0$ for stability; $k\ge d$ for families containing a continuously parametrized $d$-cell (scoped — constrained families can identify below rank $d$, F6(iii)); archive rank conditions; **budget consequences of $k\le5$:** exact stable recovery only for task-dimension $d\le5$; generic-design global identification for nonlinear analytic families only for $d\le2$; for $d\in\{3,4,5\}$ nonlinear, identifiability must be verified per support configuration (no generic guarantee).
7. **Failure conditions.** The certified catalogue in `failure_modes.md` (structural, configurational, coverage, stability, representational, evaluational).
8. **Degenerate cases.** $k=0$: $\mathbb A=$ baseline, minimax $\tfrac12\operatorname{diam}T_x$; query in support: floor $\varepsilon$ (F4); singleton family: everything identifiable at $k=0$, $\Delta\equiv0$ — the pure-baseline case; empty family: excluded by convention (no realizable data, operator nowhere defined); repeated/coincident support points: rank drop; overdetermined consistent data $k>d$: realizability check doubles as misspecification detector; forced-zero archive case (C13).
9. **Falsifiable mathematical predictions** — all worst-case; each is an explicit adversarial protocol with computable constants:
   - **P1 (radius floor).** For computable $\omega$: construct two members at trace distance $\le2\varepsilon$ with query gap near $\omega(2\varepsilon)$, feed both the midpoint data; *any* system errs $\ge\approx\tfrac12\omega(2\varepsilon)$ on one. Falsifier: beating the floor on both.
   - **P2 (capacity ceiling).** Linear task family of dimension $d>k$: no system reliably recovers more than $k$ continuous-linearly-independent task functionals from $k$ observations; a system appearing to must be importing member-level assumptions — which the protocol forces it to declare.
   - **P3 (sensitivity sum).** If constants lie in the family, translation equivariance forces any locally Lipschitz optimal rule to have a.e. sensitivities summing to $1$ (derived from OP-7(ii), $\alpha=1$); test by finite differences along $\mathbf1$.
   - **P4 (query-dependence curves).** Error at $x$ under adversarial members tracks the computable certificate curve ($\varepsilon\Lambda_*(x)$ / $P_D(x)$ / envelope width) — not any query-independent constant.
   - **P5 (support placement).** At fixed $k\le5$, the worst-case error ratio between two node configurations equals the ratio of their exact (computable, small-$k$) Lebesgue-type constants; adversarial member: the Lebesgue extremal.
   - **P6 (off-coverage collapse).** Apply the F18 modification at an uncovered $x$: all archive and support data unchanged, query value arbitrary — any system's claimed accuracy off $U$ is refuted by construction unless it declares the member-level assumption doing the work.
   - **P7 (the $2d+1$ phenomenon).** $\sin(\theta x)$ tasks: at $k=2$, for *every* support $\{a,b\}$ the explicit pair $\theta=\frac{\pi}{2b}-\frac{\pi}{a},\ \theta'=\frac{\pi}{2b}+\frac{\pi}{a}$ yields identical supports and different query values — irreducible ambiguity; at $k=3$ with pairwise-irrational ratios it vanishes. Systems must show exactly this support-size cliff.
   - **P8 (fragility).** Tanh-type families: exact adaptation at $\varepsilon=0$, error diverging as the parameter grows for any fixed $\varepsilon>0$ — protocol: drive $|\theta|$.
   - **P9 (discontinuity spikes).** Two-member boundary family of MP-4: any continuous system carries sup-error $\ge$ half the jump localized at transition supports.
   - **P10 (non-monotone updates).** Nested-section protocol ($[0,10],[9,10],[9,9.2]$): optimal estimates move non-monotonically while certificates tighten; systems constrained to monotone updates deviate from conditional minimaxity by a computable margin.

---

## 7. Summary of the principle

*What is meta-learned* is a constraint — the closed, projectively consistent window system, the only family-level object that exists, is gauge-free, and is learnable, and only on covered ground. *What adaptation does* is intersect that constraint with $k\le5$ pieces of evidence, which buys at most $k$ continuous dimensions of task identity and provably nothing off the support without the constraint. *What prediction is* is the unique conditionally-minimax selection from what remains, priced by a radius that is itself part of the output. *Why adaptation is query-dependent* is that support and query meet inside one geometric object — the joint window — so validity, read-out, and certificate all vary with the query as theorems. Everything a future system must approximate is fixed by this; nothing about how to approximate it is.
