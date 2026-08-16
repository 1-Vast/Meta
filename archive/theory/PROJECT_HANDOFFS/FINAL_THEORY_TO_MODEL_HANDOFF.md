# FINAL THEORY-TO-MODEL HANDOFF

> **Frozen complete specification. 2026-08-03.** This document compresses Phases 0–6 into the minimum complete mathematical specification a future researcher needs. It derives nothing new and opens nothing new; every statement is a frozen, adversarially-refereed result of the corpus, cited to its origin. All guarantees are **worst-case (minimax), distribution-free, per-query (sup-loss)**.
>
> **Corpus map.** Phase 0 raw treatise + Phase 1 foundation (`F`-, `C`-numbers). Phase 2 refinement (`CP`-, `OP`-, `F`-numbers). Phase 3 meta-principle (`MP`-, `P`-numbers). Phase 4 differentiable operator (`DM`-, `IB`-, `A`-axioms, `NP`-numbers). Phase 5 realizability (`CR`-numbers). Phase 6 conditional/imperfect archive (`CI-A/B/C/D`-, `T`-numbers). Section references point into `00_raw_outputs/` … `06_conditional_imperfect_archive/` and `RESEARCH_HANDOFF.md`.
>
> **Epistemic tags used throughout:** **[P]** proved · **[C]** conditional on stated assumptions · **[X]** impossible · **[O]** open but non-blocking.

---

## 1. The original abstract problem

A family $\mathcal F=\{f_a:\mathcal X\to\mathbb R\}$ on a set $\mathcal X$ with no assumed structure. Three data strata:
- **Archive:** many prior members, each observed on its own finite design $D_a\subset\mathcal X$, values possibly noisy ($|y_{ai}-f_a(x_{ai})|\le\delta_a$); optional auxiliary quantity $c_a\in C$ per member; covered set $U=\bigcup_a D_a$.
- **New-member support:** for a new member $f_b$, its auxiliary $c_b$ known, plus $\le k\le5$ evaluations $\tilde y$ on a design $D$ with $\|\tilde y-f_b|_D\|_\infty\le\varepsilon$ ($\varepsilon$ known; $\varepsilon=0$ = exact).
- **Query:** $f_b(x)$ required at $x\notin D$ — *query-dependent inference*.

**Task:** determine, for each query, whether it is answerable, with what accuracy, and with what certificate — deriving (never positing) the structure a learning system would need.

---

## 2. The central impossibility result

**Theorem 1 (exact minimax). [P]** For any nonempty $\mathcal F$, design $D$, query $x$, noise $\varepsilon$, the minimax error of *every* estimator (randomization included) is exactly
$$\inf_\Phi\sup_{f,\tilde y}\,|\Phi(\tilde y)-f(x)|\;=\;\tfrac12\,\omega_{x,D}(2\varepsilon).$$
**Corollary (strongest impossibility, F4/C1). [P]** For $\mathcal F=\mathbb R^{\mathcal X}$ (or any family with pairs agreeing on $D$ but differing arbitrarily at $x$): error is $\varepsilon$ on $D$ and $+\infty$ everywhere else. **Information exists exactly on the design.** Membership in a family is worthless unless the family constrains the joint trace geometry; the impossibility recurs at the archive level (off the covered set, F18/CI-C).

---

## 3. The single primitive underlying the theory

**The trace modulus.** For the joint window $T_{D,x}=\{(f|_D,f(x)):f\in\mathcal F\}\subseteq\mathbb R^{k+1}$,
$$\omega_{x,D}(t)=\sup\{|f(x)-g(x)|:f,g\in\mathcal F,\ \|f|_D-g|_D\|_\infty\le t\}.$$
Everything in the program is a statement about this one object and its sections $S_\varepsilon(\tilde y)$. Minimal sufficiency (CP-1/CP-2, **[P]**): the inference problem at $(D,x)$ is a function of the trace set alone — with the $\varepsilon=0$ clause load-bearing (at fixed $\varepsilon>0$ the trace set is pinned only up to dense-in-section modification). The three Phase-6 additions are three reductions to this same primitive.

---

## 4. The canonical conditional operator

**Definition (CI-D1). [P under §8 hypotheses]**
$$A_c(\text{archive},c_b,S_b,x,\varepsilon)\;=\;\Big(\operatorname{cen}\,S_\varepsilon(\tilde y\mid c_b),\ \tfrac12\operatorname{diam}S_\varepsilon(\tilde y\mid c_b)\Big),\qquad \operatorname{cen}(S)=\tfrac12(\inf S+\sup S),$$
the **center + radius** of the section of the **union-fiber family** $\big(\mathcal F_{\mathrm{union}}\big)_{c_b}$, where auxiliary information enters as a **fiber** (Part I), archive noise as a **union** $\mathcal F_{\mathrm{union}}=\bigcup_{W\in\mathcal W_{\mathrm{arch}}}\mathcal F_W$ (Part III), and irregular patterns via the identification conditions determining $\mathcal W_{\mathrm{arch}}$ (Part II).

Canonicity (frozen): minimax-optimal by composition (CI-D2); the *unique conditionally-minimax* rule on realizable data when $\omega(2\varepsilon)<\infty$ (OP-3); a **partial** operator (three flag sources, §8); radius valued in the compactified $[0,+\infty]$; **nonlinear** in general — for linear families first at $k=3$ (the two-point midrange is a $k=2$ coincidence; OP-4b). No unary form $A(S,x)$ exists (OP-10). The family argument is minimally-sufficiently the trace set.

---

## 5. Exact definitions of all inputs

| Input | Definition | Origin |
|---|---|---|
| **archive** | $\{(g_a,D_a,\delta_a,c_a)\}$: prior members, finite designs, per-member value tolerances, auxiliary labels; induces $\mathcal W_{\mathrm{arch}}=\{$declared-class window systems $(\delta_a)$-consistent with the (labeled) archive$\}$ | CI-C1, CI-D1 |
| **$c_b\in C$** | known auxiliary of the new member; acts by fiber restriction $\mathcal F\to\mathcal F_{c_b}=\{f:(c_b,f)\in\mathcal F^+\}$; $C$ arbitrary (set membership only) | CI-A1 |
| **$S_b$** | the finite labeled support: design $D=\{x_1,\dots,x_k\}$, $k\le5$; data $\tilde y\in\mathbb R^k$; noise level $\varepsilon\ge0$ (closed ball); enters *only* through the section | SA1–SA4, IB-2 |
| **$x$** | query point $\notin D$; selects window and validity | F18, MP-6 |
| **$\varepsilon$** | known noise level; an **explicit argument** (no single rule optimal across levels, F1 Rem. 1.4) | A4 |

Derived objects: section $S_\varepsilon(\tilde y\mid c)=\{t:\exists(u,t)\in T_{D,x}(c),\ \|u-\tilde y\|_\infty\le\varepsilon\}$; conditional modulus $\omega_{x,D}(t\mid c)$; radius of information $\rho=\tfrac12\omega(2\varepsilon\mid c)$.

---

## 6. The information ledger [P]

| Source | Supplies | Certificate |
|---|---|---|
| **Archive** | the consistent set $\mathcal W_{\mathrm{arch}}$ of augmented, **size-$\le k{+}1$-truncated**, **covered-region** window systems — identified under Part II's pattern conditions, to residual ambiguity $h$ under Part III | CI-B1/B5, CI-C1, DM-3 |
| **Auxiliary $c_b$** | fiber selection within each candidate ($\mathcal F_W\to(\mathcal F_W)_{c_b}$) **plus** misspecification detection (empty fiber) | CI-A1, CI-A3 |
| **Support $S_b$** | the section cut — **at most $k$ continuous dimensions** of member identity ($2k$ inequalities) | F20/CP-3, DM-3 |
| **Query $x$** | selection of the window and the validity flag (covered? in the identifiable region of each $W$?) | F18, O4 |
| **Irreducible ambiguity** | exactly $\tfrac12\,\omega^{\mathrm{union}}_{x,D}(2\varepsilon\mid c_b)$ | Theorem 1 on the union-fiber family |

The five reductions compose without cross-terms under §8 (union commutes with fiber; conditioning commutes with sectioning; the new label does not feed back into archive consistency).

---

## 7. The mathematical objects a future system must approximate

Four objects (Phase 4/5), plus the Phase-6 union-fiber structure they now range over:

- **O1 — class representation** $\Phi(\mathbb T)$: the **closed, projectively consistent, size-$(k{+}1)$-truncated window system on the covered region** — now the *consistent set* of such systems, $\mathcal W_{\mathrm{arch}}$, in the augmented (labeled) class. Operationally: **conservative outer envelope pairs** $(\inf S,\sup S)$ per configuration. Latent coordinates carry no meaning (gauge); only the decoder image does.
- **O2 — task adaptation** $U(z,S_b,\varepsilon)$: the section-determining intermediate (labeled sample, $\varepsilon$-aware per DM-5 strata, flag-raising on empty sections).
- **O3 — query readout** $R(z_S,x)$: center + compactified one-sided radius; query-coupled.
- **O4 — coverage/validity/realizability decision map**: one-sided (false inclusion = fabrication; false exclusion = conservative), definable-flag-describable for tame classes.

**What must NOT be approximated away:** the certificate, the partiality flags, the discontinuities, the $\varepsilon$-dependence, the capacity ceilings, the gauge-invariance of every reported quantity, and — added at Phase 6 — the fiber (auxiliary) channel and the union (archive-noise) enclosure.

---

## 8. Necessary constraints inherited from Phases 4–6

The composite must satisfy (each a theorem; violation is certified error):

- **Invariances:** permutation symmetry of the support (composite-level theorem, factor-level WLOG — IB-1/AP-1); joint affine equivariance including reflections, $\varepsilon\mapsto|\alpha|\varepsilon$ (IB-5/AP-4); gauge invariance — no reported quantity may depend on latent coordinates or basis choice (IB-4/CR-5).
- **Uncertainty:** output the pair (center, radius $\in[0,+\infty]$); **outer** (one-sided) envelopes — under-approximation is a false certificate; radii **nondecreasing in $\varepsilon$** (A4) and in each $\delta_a$ (CI-C4), **nonincreasing under support refinement** (A5); **the center is exempt from all monotonicity** (OP-8).
- **Regularity:** no imposed continuity across certified transitions (MP-4/CR-3 — soft/randomized selection does not escape the half-jump); no monotone-center updates; center $1$-Lipschitz in the reported envelopes (IB-10); reproduction on exact identifiable data (IB-12); sensitivity-sum on constants-including classes (IB-11).
- **Structure:** the branch normal form — finitely many $C^r$ branches on definable cells jointly in class/support/query/$\varepsilon$, discrete selector, margin-guarded composition (CR-1/CR-4); flags exhaust the undefined locus (CR-8); certificate inflation is the lawful alternative to discreteness (CR-6).
- **Declaration:** the model **must declare** its closure class, its tameness route (definable or closed-form), and its stability certificate (modulus locally bounded) — *no closure assumption ⇒ no valid radius* (M1).
- **Phase-6:** compute $\mathcal W_{\mathrm{arch}}$ in the **augmented (labeled) class**; treat auxiliary conditioning as fiber restriction with the empty-fiber detection branch; report the **union** enclosure under archive noise; carry the DM-5 stratum membership as a flagged piecewise object.

Full lists: `future_model_constraints.md` M1–M12; `adaptation_operator_properties.md` AP-1…AP-8; `combined_conditional_operator.md` ledger.

---

## 9. Capacity and identifiability limits

- **Support capacity. [P]** At most $k$ continuous dimensions of member identity from $k$ observations (F20/CP-3); $k\le5$ is the standing budget. Exact identifiability at $(D,x)$ iff $\omega_{x,D}(0)=0$ (F2); stability iff $\omega(0^+)=0$ (F3) — strictly stronger (tanh, C3).
- **Task dimension at $k\le5$. [P/C]** Exact stable recovery only for $d\le5$; generic-design global identification of nonlinear analytic families only for $d\le2$ ($2d{+}1$ threshold, sharp at $d=1$ via $\sin\theta x$, C4); nonlinear $d\in\{3,4,5\}$ requires per-configuration verification. Linear families: unisolvent designs of size $d$ exist on any set; constrained families can identify below full rank (F6(iii)).
- **Design-independence. [X for $d\ge2$ continuous]** On triod-containing domains no size-$d$ design is universally unisolvent for a $d\ge2$-dim continuous family (Mairhuber–Curtis–Sieklucki, F10) — identifiability is jointly a property of family and design.
- **Latent dimension (meta). [P]** For the exactly-$d$ class over $N$ covered points ($\cong\operatorname{Gr}(d,N)$): floor $d(N-d)$ in continuous **and** definable categories; Whitney ceiling $2d(N-d)$; topological excess where obstructed ($\operatorname{Gr}(1,3)=\mathbb{RP}^2$ needs $4>2$, full-class hypothesis); no finite continuous latent for infinite-dimensional classes; **no globally continuous basis-style gauge** (Stiefel–Whitney, CR-5). Minimal member summary is query-relative, $\dim\le\operatorname{rank}G$ (CP-3). **Budget truncation:** only size-$\le(k{+}1)$ windows are learnable (DM-3).
- **Archive counting. [P]** Irregular pattern identifies $V|_U$ only if $\sum_a(k_a-d)_+\ge d(N-d)$ (CI-B1) and every point has $\ge d$ observers (CI-B2, non-pivot genericity); connectivity is **not** sufficient (CI-B3).

---

## 10. Approximation and stability limits

- **Master decomposition. [P]** prediction error $\le\ \tfrac12\,\omega(2\varepsilon+2h)+h\ +\ \eta_s\ (+\ \tfrac12 J$ at transitions for continuous models$)$: representation error $h$ enters **inside the modulus** (nonlinear, possibly infinite — C3); selection slack constant $1$ (tight, MP-5); half-jump sharp and localized (MP-4).
- **Entropy floor. [P]** $N(\mathcal C,2h)\lesssim(Lr/h)^m$ — impossibility bounds the **pair** $(m,L)$, never dimension alone (CR-9/DM-8).
- **Positive convergence. [C]** Declared tame + stable + compact-cell scope ⇒ finite-parameter $C^1$-convergent certified approximation via branch form + Bernstein density + outer rounding (CR-6); effective for semialgebraic classes.
- **Jump/instability impossibilities. [X]** No uniform continuous convergence across genuine jumps (MP-4/CR-3). Instability trichotomy (CR-7): global instability forces a never-vanishing void band; **local instability at realizable data forces radius $+\infty$ under any representation error $h>0$** (unconditional, witness-backed); exactly-representable unstable classes are honest but vacuous.
- **Sample/task floors. [P]** Excess-information floor $\sum(k_a-d)_+\ge d(N-d)=\dim\operatorname{Gr}(d,N)$ (DM-7, core pattern; CI-B1 general necessity); per-task support $k\ge d$ ($2d{+}1$ nonlinear-generic).

---

## 11. Archive irregularity and archive-noise results

**Irregular archives (Part II). [P/C/O]** The object is a **sheaf-like consistency system**: local windows with a gluing law.
- Gluing across **unisolvent** overlaps is exact (CI-B4a). **[P]**
- Rank-$<d$ overlaps produce **holonomy** — local consistency without global sectionability (explicit $d=2$ witness; holonomy iff transported partial identifications fail to span the member space around a cycle). **[P]**
- Sufficient: connected patches covering $U$, each locally F17-identified to full dimension $d$, all overlaps unisolvent ⇒ identified, cycles automatic (CI-B5). **[C]**
- Size-$(d{+}1)$ windows determine $V|_U$ (CI-B6, DM-3 coherence). **[P]** $d=1$: identifiable iff the **nonzero-value** incidence subgraph is connected (CI-B7). **[P]**
- Exact necessary-and-sufficient combinatorial characterization of **unique** completability for general $d$: **[O non-blocking]** — CI-B1 (necessity) and CI-B5 (sufficiency) sandwich it; literature (Pimentel-Alarcón–Boston–Nowak) characterizes only *finite* completability.

**Imperfect archives (Part III). [P/C/X/O]**
- **Union reduction (CI-C1). [P]** The hull-of-union operator is minimax-optimal with **valid outer certificates at any $\delta$** — **the OPEN archive-noise constant is bypassed for validity**, demoted to a size question.
- **Perturbation (CI-C2). [C]** $h\le C\delta$ under $\sigma_0$-conditioning and $\delta\le\delta_0(\sigma_0,M,d,\text{transport})$; $C\sim C_{\mathrm{step}}^{L}$ (exponential in transport depth); then query inflation $\le\tfrac12\omega^{\mathrm{surr}}(2\varepsilon+2h)+h$.
- **No uniform bound (CI-C3). [X]** Value-diameter swings $\Theta(\delta/\sigma_{\min})$ and $\kappa^L$ chains preclude any conditioning-free bound; the missing regularity is exactly the $\sigma_0$ floor.
- Radii nondecreasing in $\delta$ (CI-C4). **[P]** Sharp exact constant $C$: **[O non-blocking]**.

---

## 12. Auxiliary-information results

- **Useful-iff (CI-A3). [P]** Auxiliary information is useful at $(D,x)$ **iff it changes the joint window** $T_{D\cup\{x\}}$. (Proof uses the $\varepsilon=0$ clause; covers pure misspecification-detection.)
- **Conditional ambiguity reduction. [P]** $\omega_{x,D}(t\mid c)\le\omega_{x,D}(t)$ always; the certified information supplied by $c_b$ is $\Gamma_c=\tfrac12(\omega_{x,D}(2\varepsilon)-\omega_{x,D}(2\varepsilon\mid c_b))\ge0$, zero iff useless in the worst case. Synergy witness: $c$ + one evaluation identify a value neither identifies alone (CI-A, determination 3).
- **Incorrect-auxiliary substitution (CI-A5).** Conditioning on a wrong $c'$ either **fires the realizability flag** (= emptiness of $S_\varepsilon(\tilde y\mid c')$) or errs silently; the silent error is $\le\omega_{x,D}(2\varepsilon\mid c',c_b)$ **[P]**. Harmless-for-all-data iff $T_{D\cup\{x\}}(c')=T_{D\cup\{x\}}(c_b)$ **[P]**. **No distance on $C$ grades the harm without declared structure on $C$** **[X]**.
- Reduction scope: the fiber view is primitive; augmentation preserves the exactly-$d$ class **iff** $c=L(f)$ is linear **[C]**.

---

## 13. Decisive falsification tests

All are worst-case adversarial protocols with computable constants; benchmark averages neither confirm nor refute (F3).

**Base + meta (P1–P10):** P1 radius floor · P2 capacity ceiling ($k$ dims) · P3 sensitivity sum · P4 query-degradation curves · P5 support-placement ratios (exact small-$k$ Lebesgue constants) · P6 off-coverage collapse · P7 the $2d{+}1$ cliff ($\sin\theta x$) · P8 noise fragility (tanh) · P9 discontinuity spikes · P10 non-monotone updates.

**Differentiable/meta-object (NP-1–NP-6):** NP-1 budget truncation (perturb size-$(k{+}2)$ windows) · NP-2 excess-information floor · NP-3 meta-transition spikes (archive-side) · NP-4 $\varepsilon$-sweep monotonicity · NP-5 partiality surfacing (engineered inconsistent support) · NP-6 comparability (induced-operator metric only; no latent-space comparison).

**Conditional/archive (T1–T5):** T1 replace $c_b$ by unrelated $c'$ → realizability-flag rate / interval shift · T2 remove $c_b$ → radius increase by $\Gamma_c$ · T3 foreign support → flag or center shift · T4 break overlap structure → residual-ambiguity dimension (radius on affected points) · T5 increase $\delta$ → radius nondecreasing ($O(\delta)$ conditioned; $\Theta(\delta/\sigma_{\min})$ / $\kappa^L$ near degeneracy).

Each test names a certificate quantity forced to change **iff** its source is genuinely informative — decisive in both directions.

---

## 14. Epistemic status of the whole program

- **Proved [P]:** Theorem 1 and the minimax framework; F2/F3 identifiability–stability boundaries; F6 rank; F7 linear minimax; F10 design-independence obstruction; F15/F16 relative values & gauge; F17/F18 archive; F19/F20 summary dimension; CP-1/2/3; OP-1…OP-10; MP-1/2/4/6; DM-1/2/3/6/7 (and DM-8 under its hypotheses); CR-1/2/3/4/5/6 and the CR-7 instability facts; CI-A1/A2/A3; CI-B1/B3/B4a/B4b/B6/B7; CI-C1/C4; CI-D1/D2; T1–T5.
- **Conditional [C]:** F12/F13 generic identifiability (bounded, connected, real-analytic separating); MP-2/MP-3 (bounded prior sections / gauge anchor); DM-5 $\varepsilon$-strata; CR-6 positive convergence (tame + stable + effective); CI-A augmentation transfer (linear label); CI-B2 (non-pivot), CI-B5 (cover + local F17 + unisolvent overlaps); CI-C2 ($\sigma_0$-conditioning, $\delta\le\delta_0$); CI-D2 (well-specification, no-coupling, augmented-class $\mathcal W_{\mathrm{arch}}$).
- **Impossible [X]:** off-design determination without constraints (F4); design-independence for $d\ge2$ continuous (F10); no finite continuous summary for infinite-dimensional classes (F19); no unary operator (OP-10); discontinuity un-smoothable (MP-4/CR-3); local-instability certificate void (CR-7); no uniform archive-noise bound without conditioning (CI-C3); no distance-graded auxiliary harm without structure on $C$ (CI-A5).
- **Open but non-blocking [O]:** sharp archive-noise constant $C$ (validity bypassed by CI-C1); exact unique-completability characterization for general $d$ (sandwiched by CI-B1/B5); exact optimality under member/archive coupling (validity holds; linear class uncoupled); joint-loss constants, $\varepsilon$-adaptive rules, general off-coverage dichotomy, anchored-projection exact class (all scope-discharged). Distributional guarantees are **out of scope by design**, not open.

**No open item blocks construction, certification, or falsification of the operator.**

---

## 15. Final stopping verdict

$$\boxed{\textbf{THEORY\_COMPLETE\_FOR\_HANDOFF}}$$

The program that began with *"with no assumptions the data determine nothing off the design"* ends by absorbing its three most realistic complications — auxiliary information, irregular archives, imperfect archives — into the same single primitive via three reductions to one minimax identity: conditioning (fiber), gluing (sheaf), union. The operator, its certificate, its ledger, and its falsification tests survive each relaxation intact.

---

## 16. MATHEMATICAL INTERFACE FOR THE NEXT RESEARCHER

**What must be approximated.**
The restriction of the canonical conditional operator $A_c$ to a **declared, tame, stability-certified class**, as four objects: **O1** the closed, size-$(k{+}1)$-truncated, covered-region window-system *consistent set* $\mathcal W_{\mathrm{arch}}$ in the augmented (labeled) class, surfaced as conservative **outer** envelope pairs; **O2** the section-determining task-adaptation map; **O3** the query-coupled readout emitting (center, radius $\in[0,+\infty]$); **O4** the coverage/validity/realizability decision map. The output at each query is the center of the union-fiber section and its radius.

**What information each input must contribute.**
- *Archive* → the consistent set of window systems on the covered region, to its residual ambiguity $h$ (identified only where Part II's pattern conditions hold; honestly widened where they do not).
- *Auxiliary $c_b$* → fiber selection within each candidate, plus misspecification detection via the empty fiber.
- *Support $S_b$* → the section cut, contributing at most $k$ continuous dimensions of member identity; must carry locations, values, and $\varepsilon$.
- *Query $x$* → window selection and the validity flag.
- *Nothing else may contribute anything*: any apparent information beyond these is either fabrication (off-coverage, unbounded section) or an undeclared class assumption.

**What properties the approximation must preserve.**
Permutation symmetry, joint affine equivariance (including reflections), gauge invariance of all reported quantities; the (center, radius) pair with outer one-sided semantics and compactified radius; radii monotone in $\varepsilon$ and $\delta$, nonincreasing under support refinement, with the center exempt; reproduction on exact identifiable data; $1$-Lipschitz center-in-envelopes and the sensitivity-sum; the branch/flag/margin structure across certified discontinuities; the declared closure class, tameness route, and stability certificate; the capacity ceilings ($k$ task dimensions, size-$(k{+}1)$ windows) and the latent-dimension floors.

**What behaviors would falsify the intended mechanism.**
Any violation of P1–P10, NP-1–NP-6, or T1–T5; specifically: beating the radius floor $\tfrac12\omega(2\varepsilon)$ in the worst case; recovering more than $k$ task dimensions; emitting finite radii off coverage or at unstable-local configurations; certificate under-coverage (false outer enclosure); $\varepsilon$-insensitive radii; attaching meaning to latent coordinates or comparing latents across systems; a radius unchanged when $c_b$ is removed yet the auxiliary is genuinely window-changing (T2); no flag when a wrong $c'$ or foreign support is substituted (T1/T3); ambiguity unchanged when overlap structure is broken at a point below $d$ observers (T4); radius insensitive to genuinely binding $\delta$ (T5).

*This document states what must be approximated, what each input must contribute, what must be preserved, and what would falsify the mechanism. It states nothing about how any of these is to be implemented. It freezes the complete theoretical program.*
