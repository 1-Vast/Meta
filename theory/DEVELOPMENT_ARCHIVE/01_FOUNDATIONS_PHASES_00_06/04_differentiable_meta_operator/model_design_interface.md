# Model Design Interface — Handoff for the Architecture Phase

> **Status:** Phase-4 handoff, 2026-08-02. Audience: a researcher who will design a deep meta-learning system. This document answers **"what mathematical objects must a future architecture approximate?"** — it is an interface specification, not an architecture. Sources: the frozen corpus and the Phase-4 files (`latent_operator_theory.md`, `adaptation_operator_parameterization.md`, `learnability_conditions.md`, `meta_learning_inductive_bias.md`; DM-numbers, IB-numbers). Everything below is refereed mathematics; nothing below chooses a module, a loss implementation, or a domain.

---

## 1. The four objects to approximate

A future system approximates the canonical operator $\mathbb A(\mathbb T;S,x,\varepsilon)=(\text{center},\text{radius})$ through **four** mathematical objects (the fourth is not optional — refereed plan audit):

**O1 — The class representation** $z=\Phi(\mathbb T)$.
Target: the **closed, projectively consistent, size-$(k{+}1)$-truncated window system of the declared class on the covered region** (DM-3: nothing finer is learnable). Constraints: latent dimension floors $\dim$-of-class in *both* the continuous and definable categories, with topological excess where obstructed (DM-2; the $\operatorname{Gr}(1,3)=\mathbb{RP}^2$ example: floor 2, minimum 4, full-class hypothesis); constructive equivariant option (projector embedding) at higher dimension — the equivariance/width tradeoff is a fact to budget, not to wish away; decoder range must be closed, consistent, realizable (MP-1 caveats), and **outer-admissible** relative to the declared class (inner-approximating decoders make every certificate false). Latent coordinates carry no meaning (gauge); only the decoder image does.

**O2 — The task-adaptation map** $z_S=U(z,S,\varepsilon)$.
Target: the section-determining intermediate. Constraints: labeled-sample dependence (locations *and* values — IB-2); symmetric WLOG (IB-1); width bounded below by the query-relative rank $r_Q$ and *not* matching it in general (DM-4: the Lipschitz $\Theta(k)$ witness; the F19 sandwich governs); $\varepsilon$ **enters $U$ in general** — deferral to the readout is valid only in declared strata (DM-5: surjective-trace linear; Lipschitz — and *fails already for overdetermined linear designs*, which are typical at $k=5$, $d\le2$); carries data-dependent certificate content outside stratum (i); emits the misspecification flag on empty sections rather than a state (IB-8).

**O3 — The query-conditional readout** $(\hat y,\hat\rho)=R(z_S,x)$.
Target: center and radius of the section. Constraints: query coupling (validity region, $x$-dependent sensitivity profile, query-dependent certificate — IB-3); radius valued in the **compactified** $[0,+\infty]$ with $+\infty$ a required output; one-sided (outer) envelope semantics; selection coherence (center $1$-Lipschitz in reported envelopes — IB-10); sensitivity-sum on constants-including classes (IB-11); reproduction on exact identifiable data (IB-12); no continuity assumption across certified transitions (IB-9); no monotone-update structure (IB-7/OP-8).

**O4 — The coverage/validity decision map (DM-10).**
Target: the flags "is $(D\cup\{x\})$ within the covered region?", "is the query in the validity region?", "is the support realizable at level $\varepsilon$?". This object has its **own one-sided semantics** — a false *inclusion* is fabrication (F18/A10); a false *exclusion* is merely conservative — and its own discontinuity structure: at the coverage boundary the radius jumps to $+\infty$, where continuous-surrogate error bounds become vacuous unless the codomain is compactified. For tame classes the flag sets are definable and finitely describable (DM-6(d)). No architecture may fold O4 silently into O3: the decision is a separate mathematical object with a separate error semantics.

---

## 2. The training-signal statement

What any training procedure can pin, at most:

$$\textbf{the closed, size-}(k{+}1)\textbf{-truncated window quotient on the covered region, modulo decoder gauge.}$$

- Nothing off the covered region (F18); nothing beyond the truncation (DM-3); nothing about latent coordinates (gauge); nothing below closure at $\varepsilon>0$ (CP-2).
- The archive requirement to pin even this much, in the exactly-$d$ linear class under the core pattern: $n\ge d$ diverse tasks and clipped excess observations $\sum_j(k_j-d)_+\ge d(N-d)=\dim\operatorname{Gr}(d,N)$ — **only observations beyond a task's own degrees of freedom carry meta-information** (DM-7; beyond-core patterns OPEN).
- **Correctness metric (declared):** sup over realizable covered configurations, on compact subsets of open cells, of the error against the canonical operator — with envelope slack measured one-sidedly. No average-case surrogate certifies it (F3).

---

## 3. Error budget

| Source | Transfer to prediction error | Status |
|---|---|---|
| class-representation error $h$ (operational topology) | $\tfrac12\,\omega(2\varepsilon+2h)+h$ — **inside the modulus; nonlinear; possibly $\infty$** (C3) | proven; composition, not a sum |
| envelope/selection slack $\eta_s$ | additive, constant $1$ (tight) | proven (MP-5) |
| continuous smoothing at transitions | $\ge\tfrac12\cdot$ local jump, localized on lower-dimensional definable sets | proven (MP-4, DM-6) |
| entropy floor of an $m$-dimensional Lipschitz decoder | $h\gtrsim$ entropy numbers of the class — two-sided information loss | proven under its hypotheses (DM-8) |
| archive noise | — | **OPEN** (handoff §9.1); no additive bound may be assumed |

---

## 4. Falsification hooks

Inherited: **P1–P10** (radius floor; capacity ceiling; sensitivity sum; query-degradation curves; support-placement ratios; off-coverage collapse; the $2d{+}1$ cliff; noise fragility; discontinuity spikes; non-monotone updates). New at this phase:

- **NP-1 (budget truncation; DM-3).** Perturb windows of size $k+2$ leaving all smaller windows fixed: the operator — hence any legitimate training signal — is unchanged. A model whose behavior differs is representing its prior, not the data.
- **NP-2 (excess information; DM-7).** Archives of core-pattern tasks with clipped excess below $d(N-d)$: identification must fail with realized ambiguity; apparent success identifies an undeclared class restriction.
- **NP-3 (meta-transition spikes; distinguished from P9).** Vary the *archive* across a window-system transition at fixed support: continuous class-representations exhibit localized error spikes at the transition — P9's phenomenon one level up.
- **NP-4 ($\varepsilon$-sweep; IB-13.5).** Radii must be nondecreasing in the declared $\varepsilon$; an $\varepsilon$-free model is provably suboptimal at some noise level.
- **NP-5 (partiality surfacing; IB-13.6).** Engineered inconsistent supports must raise the flag, not a projected prediction.
- **NP-6 (comparability; DM-9).** Two trained systems are compared only in the induced-operator metric ($\sup$ Hausdorff distance between reported envelope intervals over realizable covered configurations — a metric on the double gauge quotient). Latent-space comparisons are inadmissible in evaluation.

---

## 5. The closing answer

**What must a future deep meta-learning architecture approximate?**

The restriction of the canonical operator to a **declared, tame class** — tame by uniform definability or by closed form (DM-6; the learnability line coincides with the identifiability line) — factored as the four objects above:

$$\Phi:\ \text{class}\to z\quad(\text{O1: truncated closed window quotient, outer-admissible decoder}),$$
$$U:\ (z,S,\varepsilon)\to z_S\quad(\text{O2: labeled, symmetric-WLOG, }\varepsilon\text{-aware, flag-raising}),$$
$$R:\ (z_S,x)\to(\hat y,\hat\rho)\quad(\text{O3: query-coupled center + compactified one-sided radius}),$$
$$V:\ (z,D,x,\tilde y,\varepsilon)\to\text{flags}\quad(\text{O4: coverage/validity/realizability, one-sided}),$$

under the constraints IB-1…IB-14, within the error budget of §3, trainable at most to the signal of §2, and falsifiable by P1–P10 and NP-1…NP-6.

What must it **not** approximate away: the certificate, the partiality, the discontinuities, the $\varepsilon$-dependence, the capacity ceilings, and the gauge-invariance of everything it reports. The mathematics is indifferent to how these objects are realized; it is exact about what they are.
