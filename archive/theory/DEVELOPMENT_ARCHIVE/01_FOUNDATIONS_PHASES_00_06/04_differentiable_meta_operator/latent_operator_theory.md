# Latent Operator Theory — Finite Representation of Window Systems

> **Status:** Phase-4 derivation, 2026-08-02. Sources: the frozen corpus (F-, C-, CP-, OP-, MP-numbers; `../RESEARCH_HANDOFF.md`). New results carry DM-numbers; all were adversarially refereed (verdicts and corrections incorporated — notably the operational-topology requirement, category labels on every dimension bound, and the full-class hypothesis on the projective-plane obstruction). No neural embeddings are assumed anywhere; "representation" is a pure mathematical object.

**Question.** How can an infinite window system $\mathbb T=(T_S)_{S\ \text{finite}\subset U}$ be represented by a finite-dimensional learnable object $z_T$?

---

## 1. The topology must come first (DM-1)

Dimension theory is meaningless until the space of window systems carries a topology — and the naive choice fails: windows of linear families are unbounded sets, so Hausdorff distance between distinct subspaces is $+\infty$, inducing the **discrete** topology, under which every injection is continuous and every lower bound evaporates.

**Definition (operational topology).** Topologize a class $\mathcal C$ of window systems by the **initial topology of the envelope pairs**: the coarsest topology making $(\mathbb T,D,\tilde y,x,\varepsilon)\mapsto(\inf S_\varepsilon(\tilde y),\sup S_\varepsilon(\tilde y))$ continuous in $\mathbb T$ at every realizable covered configuration. This is the topology *operationally relevant to learning*: two systems are close iff they induce nearby operators.

**Lemma DM-1 (identification for the exactly-$d$ linear class; refereed).** Let $U$ be finite, $|U|=N\ge d$, and $\mathcal C_d$ = window systems of exactly-$d$-dimensional subspaces $V\subseteq\mathbb R^U$. The map $V\mapsto\mathbb T(V)$ is a bijection onto $\mathcal C_d$ (the window at $U$ *is* $V$), and via orthogonal projection matrices $V\mapsto P_V$, $\mathcal C_d$ with the operational topology is homeomorphic to the Grassmannian $\operatorname{Gr}(d,N)$ with its usual (principal-angle) topology. *Proof sketch:* the envelope pair is continuous in $P_V$ (suprema over slabs of continuously varying polytopes), giving continuity of $\mathrm{id}:(\operatorname{Gr},\text{gap})\to(\mathcal C_d,\text{operational})$; the operational topology is Hausdorff because envelope data separate distinct subspaces; a continuous bijection from a compact space to a Hausdorff space is a homeomorphism. $\square$

All statements below carry a **category label**: (cont) = continuous category; (def) = definable-but-possibly-discontinuous category; (meas) = merely measurable. This is forced by the corpus's own tension: MP-4 forbids demanding continuity of the operator, while C8 (Borel collapse) kills dimension claims without *some* regularity. The o-minimal category is the reconciliation: **definable injections preserve dimension** even when discontinuous, so definable-category floors survive where continuous-category obstructions do not.

---

## 2. What a representation is

**Definition.** A *representation* of a class $\mathcal C$ is a pair $(z,\beta)$: an encoder $z:\mathcal C\to\mathbb R^m$ and a decoder $\beta:\mathbb R^m\to(\text{closed, projectively consistent window systems})$ with $\beta(z(\mathbb T))=\mathbb T$ (exact) or within error $h$ in the operational topology (approximate). Exactness has **two** directions (refereed correction): injectivity of $z$ modulo decoder-gauge, *and* decoder-range surjectivity onto the class quotient — a decoder that cannot express some member of the declared class is inexact no matter how good the encoder.

Three standing constraints on $\beta$'s range, inherited from the corpus:
1. **Closure:** at $\varepsilon>0$, window systems are data-determined only up to closure (CP-2); exactness means exactness of the *closed* system.
2. **Consistency and realizability:** the range must consist of projectively consistent systems, and consistency is necessary but **not sufficient** for realizability by a family (MP-1; Waterhouse) — countable coverage or compact windows restore sufficiency.
3. **Outer admissibility (certificate direction):** if the decoder range *inner*-approximates the declared class, reported radii under-cover and every certificate is false — the corpus's inner/outer tension recurring one level up. Admissible decoders outer-approximate the declared closure class, or radii must be evaluated against the declared class rather than the learned range.

---

## 3. Dimension theory of the meta-latent (DM-2)

**Theorem DM-2 (refereed: confirmed).** For the exactly-$d$ linear class $\mathcal C_d\cong\operatorname{Gr}(d,N)$, $\dim\operatorname{Gr}(d,N)=d(N-d)$:

(a) **(cont) Floor.** Any continuous injective $z:\mathcal C_d\to\mathbb R^m$ requires $m\ge d(N-d)$ (chart + invariance of domain). **The meta-latent dimension must grow with coverage**: linearly in the number of covered locations.

(b) **(def) The floor survives discontinuity.** In the definable category the same floor holds — definable injections preserve o-minimal dimension — so escaping to discontinuous-but-tame encoders does not reduce $m$ below $d(N-d)$. Only (meas) encoders collapse dimension (C8), and they are excluded by any stability requirement.

(c) **(cont) Ceiling.** $m=2d(N-d)$ always suffices (strong Whitney embedding of compact smooth manifolds). A *constructive, gauge-free, equivariant* alternative: the projector embedding $V\mapsto P_V\in\operatorname{Sym}(N)$, $m=\tfrac{N(N+1)}2$ — larger than Whitney's bound, exposing a genuine derived tension: **equivariant representations may require strictly more dimensions than non-equivariant ones** (recorded as an inductive-bias tradeoff in `meta_learning_inductive_bias.md`).

(d) **(cont) The floor is not attained in general — the meta-tripod.** $\operatorname{Gr}(1,3)\cong\mathbb{RP}^2$: dimension $2$, embeds in $\mathbb R^4$ (Veronese-type map), admits **no continuous injection into $\mathbb R^3$** (Alexander duality: an embedded closed surface in $S^3$ separates, hence is orientable; compactness upgrades no-embedding to no-continuous-injection). So for $d=1,N=3$ the minimal continuous latent dimension is exactly $4>2$. *Hypothesis (refereed):* this requires the class to be **all** of $\operatorname{Gr}(1,3)$, including the degenerate coordinate-vanishing lines (the C13-type configurations); excluding them leaves an open subsurface (a Möbius band) which does embed in $\mathbb R^3$. In the (def) category the obstruction disappears entirely (cut along a curve); it is a continuous-category refinement.

(e) **Unbounded classes.** For classes whose window systems are genuinely infinite-dimensional — e.g. containing continuously injected copies of $\operatorname{Gr}(d,N)$ for every $N$ (window-wise topology; refereed wording) — no finite $m$ suffices in (cont) or (def). *Scope (refereed correction):* this does **not** apply merely because $|U|=\infty$: exactly-$d$ subspaces of a fixed finite-dimensional ambient span form a finite-dimensional Grassmannian regardless of coverage. The obstruction is infinite-dimensionality of the class, not infinitude of the domain.

---

## 4. What is representable is bounded by what is learnable (DM-3)

**Theorem DM-3 (budget truncation — new; refereed via plan audit).** At support budget $k$ (and scalar queries), the canonical operator depends on $\mathbb T$ only through the joint windows of size $\le k+1$. Projective consistency determines small windows from large but **not conversely**; two window systems agreeing on all windows of size $\le k+1$ induce identical operators at budget $k$, while their larger windows may differ (already at $N=3$, $k+1=2$: the three pairwise windows do not determine the triple window). **Consequently the learnable meta-object is exactly the closed size-$(k+1)$ truncation of the window system on the covered region** — a hard ceiling on every training signal, independent of model capacity. (Falsification protocol NP-1, `model_design_interface.md`: perturb windows of size $k+2$ leaving all smaller windows fixed; the operator, and hence any legitimate training signal, is unchanged.)

Representation targets should therefore be truncated systems; encoding structure beyond size $k+1$ is encoding the unlearnable.

---

## 5. Information loss

**Forward direction (established; nonlinear, possibly infinite).** A representation error $h$ (operational topology) inflates sections into $S_{\varepsilon+h}\pm h$, giving prediction error $\le\tfrac12\omega_{x,D}(2\varepsilon+2h)+h$. The composition is **not additive** — the window error enters *inside* the modulus — and for identifiable-but-unstable families (C3) it is infinite for every $h>0$. Certificate semantics survive only under one-sided (outer) decoding (§2.3).

**Converse direction (DM-8 — entropy floor; stated under its hypotheses).** If the decoder is $L$-Lipschitz on a bounded latent domain of dimension $m$, then the class is covered within error $h$ by the image of an $h/L$-net of the latent ball, so $N(\mathcal C,2h)\lesssim(Lr/h)^m$; classes whose metric entropy exceeds every such bound force $h$ bounded away from $0$ for every finite $m$. For Lipschitz-type classes (entropy of infinite-dimensional balls) this recovers the no-finite-representation theorem *quantitatively*: information loss is bounded below by entropy numbers, making "information loss" two-sided. (Hypotheses — Lipschitz decoder, bounded domain — are load-bearing; without them C8-type collapse reappears.)

---

## 6. Identifiability of the latent

1. **Second-level gauge.** If $\beta\circ g=\beta$ for a transformation $g$ of latent space, no data distinguish $z$ from $g(z)$: the latent is identifiable at most modulo the decoder-stabilizer. **Latent coordinates are meaningless**; only $\beta(z)$ — the induced window system — carries identifiable content. (The base-level theorem CP §4.5, lifted.)
2. **Only through covered, truncated, closed windows.** The archive pins $\beta(z)$ exactly on the covered region under F17-type conditions, only up to closure at $\varepsilon>0$ (CP-2), and only up to size-$(k+1)$ truncation (DM-3). The identifiable content of the latent is precisely the **closed, budget-truncated, covered-window quotient**.
3. **Comparability (DM-9).** The induced distance between two representations, $\sup$ over realizable covered configurations of the Hausdorff distance between reported envelope intervals, is a pseudometric descending to a metric on the double quotient (decoder-gauge $\times$ window-gauge). Cross-representation comparison is legitimate only in this metric; any latent-space metric is non-identifiable and carries no meaning (formalized as protocol NP-6).

---

## 7. Summary: necessary properties, sufficient conditions

| Requirement on $(z,\beta)$ | Status | Source |
|---|---|---|
| $m\ge\dim$ of the class (cont **and** def categories) | necessary | DM-2(a),(b) |
| $m\ge$ embedding dimension when topology obstructs (full-class hypothesis) | necessary (cont) | DM-2(d) |
| $m=2\cdot\dim$ | sufficient (cont, compact smooth classes) | DM-2(c) |
| decoder range: closed, consistent, realizable, outer-admissible | necessary for valid certificates | §2 |
| target = size-$(k+1)$ truncated covered quotient | exactly the learnable content | DM-3, §6.2 |
| finite $m$ for infinite-dimensional classes | impossible (cont/def); entropy-quantified | DM-2(e), DM-8 |
| latent interpretability | none (gauge); only $\beta(z)$ identifiable | §6.1 |

**Degenerate cases.** $N\le d$: $\operatorname{Gr}$ is a point, $m=0$ suffices — the class is rigid and the archive is informationless but also unnecessary. $d=0$: the zero family; everything forced. Classes excluding degenerate subspaces: dimension floors persist, topological obstructions may not (DM-2(d) hypothesis).
