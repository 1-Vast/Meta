# Object Redesign (§1)

> **Status:** Phase-17 (meta-operator reconstruction), 2026-08-03. Phases 0–15 unmodified; Phase-15's parameterization is **superseded, not patched** — the four audited failures (`../16_final_deep_operator_audit/FINAL_VERDICT.md`: invalid cube parameters; constraint-breaking interpolation; accuracy-indexed sieve; non-closed continuous value class) are treated as symptoms of a wrong object choice, and the object is re-derived from the frozen theory. New results carry **MR-** numbers, tagged **[proved] / [declared] / [rejected]**. Role discipline: mathematics only — no networks, no models, no pipelines.

---

## 1. What the frozen theory forces the learnable object to be

Fixed by Phases 0–7 (not renegotiable): identification produces a set; decisions require a declared completion of the dominance order (DE-P); population information may only weight the inside of identified sets (DE-H2/H3); every emission carries certificates that the learned layer must not manufacture (DE-L5). Fixed by Phases 8–13: what the decision operator consumes from the learned layer is a **constraint class of laws with confidence and rung tags**, per (context, query, specification).

The key structural observation that dissolves the Phase-15 failures:

**Proposition MR-1 (the value space is a description polytope). [proved]**
Under a declared finite skeleton (atlas/grid; §2), every value the decision layer consumes is determined by a finite **band vector** $b$ (event endpoint pairs), and the set of *valid* band vectors
$$\mathbb B\ =\ \big\{\,b\ :\ \text{ranges, order }l\le u,\ \text{coherence ties, band monotonicity (Route B), and nonemptiness of }K(b)\,\big\}$$
is a **compact convex** subset of a finite-dimensional space — nonemptiness included: if $P_1\in K(b_1)$ and $P_2\in K(b_2)$ then $\lambda P_1+(1-\lambda)P_2\in K(\lambda b_1+(1-\lambda)b_2)$ (all constraints are affine in $(b,P)$ jointly), so feasibility survives convex combination; $\mathbb B$ is the projection of a lifted polytope in $(b,P)$-space. Full construction and the Route-B closure repair: `structure_preserving_parameterization.md`. $\square$

## 2. The candidate objects, decided

**(a) Point-valued operator (into laws). [rejected]** A single law per index fabricates: it erases the residual ambiguity that Phase 7 proved essential and DE-H2/L5 protect. Rejected on frozen grounds, independent of realizability.

**(b) Set-valued correspondence. [adopted — as semantics]** The decision layer *does* consume sets. But parameterizing arbitrary set-valued maps directly caused every Phase-15 failure: sets were encoded by unconstrained coordinates, and set-space has no linear structure for interpolation.

**(c) Constrained latent representation. [rejected]** Coordinates without validity predicates are meaningless or fabricating (MC-7, frozen gauge unidentifiability). Rejected — unchanged from Phase 9.

**(d) Probabilistic operator (a law over operators). [rejected]** Randomizing the operator adds an unidentifiable outer layer the decision operator cannot consume without collapsing it to its mean class; it supports no additional few-shot adaptation and inflates the object. Rejected as non-minimal.

**(e) The minimal object. [adopted]**
$$\boxed{\ \textbf{A point-valued map into the valid-description polytope }\mathbb B\ \text{— set-valued in semantics, convex-point-valued in representation.}\ }$$
The correspondence (b) is carried by the *denotation* $b\mapsto K(b)$; the learnable map $H\mapsto b$ is point-valued into a compact convex polytope. This is the weakest relaxation that keeps everything the frozen theory demands and makes realizability a convexity fact:

**Theorem MR-2 (why this object dissolves all four failures). [proved]**
(i) *Invalid parameters cannot exist:* parameters range over (products of) $\mathbb B$ itself — validity is the type, not a hoped-for property of a cube.
(ii) *Interpolation preserves constraints:* $\mathbb B$ is convex, so every convex combination of valid descriptions is valid — the exact operation Phase-15's multilinear interpolation performed *on the wrong space* (coordinate cubes) is safe on $\mathbb B$.
(iii) *No accuracy sieve inside a deployment:* the target of learning is a point of a fixed finite-dimensional compact set per index, and the canonical operator itself is included as a fixed decoder component (§2/§3) — the family need not grow with any $\varepsilon$.
(iv) *Closed feasible representation for continuous outputs:* the Route-B value class is redefined with the closed-constraint convention (lower bounds on closed intervals, upper bounds on open ones), making $K(b)$ $W_1$-closed — the audit's $\delta_{t+1/n}$ witness lands inside the set (proof in §3). $\square$

**Few-shot adaptation support [verified].** The object supports adaptation exactly as the frozen theory prescribes: the current task's identified set enters by likelihood-free support restriction of $K(b)$ (DE-H2); the context/query/specification indexing survives (LC-5/7); the confidence/rung side-channels are canonical observables. Nothing in the redesign touches the adaptation semantics — only the representation of the population object changed.
