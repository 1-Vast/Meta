# Learnability Conditions

> **Status:** Phase-4 derivation, 2026-08-02. Question: when can the abstract operator $\mathbb A$ be approximated by a finite-parameter model? All results refereed; scoping corrections incorporated (uniformly-definable-family hypothesis; clipped excess sums; common-core scope; nonlinear error composition).

**Scope declaration (required by the corpus).** Everything here is worst-case and distribution-free; correctness of an approximation means closeness to the canonical operator in the metric *sup over realizable covered configurations, on compact subsets of open cells* (see DM-6). No PAC or distributional sample-complexity statements appear — not as an omission but as a theorem-respecting choice: the corpus contains no distributional assumptions, and benchmark averages neither confirm nor refute any statement in it (failure-mode F3). Multi-query claims are per-query (sup-loss exactness only, F1 Rem. 1.1/1.2). The archive is assumed **exact**; the noisy-archive constant is OPEN (handoff §9.1) and no statement below silently assumes it.

---

## 1. The regularity theorem: definable classes are piecewise-tame (DM-6; refereed: confirmed)

**Hypothesis (load-bearing, refereed).** The class is a **uniformly definable family**: $f_\theta(x)=F(\theta,x)$ with $F$ definable in an o-minimal expansion of the real field ($\Theta\subseteq\mathbb R^p$, $\mathcal X\subseteq\mathbb R^q$ definable; boundedness *not* required).

**Theorem DM-6.**
(a) The envelope maps $(x_1..x_k,\tilde y,x,\varepsilon)\mapsto\inf S_\varepsilon(\tilde y),\ \sup S_\varepsilon(\tilde y)$ are definable partial functions into $\mathbb R\cup\{\pm\infty\}$ (the section is a first-order projection; suprema of definable families are definable).
(b) For every finite smoothness order $r$: the domain admits a **finite** decomposition into definable cells on each of which the envelopes are $C^r$; the envelopes are continuous at every point of the finitely many open (full-dimensional) cells, and the complement of the open cells — where all discontinuities live — is a definable set of **lower dimension**.
(c) On any compact subset of any open cell, center and radius are uniformly approximable to arbitrary accuracy by finite-parameter continuous families (polynomials suffice — Stone–Weierstrass; this is an existence statement, not a method). Globally uniform continuous approximation remains impossible wherever a genuine jump exists (MP-4).
(d) The $\pm\infty$ loci (partiality) are themselves definable: the operator's domain flags are finitely describable.

**Two distinct tameness routes (refereed correction).** DM-6 does *not* cover the corpus's flagship nonparametric regimes: Lipschitz classes and RKHS balls are not definable families. Their envelopes are piecewise-explicit by **closed form** (McShane envelopes; Golomb–Weinberger interval) — a second, independent route to piecewise regularity with *different uniformity properties in $k$* (DM-6's cell count is non-uniform in $k$; the closed forms are explicit for every $k$). A learnability claim must state which route it invokes.

**Consequences.** Finite-parameter approximability of $\mathbb A|_{\mathcal C}$ holds — cell-wise, away from the lower-dimensional transition set — exactly when the class is tame by one of the two routes. Untamed smooth classes (C5-type flat families) offer no such structure: their transition sets can have interior, and no genericity or cell-wise approximation statement survives. **The learnability dividing line coincides with the identifiability dividing line: definability/quasianalyticity.**

---

## 2. Approximation error

The budget is a **nonlinear composition**, not a sum (refereed correction):

$$\text{prediction error}\ \le\ \underbrace{\tfrac12\,\omega_{x,D}\!\big(2\varepsilon+2h\big)+h}_{\text{representation error }h\text{ through the modulus}}\ +\ \underbrace{\eta_s}_{\text{selection slack, constant }1\text{ (tight)}}\,;\qquad \text{plus}\ \underbrace{\ge\tfrac12\,\text{jump}}_{\text{smoothing error at transitions}}\ \text{for continuous models.}$$

- The window term sits *inside* the modulus and can be infinite while $h$ is arbitrarily small (C3): certificate stability under meta-level error is governed by $\omega$, and no model class escapes this.
- The selection constant $1$ is sharp (MP-5).
- At **coverage/validity boundaries** the radius jumps to $+\infty$: the half-jump lower bound for continuous surrogates becomes vacuous/infinite there unless the codomain is compactified to $[0,+\infty]$ — the honest formulation makes the radius head a $[0,+\infty]$-valued map and treats the boundary as a definable flag set (DM-6(d)), not a regression target.
- One-sided semantics throughout: only outer envelope slack preserves certificates; symmetric error notions produce false certificates.

---

## 3. Sample complexity: the excess-information theorem (DM-7; refereed: confirmed as scoped)

**Setting: the exactly-$d$ linear class under F17's common-core pattern** ($n$ tasks, every design contains the core $X_0$, $|X_0|=d$; archive matrix rank $d$; covered set of size $N$).

(a) **Tasks at $\le d$ points teach nothing.** A task observed at $k_j\le d$ points (any fixed design; "generic" unneeded — refereed) can exclude only a null set of candidate subspaces; it cannot reduce the consistent set to a point. *Scope:* within the per-point-rank-$d$ subclass — against rank-degenerate candidates such a task can be refuting.
(b) **Necessity.** Identification of the full window system requires every non-core point to be observed by $\ge d$ tasks (F17's per-point rank iff needs $\ge d$ rows). Counting task–point incidences under the core pattern:
$$\sum_j\big(k_j-d\big)_+\;\ge\;d\,(N-d).$$
(The clip $(\cdot)_+$ is required — unclipped sums would be falsely relaxed by small tasks; refereed.)
(c) **Dimensional tightness.** $d(N-d)=\dim\operatorname{Gr}(d,N)$: the required *excess* observations — beyond each task's own $d$ degrees of freedom — equal exactly the dimension of the meta-object. **Only observations beyond a task's own degrees of freedom carry meta-information.**
(d) **Sufficiency at the same count.** Core-pattern configurations achieving equality exist (e.g. $d$ tasks each observing all of $U$), conditional on F17's rank conditions.
(e) **Scope and the general case (refereed corrections).** Necessity is confined to the common-core pattern (beyond-core identification is the corpus's OPEN item, handoff §9.4). For general meta-parametric classes the naive bound $\sum(k_j-d)_+\ge M-g$ is **false** (counterexample: classes whose per-task parameters do not affect outputs — one task with $k_1=M$ observations identifies an $M$-dimensional meta-parameter). The provable general bound is only $\sum_jk_j\ge M-g$ for continuous recovery of the meta-parameter modulo a gauge of dimension $g$ (invariance of domain, under the hypothesis that the gauge quotient is locally Euclidean of dimension $M-g$). The per-task surcharge $d$ is class-dependent.

**Task-diversity requirements.** $n\ge d$ tasks with rank-$d$ archive matrix (necessary within the class, F17); a single-task archive with $d\ge2$ leaves realized ambiguity of dimension $\ge d-1$; archive rank certifies the class dimension from below only.

---

## 4. Support-size requirements (per task)

Established, restated at this phase's budget: $k\ge d$ for families containing a continuously and injectively parametrized $d$-cell (F14(i)); constrained families can identify below full rank (F6(iii)); generic-design *global* identification for analytic separating families needs $k\ge2d+1$ generically, attained sharp at $d=1$ (C4); at $k\le5$: exact stable recovery only for $d\le5$, generic-global guarantees only for $d\le2$, per-configuration verification required for nonlinear $d\in\{3,4,5\}$.

---

## 5. Identifiability of latent states

1. Latents are pinned at most modulo decoder-gauge, and only through the **closed, budget-truncated (size $\le k+1$; DM-3), covered-window quotient**. Training signals of any kind cannot do better — this is the ceiling of the learnable, not of the learner.
2. **Cross-model latent comparison is meaningless**; only induced operators are comparable, in the DM-9 metric ($\sup$ over realizable covered configurations of Hausdorff distance between reported envelope intervals — a metric on the double gauge quotient). Protocol NP-6 operationalizes this.
3. Finite-radius targets can live on a null set of queries (measure-zero validity regions, `adaptation_operator_parameterization.md` §7): a learnable regression signal for the *center* exists only where the radius is finite; elsewhere the only correct signal is the flag.

---

## 6. Summary — when is $\mathbb A$ finite-parameter approximable?

| Condition | Role |
|---|---|
| declared class, tame by definability (DM-6) or closed form | necessary for cell-wise uniform approximation with finitely describable flags |
| representation dimension $\ge$ class dimension (cont & def categories), with topological excess where obstructed | necessary (DM-2) |
| archive: core pattern, $n\ge d$, clipped excess $\ge d(N-d)$ | necessary & sufficient (within scope) for pinning the meta-object (DM-7) |
| per-task support $\ge d$ ($2d+1$ nonlinear-generic) | necessary for task identification (F14, F12) |
| outer (one-sided) approximation of envelopes | necessary for valid certificates |
| discontinuity accommodation (cell-wise targets; compactified radius) | necessary — continuous surrogates carry irreducible transition error (MP-4) |
| stability of the family ($\omega$ finite near $0$) | necessary for *any* nonzero-error tolerance to the meta-representation (C3 otherwise) |
