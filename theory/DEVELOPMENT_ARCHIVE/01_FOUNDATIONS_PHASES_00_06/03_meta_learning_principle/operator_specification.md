# Operator Specification — The Mathematical Interface

> **Status:** Phase-3 interface, 2026-08-02. This file answers one question for a future researcher: **"What mathematical object must a model approximate?"** — never "how should a model be built." Sources: F-numbers (`../02_theory_refinement/theorem_formalization.md`), OP/CP-numbers (`../02_theory_refinement/`), MP-numbers (`meta_learning_abstraction.md`). All requirements below are theorems or refereed propositions of the corpus; none is a design preference.

---

## 1. The approximation target

**Primary target: the closed window system**
$$\mathbb T=\big(T_S\big)_{S\ \text{finite}\subset U},\qquad T_S=\{f|_S:f\in\mathcal F\}\subseteq\mathbb R^{S},$$
on the covered region $U$ of the archive. Everything else is *derived*: by CP-1/CP-2 and OP-9 the entire inference problem — sections, risks, optimal rules — is a function of the windows and of nothing finer (minimal sufficiency, with the $\varepsilon=0$ clause). A model that approximates $\mathbb T$ has approximated everything the theory permits knowing at the meta level.

**Equivalent surface form: the envelope pair.** For each support configuration $D$, query $x$, noise $\varepsilon$, and support values $\tilde y$, the windows induce the section $S_\varepsilon(\tilde y)$ and thence the pair
$$\Big(\underline E,\ \overline E\Big)(D,\tilde y,x,\varepsilon)\;=\;\Big(\inf S_\varepsilon(\tilde y),\ \sup S_\varepsilon(\tilde y)\Big),$$
from which prediction and certificate follow: $\mathbb A=\tfrac12(\underline E+\overline E)$, radius $=\tfrac12(\overline E-\underline E)$. **A model must approximate the pair, not the center alone:** the point-prediction operator loses the certificate, and distinct window systems can share centers while differing in radii. Approximating $(\underline E,\overline E)$ on encountered configurations is equivalent to approximating the sections, which is the usable face of $\mathbb T$.

**What is *not* a target.** The member set of the family (unidentifiable even in principle — CP §2.3); the parametrization/gauge (CP §4.5); windows off the covered region (F18: any extension off $U$ is an assumption, to be declared, not learned).

---

## 2. Signature

$$\mathbb A:\quad \underbrace{\mathbb T}_{\text{meta}}\ \times\ \underbrace{\Big(\bigcup_{k'\le k}(\mathcal X\times\mathbb R)^{k'}/\mathfrak S_{k'}\Big)}_{\text{support set }S_t}\ \times\ \underbrace{\mathcal X}_{\text{query }x}\ \times\ \underbrace{[0,\infty)}_{\varepsilon}\ \dashrightarrow\ \underbrace{\mathbb R\times[0,\infty]}_{(\text{center},\ \text{radius})}$$

- The support enters as a **finite set** (permutation invariance is a theorem, OP-7(i), so the quotient is exact, not a modeling choice).
- $\varepsilon$ is an **explicit argument**: no single rule is optimal across noise levels (F1 Rem. 1.4); an interface without $\varepsilon$ cannot represent the optimal operator.
- The arrow is **partial** ($\dashrightarrow$), with three structural sources of undefinedness that the interface must *surface, never silently extrapolate*:
  1. **unrealizable data** — $S_\varepsilon(\tilde y)=\emptyset$ (support inconsistent with every member at level $\varepsilon$); doubles as a misspecification detector; the recorded convention (project to the realizable set) costs a factor of two in the modulus argument (treatise §10.4);
  2. **off-coverage** — $(D\cup\{x\})\not\subset U$: the windows are unknown there and unknowable (F18);
  3. **unbounded sections** — covered but unidentifiable queries (e.g. $\phi(x)\notin\operatorname{row}(G)$ in the unconstrained linear regime): the center is undefined and radius $=+\infty$ is the **legal, required output**.

---

## 3. Axioms (necessary conditions any admissible approximation must satisfy)

Each axiom is a certified property of the true operator; violating it certifies error.

- **A1 (Permutation invariance).** Invariance under reordering of support pairs. *(OP-7(i).)*
- **A2 (Affine equivariance).** For $\alpha\ne0,\beta$: transforming all values ($T\mapsto\alpha T+\beta$, $\tilde y\mapsto\alpha\tilde y+\beta\mathbf1$, $\varepsilon\mapsto|\alpha|\varepsilon$) maps center $\mapsto\alpha\cdot\text{center}+\beta$, radius $\mapsto|\alpha|\cdot$radius — including reflections $\alpha<0$. *(OP-7(ii).)*
- **A3 (Gauge invariance).** Dependence on the family only through windows; reparametrization acts trivially. *(CP-1.)*
- **A4 ($\varepsilon$-monotonicity).** Sections nested nondecreasing in $\varepsilon$; radius nondecreasing in $\varepsilon$. *(Definition of $S_\varepsilon$.)*
- **A5 (Support monotonicity of the guarantee).** Adding support points shrinks sections: radius nonincreasing under support refinement. **The center is exempt** — estimates are certifiedly non-monotone (OP-8; nested $[0,10],[9,10],[9,9.2]$ → centers $5,9.5,9.1$). An interface must not impose monotone updates.
- **A6 (Projective consistency of the meta-object).** $\operatorname{proj}_S\widehat T_{S'}=\widehat T_S$ for $S\subseteq S'$ on represented windows. Consistency is necessary; it is **not sufficient** for realizability by a family (MP-1(iii), Waterhouse) — countable coverage or compact windows restore sufficiency (MP-1(i),(ii)); at $\varepsilon>0$ windows are determined only up to closure, so the represented system should be closed.
- **A7 (Reproduction).** On exact realizable data at identifiable configurations: $\mathbb A(\mathbb T;f|_D,x,0)=(f(x),0)$. *(F2 — the cheapest necessary test of any candidate.)*
- **A8 (Lipschitz selection).** The center is $1$-Lipschitz in the section w.r.t. Hausdorff distance. *(OP-5; constant $1$ sharp.)*
- **A9 (Sensitivity sum — derived, conditional).** When constants lie in the family, translation equivariance ($\alpha=1$ case of A2) forces any locally Lipschitz realization of the operator to have a.e. support-value sensitivities summing to $1$; in the linear regime the sensitivity profile at $x$ is exactly $w(x)=(G^+)^\top\phi(x)$ on the validity region. *(F8, OP-7, MP-6.)*
- **A10 (Partiality surfacing).** The three undefined regions of §2 must be reported as such. A center emitted where the true radius is $+\infty$, or off coverage, is fabrication, not approximation.

---

## 4. Certificate semantics: one-sided error

The radius is a certificate, and certificates have **one-sided semantics**:

**Outer-enclosure requirement.** An admissible approximation must satisfy $\widehat{\underline E}\le\inf S_\varepsilon(\tilde y)$ and $\widehat{\overline E}\ge\sup S_\varepsilon(\tilde y)$ wherever it reports at all. Then the true value always lies in the reported interval: the certificate is *valid* and merely conservative. Any under-approximation is a **false certificate**: there exists a consistent member whose value escapes the reported interval (immediate from the definition of the section). Symmetric (Hausdorff-ball) error notions destroy certificate semantics; the interface therefore demands conservative envelopes, with tightness a quality metric, validity a requirement.

**Inner/outer tension at the meta level.** A finite archive of $n$ tasks yields an **inner** approximation of each window (observed traces are a *subset* of $T_S$). Inner approximations under-cover sections and thus produce false certificates. The conversion from inner to outer is exactly what a **model-class closure assumption** buys (e.g. F17's exactly-$d$ linear class: $n\ge d$ spanning tasks pin the window *equal to* a subspace). The interface makes this explicit: *no closure assumption, no valid radius* — a model must either declare its class or output only uncertified centers.

---

## 5. Error calculus

Two factors, with sharply different regularity — the composition, not either factor alone, is the honest calculus:

1. **Selection factor (benign, constant 1 — tight).** If the surrogate section satisfies $d_H(\widetilde S,S)\le\eta$ (both nonempty bounded), then $|\operatorname{cen}\widetilde S-\operatorname{cen}S|\le\eta$ and the surrogate's true conditional worst-case error is $\le\tfrac12\operatorname{diam}S+\eta$; tight ($S=[0,1]$, $\widetilde S=[\eta,1+\eta]$). *(The error-transfer lemma — labeled **MP-5** and defined here; selection step OP-5; refereed with tightness.)*
2. **Window factor (governed by the modulus — can be infinite).** A Hausdorff-$h$ error in the joint window inflates sections into $S_{\varepsilon+h}\pm h$, giving prediction error
$$\le\ \tfrac12\,\omega_{x,D}\big(2\varepsilon+2h\big)\;+\;h,$$
which for identifiable-but-unstable families (C3, tanh) is **infinite for every $h>0$**: certificate stability under meta-level error is itself controlled by the trace modulus, and no interface can promise better.
3. **Archive-noise constant: OPEN.** With inexact archives, only the factorization and a first-order expectation are established; the quantitative constant is explicitly open (`operator_formulation.md` §6). The interface must carry this flag rather than assume an additive bound.

**Obstruction to continuous surrogates (MP-4).** The true operator is discontinuous in the support values at section-topology transitions; any continuous approximator has irreducible sup-error $\ge$ half the local jump. Approximation targets should therefore be the *envelopes* (whose jumps are one-sided and locatable at realizability boundaries) with declared discontinuity locations, or accuracy claims must exclude transition neighborhoods.

---

## 6. Admissibility and correctness

**Definition (admissible approximation).** A tuple $(\widehat{\mathbb T},\widehat{\underline E},\widehat{\overline E})$ together with domain flags, satisfying A1–A10 and the outer-enclosure requirement on its declared domain, with a declared model class converting inner archive evidence to outer windows.

**Definition (correct at level $(\eta_{\mathrm w},\eta_{\mathrm s})$).** Window error $\le\eta_{\mathrm w}$ (Hausdorff, per represented window) and envelope slack $\le\eta_{\mathrm s}$ (excess of the reported over the true envelopes). Then by §5 the prediction error is bounded by $\tfrac12\omega(2\varepsilon+2\eta_{\mathrm w})+\eta_{\mathrm w}+\eta_{\mathrm s}$ at identifiable configurations, and every reported certificate is valid.

**What "success" can never mean.** Beating the radius floor $\tfrac12\omega(2\varepsilon)$ in the worst case (F1); extracting more than $k$ continuous task dimensions from $k$ observations (F20/CP-3); emitting certified values off coverage (F18); or dispensing with $\varepsilon$ (F1 Rem. 1.4). A reported violation of any of these identifies an undeclared assumption, not a transcended theorem.

---

## 7. Interface summary

| Question | Answer |
|---|---|
| What must be approximated? | The closed, projectively consistent window system on the covered region — surfaced as conservative envelope pairs $(\underline E,\overline E)$ per configuration. |
| What is the prediction? | Center of the reported interval; the radius ships with it. |
| What is learnable from the archive? | Windows on $U$, exactly under F17's conditions; nothing off $U$ (F18). |
| What does the support buy? | A section cut of codimension $\le k$; at most $k$ continuous task dimensions. |
| What must be exposed, not hidden? | $\varepsilon$; the three partiality sources; the model-class closure assumption; discontinuity locations. |
| What invariances are theorems? | A1–A5, A7–A9 (and A6 as a consistency law of the meta-object). |
| What error calculus applies? | Constant-1 selection factor; modulus-governed window factor (possibly infinite); archive-noise constant open. |
| What can falsify a candidate? | Any violation of A1–A10, any false certificate, or any of predictions P1–P10 (`meta_learning_abstraction.md` §6.9). |
