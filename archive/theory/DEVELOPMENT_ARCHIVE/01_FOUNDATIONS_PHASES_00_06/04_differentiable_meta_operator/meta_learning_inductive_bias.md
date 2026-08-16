# Meta-Learning Inductive Bias — Derived Constraints

> **Status:** Phase-4 derivation, 2026-08-02. This file derives what any valid model of the canonical operator must preserve. **No modules are proposed; no networks are designed.** Every constraint is a theorem of the corpus or a refereed DM-result; violating a constraint certifies error, independently of how the model is built. Constraints are numbered IB-1 … IB-14; each carries its certificate.

---

## Structural constraints

**IB-1 (Support-set symmetry — composite level).** The composite operator is invariant under permutations of support pairs (OP-7(i)/A1). *Refereed scoping:* this constrains the composite, not each factor of a pipeline; an order-dependent intermediate composed with an order-discarding readout is admissible. Symmetric factorizations exist WLOG.

**IB-2 (Support-location dependence).** The task state is a function of the *labeled* sample — locations and values. Certificate: the one-observation witness ($V=\operatorname{span}\{t\}$: $y=1$ at $t=1$ vs $t=2$ gives different predictions at every query $x\ne0$). Value-only pooling is provably wrong.

**IB-3 (Query coupling).** The readout must vary with the query in three certified ways: the validity region $\{x:\phi(x)\in\operatorname{row}(G)\}$ is support-dependent; the sensitivity profile is the $x$-dependent weight vector $w(x)=(G^+)^\top\phi(x)$ on the validity region; the certificate (conditional radius) is query-dependent. No query-independent head can carry the certificate (MP-6).

**IB-4 (Gauge invariance).** Model outputs may depend on the family only through windows, and on the latent only through the decoder image. Training objectives can pin **only window-level invariants** — losses defined on latent coordinates optimize the unidentifiable (CP §4.5, lifted; `latent_operator_theory.md` §6).

**IB-5 (Affine equivariance, including reflections).** Value-space transformations $f\mapsto\alpha f+\beta$ ($\alpha\ne0$, negative allowed) must commute with the model as with the operator: center $\mapsto\alpha\cdot\text{center}+\beta$, radius $\mapsto|\alpha|\cdot$radius, $\varepsilon\mapsto|\alpha|\varepsilon$ (A2). *Derived tension (DM-2(c)):* equivariant representations may require strictly more latent dimensions than the Whitney minimum — equivariance is not free, and the tradeoff is a mathematical fact, not a design preference.

**IB-6 (Capacity ceilings as bias).** Two hard ceilings: (i) at most $k$ continuous dimensions of task identity from $k$ observations (F20/CP-3; protocol P2); (ii) the meta-object is pinned only up to the **size-$(k{+}1)$ window truncation** (DM-3; protocol NP-1) — structure beyond that is unlearnable, and a model representing it is representing its own prior. Enforcing the ceilings is valid bias; exceeding them in output is fabrication.

## Uncertainty constraints

**IB-7 (Certificate preservation).** Outputs are pairs (center, radius) with radius in the **compactified** $[0,+\infty]$; radius $=+\infty$ is a legal, required value (validity boundaries). One-sided semantics: reported envelopes must outer-enclose true sections — under-approximation is a false certificate. $\varepsilon$ is an input; radii nondecreasing in $\varepsilon$ (A4) and nonincreasing under support refinement (A5) — while the **center is exempt from monotonicity** (OP-8): enforcing monotone updates contradicts optimality (protocol P10).

**IB-8 (Partiality surfacing — unified A10).** Three undefined regions must surface as flags, never as silent extrapolation: (i) unrealizable support data — the misspecification detector; silently projecting instead costs the factor-two convention penalty (treatise §10.4) and destroys the detector; (ii) off-coverage configurations (F18); (iii) unbounded sections at covered-but-unidentifiable queries. The flag maps are finitely describable for tame classes (DM-6(d)).

**IB-9 (Discontinuity accommodation).** The operator is discontinuous at section-topology transitions (MP-4) and at coverage/validity boundaries (radius jump to $+\infty$). A valid model either carries piecewise structure aligned with the (lower-dimensional, definable) transition sets, or declares transition-band error; a globally continuous model owns an irreducible half-jump sup-error, localized and predictable (protocols P9 base-level; NP-3 meta-level — spikes as the *archive/window system* crosses a transition, distinguished from P9's support-value transitions).

**IB-10 (Selection regularity).** The center is $1$-Lipschitz in the section (Hausdorff distance; constant sharp — A8/OP-5): models whose reported center moves more than their reported envelopes is internally inconsistent — an auditable coherence condition between the two output heads.

**IB-11 (Sensitivity coherence).** When constants lie in the family: any locally Lipschitz realization must have a.e. support-value sensitivities summing to $1$ (A9, derived from translation equivariance); in the linear regime the whole sensitivity profile is pinned to $w(x)$ on the validity region. Testable by finite differences (protocol P3).

**IB-12 (Reproduction).** On exact realizable data at identifiable configurations: output $=(f(x),0)$ (A7/F2). The cheapest necessary test of any candidate.

## Shortcut taxonomy (derived; each with its certificate and test)

**IB-13 (Certified shortcut modes to exclude).**
1. **Baseline collapse.** Predicting $\operatorname{cen}T_x$ everywhere attains the $k=0$ risk and can look adequate on benign task collections while performing zero adaptation. Test: the MP-2 separation — measured gain vs the certified gain $\tfrac12(\operatorname{diam}T_x-\omega(2\varepsilon))$; a model whose gain is zero where the certified gain is positive has collapsed.
2. **Memorization without constraint.** Interpolating the support with no family constraint (unconstrained embedding of the task) yields sections of infinite radius off the support — certificates must expose it; a model reporting finite radii there is fabricating (F4).
3. **Coverage leakage.** Apparent accuracy off the covered region is an undeclared member-level assumption (F18; protocol P6): legitimate only if declared as part of the model class.
4. **Gauge overfitting.** Latent structure carrying no window-level invariant content; exposed by the DM-9 comparability metric (protocol NP-6): two trainings equal in induced-operator metric are *the same model*, whatever their latents.
5. **$\varepsilon$-ignoring (new; refereed).** A model trained at one noise level violating A4 under an $\varepsilon$-sweep; no single rule is optimal across noise levels (F1 Rem. 1.4) — an $\varepsilon$-free model is provably suboptimal at some level (protocol NP-4).
6. **Unrealizable-data suppression (new; refereed).** Silently projecting inconsistent supports instead of flagging the empty section — destroys the misspecification detector and incurs the doubling cost (test: engineered inconsistent supports expecting the flag; protocol NP-5).

**IB-14 (Scope of evaluation).** All constraints are worst-case; average-case benchmark performance can neither establish nor refute compliance (F3). Compliance is checked by the adversarial protocols P1–P10 and NP-1…NP-6 (`model_design_interface.md`), each with computable constants.

---

## One-paragraph summary

A valid model is one that *could be* the canonical operator: symmetric in the support as a composite, location-aware, query-coupled, gauge-blind, affinely equivariant, certificate-bearing with one-sided semantics and a compactified radius, flag-raising at its three undefined regions, piecewise-regular across certified discontinuities, capacity-capped at $k$ task dimensions and size-$(k{+}1)$ windows, sensitivity-coherent, reproducing on exact data — and falsifiable by the named protocols. Everything else about the model is unconstrained by the mathematics, and nothing else is.
