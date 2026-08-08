# Research Stopping Criterion

> **Status:** Phase-5 determination, 2026-08-02. Question: is the current mathematical theory sufficient to begin engineering a learning system? The verdict below was audited against the corpus's complete open-item register (handoff §9) — the plan audit required every open item to be explicitly discharged or moved to the blocking list, and the YES to carry its three load-bearing qualifiers. It does.

---

## Verdict

**YES — for the scoped core program.** The theory is sufficient to begin engineering, for:

> **declared classes** (closure assumption stated), **tame** by the definable or closed-form route, **stability-certified** on the configurations where finite certificates are claimed, with **exact archives**, under the **worst-case, per-query (sup-loss), $k\le5$** scope — nonlinear classes admitted at $d\le2$ (generic) or with per-configuration identifiability verification for $d\in\{3,4,5\}$.

Three qualifiers are part of the verdict, not footnotes:
1. **Metric scope:** correctness is read on compact subsets of the *relative interiors of all cells* (else the center is untested where validity regions are null — stratum (i), $d>k$).
2. **Nonlinear qualifier:** $d\le2$ generic or per-configuration verified (the $2d{+}1$ budget arithmetic at $k\le5$).
3. **Effectivity qualifier:** the constructive guarantees are effective for semialgebraic (CAD-computable) classes; for general o-minimal classes, cell data exist but are not provided effectively by the corpus — an explicit assumption the engineering phase must either restrict to or supply.

**The exact mathematical objects that remain to be approximated** (nothing else remains at the theory level):

$$\textbf{O1}:\ \Phi\ (\text{class representation — closed, }(k{+}1)\text{-truncated covered-window quotient; outer-admissible decoder; DM-2 dimension budget}),$$
$$\textbf{O2}:\ U\ (\text{task adaptation — AP-1…AP-8}),\qquad \textbf{O3}:\ R\ (\text{query readout — center }+\text{ compactified one-sided radius}),$$
$$\textbf{O4}:\ V\ (\text{coverage/validity/realizability flags — one-sided, exhaustive}),$$

in the **branch normal form** of CR-1 (finitely many $C^r$ branches on definable cells jointly in class/support/query/$\varepsilon$, discrete selector, margin-guarded composition CR-4, outer rounding CR-6), under the constraint list `future_model_constraints.md` M1–M12, with convergence guaranteed by CR-6 and bounded by I-1…I-4.

**NO — for four named extensions**, each with its precise missing theorem:

| Extension | Missing theorem |
|---|---|
| Noisy archives | the quantitative subspace-perturbation constant (declared metric, norm, smallness regime, explicit $C(x,G)$ with remainder, or a matching lower bound) — the corpus's flagged OPEN item |
| Archives beyond the common-core pattern | the combinatorial characterization of identifying observation patterns (overlap chaining) |
| Distributional / average-case guarantees | any distributional foundation at all — deliberately out of scope; would require new axioms (task measure, noise law), not a patch |
| Unstable classes at unstable-local configurations | **nothing** — not open but *settled negatively*: CR-7(c) proves pointwise impossibility unconditionally (witness constructed); the engineering consequence is a prohibition (M-list), not a research gap |

---

## Discharge appendix: the complete open-item register (handoff §9)

| Open item | Disposition |
|---|---|
| 9.1 archive-noise constant | **Blocking for noisy archives** → NO-list. Non-blocking for the core program (exact-archive scope declared). |
| 9.2 exact class for anchored-projection optimality | Non-blocking refinement: the core program never relies on anchored optimality beyond its proven sufficient regime (MP-3). |
| 9.3 sharpness of $2d{+}1$ for general $d$ | Not dischargeable — absorbed into the verdict as the **nonlinear qualifier** ($d\le2$ generic / per-configuration verified). |
| 9.4 beyond-common-core archive identification | **Blocking for that extension** → NO-list. Core program uses the core pattern, where DM-7 is an iff. |
| 9.5 realization of window systems between countable and compact | Discharged for the core program: coverage is finite there (finite $U$), where realization is unconditional (MP-1). |
| 9.6 structure of the section-containment preorder | Non-blocking refinement: the radius quotient suffices for all core guarantees; the preorder refines support *selection*, not validity. |
| 9.7 joint-loss constants | Discharged by scope: the program is per-query/sup-loss, where exactness is proven (F1 Rem. 1.1); joint losses are declared out of scope, not silently used. |
| 9.8 $\varepsilon$-adaptive rules | Discharged by interface: $\varepsilon$ is an input (M4); per-$\varepsilon$ optimality is what the theory proves and what the interface demands; a single $\varepsilon$-uniform rule is a convenience question, not a correctness gap. |
| 9.9 F18 dichotomy beyond the linear class | Discharged by conservatism: O4 always flags off-coverage regardless of class (one-sided semantics); the general dichotomy would sharpen the *message*, not the *validity*, of the flag. |

No undischarged item touches the scoped core program. The claim "no missing theorem blocks the core program" is therefore exact, item by item.

---

## Closing statement of the mathematical program

The program set out to determine when a member of a function family can be inferred from $k\le5$ evaluations, and to derive — never posit — the structure a learning system would need. It ends with: an exact minimax theory whose single primitive is the trace modulus; a canonical operator forced by conditional minimaxity; a complete account of what archives can and cannot teach; dimension, sample, and capacity floors with matching constructions; a two-normal-form realizability theory with certified differentiable convergence exactly where its named impossibility theorems permit; and a constraint interface under which every future design choice is either mathematically free or certified wrong. The theory's own stopping criterion is met: what remains is approximation of four named objects under twelve named constraints — engineering — plus four named extensions, each waiting on one named theorem.
