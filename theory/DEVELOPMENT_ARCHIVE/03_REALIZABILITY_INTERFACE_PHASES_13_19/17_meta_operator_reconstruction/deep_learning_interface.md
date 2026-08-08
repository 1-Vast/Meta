# Deep Learning Interface (§5)

> **Status:** Phase-17, 2026-08-03. The exact mathematical interface a future trainable system must implement. Mathematics only: what must be approximated is stated as a function class and a factorization, never as an architecture. New result **MR-14**.

---

## The interface

**Input.** Support observations: the observable record of the current task (and, at meta-training time, of historical tasks) — designs, values, noise level, auxiliary labels. Consumed through the two frozen channels: identification (set-image) and population (multiset/statistic + context $\kappa$). Nothing latent, nothing future, no task identity.

**Output.** A structured object per declared (context, query, specification) index:
$$\big(\,b\in\mathbb B\ \text{(denoting the class }K(b)\text{)},\ \ \text{confidence},\ \ \text{rung}\,\big)\ +\ \text{the canonical certificate row (}\theta\text{-invariant)},$$
i.e. a point of the **valid-description polytope** with its side-channels — never a raw law, never an unconstrained vector, never a certificate produced by the learned component.

**Loss.** The declared convex operator-level band score of MR-11, evaluated against the task's identified query information (intervals/order-sets for censored tasks). Operator-level means: the loss consumes the emitted description, not any internal representation — internal coordinates remain gauge (MC-7) and may not be scored or interpreted.

**MR-14 (what a trainable model must approximate — the exact statement). [proved]**
By MR-3/MR-9, every valid operator in the family factors as
$$H\ \longmapsto\ \underbrace{\big(b_{\mathrm{can}}(H),\ z(H)\big)}_{\text{fixed, computable, contract-supplied}}\ \longmapsto\ \underbrace{\big(\lambda,\ w_1,\dots,w_m\big)}_{\text{the only approximable content}}\ \longmapsto\ \underbrace{(1-\lambda)\,b_{\mathrm{can}}+\lambda\textstyle\sum_j w_j\,b_j}_{\text{convex assembly in }\mathbb B\ \text{(fixed)}},$$
so the entirety of what a model must approximate is: **a continuous map from the compact statistic domain into a fixed finite-dimensional compact convex coefficient set** (mixing weight + partition weights + the learned anchor bands $b_j$ as trainable parameters of the assembly). Three consequences, each load-bearing for a future implementation and each already proved:
1. **Validity is architecture-independent.** Whatever function is put in the approximable slot, the assembly's output lies in $\mathbb B$ (MR-9): approximation error degrades *performance* ($R(\theta)$), never *validity* — no reachable state of any implementation can emit an invalid or incoherent object, because invalidity is unrepresentable in the codomain.
2. **Certificates are out of reach of training.** The canonical row and the identification channel are fixed computable components; the ledger's unconditional rows are $\theta$-invariant (MR-13(iii)). A model may not be trained, prompted, or tuned into a false certificate — the type system, not vigilance, enforces DE-L5.
3. **The learning problem is convex in the assembly parameters** (MR-12), and the only non-convexity an implementation can introduce lives in the approximable slot — where errors are performance-priced by construction.

**Failure semantics (carried).** Empty declared polytope $\mathbb B$ (deployment inconsistency), empty support-restricted class, unseen relevant fiber, untypable single-valued demand, tolerance infeasibility: the DC-A4/PM ledger flags, verbatim — the learned layer adds no new failure modes because it adds no new claim types.
