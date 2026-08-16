# Meta-Learning Compatibility (§4)

> **Status:** Phase-17, 2026-08-03. Learning defined as risk minimization over tasks — mathematical existence only, no architecture. New results **MR-11–MR-13**, tagged **[proved] / [conditional] / [declared]**.

---

## 1. The learning problem

- **Task distribution $p(T)$:** the observable task law $\Pi_{\mathrm{obs}}$ on $\mathbb T$ (marked lifts remain declared-model objects; learning touches observables only).
- **Support $S_T$:** the task's released record (its own few-shot observations) — the argument of the operator.
- **Query supervision $Q_T$:** the task's *identified* query information — the value / interval / order-set that the task's own data determine under the frozen theory (fully identified where the per-task gate holds; the forced/compatible interval or admissible-order set otherwise). $Q_T$ is observable; no latent mark is consumed.
- **Operator:** $A_\theta(S_T)\in\mathbb M$-values (a band point $b_\theta(S_T)\in\mathbb B$ per declared index), $\theta\in\Theta=[0,1]\times\mathbb B^m$.

**Definition MR-11 (operator-level loss). [declared class; properties proved]**
A declared loss $L:\mathbb B\times\mathcal Y\to[0,\infty)$, **convex and Lipschitz in the band argument**, scoring the emitted description against the task's identified query information — the canonical instance being the interval/band score family: per constrained event,
$$L\big(b,\,y\big)\ =\ \underbrace{(u-l)}_{\text{width}}\ +\ \tfrac{2}{\alpha}\,\underbrace{\mathrm{dist}\big(y,\,[l,u]\big)}_{\text{violation}},$$
summed over the declared events/grid, evaluated against $y=Q_T$ (censored tasks score against their intervals: violation measured to the *compatible* region, width unchanged — every quantity observable). Width is affine and the violation convex in $(l,u)$: $L$ is convex, Lipschitz on the compact $\mathbb B$. **[proved]**

**Definition (the learning rule).**
$$\theta^\star\ \in\ \operatorname*{arg\,min}_{\theta\in\Theta}\ \ R(\theta)\ :=\ \mathbb E_{T\sim p(T)}\ \big[\,L\big(A_\theta(S_T),\,Q_T\big)\,\big].$$

## 2. Existence — and more: convexity

**Theorem MR-12 (existence, convexity, and attainment). [proved]**
(i) $\theta\mapsto A_\theta(S_T)$ is **affine** (MR-3's decoder is a convex combination with $\theta$-affine coefficients in the band slots and bilinear only through the scalar $\lambda$; at fixed $\lambda$ it is affine, and jointly it is affine in $(\,(1{-}\lambda)\text{-rescaled slots}\,)$ after the standard reparameterization $\tilde b_j=\lambda b_j$ — the feasible set of $(\lambda,\tilde b)$ remains compact convex since $\mathbb B$ is a polytope containing the rescaling cone sections). Composed with convex $L$: $T$-wise loss convex in the parameters; expectation preserves convexity; $R$ is convex and Lipschitz on a compact convex domain.
(ii) Hence a global minimizer $\theta^\star$ **exists** (Weierstrass), the minimizing set is convex, and the optimization layer is a **convex program** — the Phase-15 "efficiency" residue shrinks to the classical solvability of convex programs, a mathematical fact rather than a declared hope. Measurability of $(\theta,T)\mapsto L(\cdot)$ (needed to type $R$) holds by MR-9(iii) + continuity of $L$. $\square$

## 3. What risk minimization learns — the honest characterization

**Theorem MR-13 (elicitation, and the certificate/preference separation). [proved / declared as marked]**
(i) *[proved — classical elicitability]* The population minimizer of the interval score at level $\alpha$ elicits **central quantile bands** of the conditional law of $Q_T$ given the operator's inputs: risk-minimization learning produces a *calibrated population band*, not the outer forced/compatible identified band — the two coincide only in degenerate cases.
(ii) *[declared resolution — the Phase-7 separation, applied]* This is not a defect but the program's own division of labor: the **learned band is decision information** (a preference/completion generator, rung-tagged, assumption-priced), while **certificates remain canonical** — the $\lambda=0$ component and the frozen identification layer carry the worst-case floors, outer envelopes, and honesty axioms, which no value of $\theta$ can alter (H1; the ledger's identification row is $\theta$-invariant by construction). A trained operator can tighten decisions; it can never tighten certificates. The loss therefore *may not* be used to emit worst-case claims, and the interface (`deep_learning_interface.md`) types this prohibition.
(iii) *[proved]* Ceiling compliance: for every $\theta$, decisions computed from $A_\theta$ are supported inside identified sets after the mandatory support restriction, and no unconditional guarantee below the frozen radius can be emitted — the Phase-7/8 theorems apply verbatim because $A_\theta(H)\in\mathbb M$ always (MR-9). Learning cannot break honesty even at adversarial $\theta$: the worst a bad parameter can do is give bad *advice*, priced by $R(\theta)$, never a false *certificate*. $\square$
