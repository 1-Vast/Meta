# Approximation Conditions (§5)

> **Status:** Phase-18, 2026-08-03. The weakest theorem under which "a neural model approximates the mathematical operator" — stated as checkable conditions with explicit constants; **no universal-approximation property is claimed for any unspecified model class**. New result **MI-13**, tagged **[proved, conditional as marked]**.

---

## The theorem

**Theorem MI-13 (weakest sufficient conditions for interface approximation). [proved]**
Let an implementation supply a map $g_{\mathrm{impl}}$ in the *approximable slot* of MR-14 — the coefficient map from the compact statistic domain $Z$ into the fixed compact convex coefficient set $\mathcal C$ (mixing weight, partition weights, anchor selections) — and emit through the **fixed** contract components (representation $r$; convex assembly; canonical side channels). Suppose:

- **(C1) Extensional invariance.** $g_{\mathrm{impl}}$'s composite with $r$ is permutation-invariant and duplicate-idempotent in $S$ — automatic if the implementation reads $S$ only through $r$ (MI-6/MI-7); otherwise it must be verified extensionally.
- **(C2) Codomain typing.** Outputs pass through the convex assembly into $\mathbb B$ — then validity, coherence, nonemptiness, and closedness hold for *every* implementation state (MR-9): approximation error can degrade only risk, never validity.
- **(C3) Uniform coefficient accuracy.** $\sup_{z\in Z}\ \|g_{\mathrm{impl}}(z)-g^\star(z)\|\ \le\ \varepsilon$, where $g^\star$ is the target coefficient map (continuous on the compact $Z$; for the trained target, the $\theta^\star$-optimal map).

Then, with the stability constants of the reconstruction ($C_A=\tfrac12\bar H$ on Route A; $C_B=D_V$ plus the declared mesh floor $2h$ on Route B):
$$\sup_{S}\ d_{\mathbb M}\big(A_{\mathrm{impl}}(S),\ A_{\theta^\star}(S)\big)\ \le\ C_{A/B}\,\varepsilon\qquad\text{and}\qquad R(\mathrm{impl})\ \le\ R(\theta^\star)\ +\ \mathrm{Lip}(L)\cdot C_{A/B}\,\varepsilon.$$
*Proof.* (C2) makes both sides well-typed elements of $\mathbb M$ everywhere; band vectors are affine images of coefficients, so coefficient deviation $\varepsilon$ moves bands by $\le\varepsilon$ (convex-combination weights); the stability theorems transfer band motion to class motion with the stated constants; the loss is band-Lipschitz. $\square$

## What is deliberately *not* claimed

**(C3) is a property of the implementer's function class, to be established by the implementer** — by citing a classical density theorem for a *specified* class on the compact finite-dimensional $Z$ (e.g. polynomials or continuous piecewise-linear functions are dense in $C(Z)$ by Stone–Weierstrass — these are cited mathematical facts for *those* classes), or by direct verification. This interface proves: *if* (C1)–(C3) hold, the operator is approximated with explicit constants, *and* (C2) holds unconditionally by the contract's construction. It does **not** assert that "neural models" satisfy (C3): that would be a universal-approximation claim about an unspecified class, which the mandate forbids and honesty does not permit. The target $g^\star$ being continuous on a compact finite-dimensional domain is exactly what makes (C3) a *reasonable obligation* — the interface's job — while its discharge is the implementation's.

## The condition set is weakest

Dropping (C1) permits two outputs for one epistemic state — ill-defined before inaccurate (MI-7). Dropping (C2) re-opens every Phase-16 invalidity (the counterexamples return verbatim). Weakening (C3) to non-uniform accuracy breaks the sup-metric conclusion by the standard construction (concentrate the error near the adversary's statistic). Each condition is therefore individually necessary for the stated conclusion; jointly they are sufficient with explicit constants — the mandated "weakest theorem". $\square$
