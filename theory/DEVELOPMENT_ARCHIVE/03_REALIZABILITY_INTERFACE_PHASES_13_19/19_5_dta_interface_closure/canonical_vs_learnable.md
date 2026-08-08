# Canonical Operator vs Learnable Approximation (Part 2)

> **Status:** Phase-19.5, 2026-08-03. Closes audit gap 2 by separating, claim by claim, what is a **theory guarantee** of the canonical construction from what is a **model approximation assumption** about an implementation — and by typing the single unproved condition (C3) as the sole, named, checkable implementation obligation. New results **DT-3–DT-4**, tagged **[proved] / [conditional] / [obligation]**.

---

## 1. The two objects

- **Canonical operator $A^\star$:** the fully specified mathematical map of DT-1(iii) with the *canonical* coefficient assignment — population bands from $z_H$, the trained-family assembly at any $\theta$, canonical side channels, support restriction, certificate row. Every ingredient is closed-form and computable from $(z_H,S,Q,\gamma)$.
- **Learnable approximation $F_\omega$:** an implementation-supplied map occupying the **approximable slot only** (MR-14/MI-13): the coefficient function from the compact statistic domain into the fixed compact convex coefficient set. Its outputs pass through the fixed convex assembly and side channels; $\omega$ indexes an *unspecified* function class.

## 2. The claim ledger — theory vs assumption

**Theorem DT-3 (guarantees that hold for every $F_\omega$, unconditionally — by type). [proved]**
For any implementation whatsoever in the approximable slot: outputs lie in $\mathbb B$ (validity, coherence, nonemptiness, Route-B $W_1$-closedness — MR-9); the certificate row, floors, flags, support restriction, and fallbacks are computed canonically from $I(S)$ and $z_H$ and are $\omega$-invariant (no reachable implementation state emits a false certificate — the honesty firewall is structural); permutation invariance and duplicate idempotence hold whenever the implementation reads $S$ through the typed representation (and are extensional obligations otherwise); degradation at empty fibers/history is the declared fallback. **These are theory guarantees: they need no assumption about $\omega$.** $\square$

**Theorem DT-4 (claims that are conditional, with their exact conditions). [conditional / obligation]**
(i) *ERM generalization [conditional — proved under tags]:* within the mathematical family (perspective-reparameterized variables, DT-0.3), empirical risk minimization over $N$ tasks generalizes to population risk at dimension-$p$ rates under (IID)/(C-IID-$\kappa$) with the missing-fiber term (MI-11 as corrected). This is a theorem *about the family*, tagged by its sampling assumptions.
(ii) *Coefficient accuracy — the single obligation **(C3)** [obligation, deliberately unproved]:* $\sup_z\|F_\omega(z)-g^\star(z)\|\le\varepsilon$ for the implementer's *specified* class. The theory proves: the target $g^\star$ is continuous on a compact finite-dimensional domain; the transfer from coefficient error to operator-metric and risk error has explicit constants (MI-13); and density in $C(Z)$ — which discharges (C3) for any $\varepsilon$ — is a classical, *citable* property of specified classes (polynomials, continuous piecewise-linear). The theory does **not** assert (C3) for an unspecified "neural" class, and no honest mathematics could: **(C3) is the entire content of "the model works", isolated, named, and checkable — nothing else about the implementation is assumed anywhere in the stack.**
(iii) *Training-attains-ERM [obligation]:* that a practical procedure reaches the (convex, in the reparameterized variables) optimum to tolerance $\gamma^{\mathrm{opt}}$ — attainability is proved in principle (net argument), efficiency is the implementer's. $\square$

## 3. The separation, in one line

$$\text{Theory: validity, honesty, certificates, fallbacks — unconditional by type; \ ERM — conditional on declared sampling.}\qquad\text{Implementation: (C3) + optimization tolerance — two named obligations, nothing tacit.}$$
