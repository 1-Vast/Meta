# Permutation Invariance (§3)

> **Status:** Phase-18, 2026-08-03. The required symmetry, proved; and the mandated classification (DeepSets-like / attention-like / kernel embedding / other) resolved as a **mathematical equivalence, not a necessity**. New results **MI-7–MI-9**, tagged **[proved]**.

---

## 1. The invariance theorem

**Theorem MI-7 (permutation invariance is exact and required). [proved]**
For every permutation $\pi$ of the observation list and every $\theta$:
$$A_\theta\big(\pi(S)\big)\ =\ A_\theta(S).$$
*Proof.* By MI-6, $A_\theta=D_\theta\circ r$; both components of $r$ are symmetric functions of the observation list — $b_{\mathrm{can}}$ because the identified set is an intersection of per-observation constraints (order-free; MI-2(i)) and the forced/compatible endpoints are functionals of that set; $z$ because counts, frequencies, and $\kappa$ (a declared function of the *record as a set*) are symmetric. Hence $r(\pi S)=r(S)$ and the claim follows. Moreover the invariance is *required*, not optional: an implementation violating it would emit two different bands for one identical epistemic state — two different claims from the same information — which contradicts well-definedness of the operator (R0 of DE-P) before any question of accuracy arises. Duplicate idempotence (MI-3.2) is proved the same way: $\cap$ and the symmetric statistics are idempotent under exact repetition within a task. $\square$

## 2. The classification question, answered as equivalence

**Theorem MI-8 (the canonical representation is already a pooling form). [proved]**
$r$ decomposes elementwise-then-pool:
$$r(S)\ =\ \rho\Big(\ \bigoplus_{s\in S}\ \phi(s)\ \Big),$$
where $\phi$ maps a single observation to per-element features (its constraint interval data; its event indicators; its context features) and $\bigoplus$ is a commutative idempotent-where-required aggregation: **min/max pooling** for the identification component (canonical endpoints are order statistics — e.g. per-location interval $[\max_i\tilde y_i-\varepsilon,\ \min_i\tilde y_i+\varepsilon]$, and forced/compatible indicators are monotone lattice operations on per-element sets) and **sum pooling** for the statistical component (counts and frequencies are literal sums of indicator features). $\rho$ is the continuous read-out (margins, clipping, normalization). So the interface's representation *is* a sum/min/max-pooled symmetric functional — constructively, with named poolings, not by appeal to a representation theorem. $\square$

**Theorem MI-9 (sufficiency of sum-decomposable classes; no architecture is forced). [proved / classical, scoped]**
(i) *Sufficiency:* since the support cardinality is bounded ($k\le5$, frozen budget), the classical sum-decomposition results for permutation-invariant functions on bounded-cardinality sets apply without the known uncountable-domain pathologies (the latent-dimension caveats arise for unbounded set sizes): every continuous invariant target on the bounded-size support domain admits a representation $\rho(\sum_i\phi(s_i))$ with finite latent dimension — and MI-8 exhibits the specific target in that form directly, so no existence theorem is even needed for *this* interface.
(ii) *Equivalence, not necessity:* "DeepSets-like" (sum-pool), "attention-like" (permutation-**equivariant** interaction followed by an invariant pool — note the min/max poolings of MI-8 are interactions expressible this way), and "kernel mean embedding" (sum-pooling in a feature space) are three carriers of the same mathematical class: symmetric continuous functionals of bounded finite sets. Each can express $r$; none is mathematically forced. The *requirements* are exactly: (a) extensional permutation invariance and duplicate idempotence (MI-7), (b) sufficiency for $r$ (MI-6), (c) the codomain typing through the convex assembly (MR-14). Any structure satisfying (a)–(c) is admissible; the choice among them is an implementation matter that this interface deliberately does not constrain. $\square$
