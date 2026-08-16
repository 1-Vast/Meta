# Existence, Identification, Learning — Separated and Retyped (§4)

> **Status:** Phase-10, 2026-08-03. Repairs the audit's tier findings: existence was mistyped (latent target from an observable law — retracted), identification lacked a sharpness statement (re-scoped to outer semantics), learning was pointwise-only (upgraded in §5). New results **LC-11–LC-14**, tagged **[proved] / [conditional] / [declared] / [retracted]**. No tier's result is used to claim another's.

---

## A. Existence

**Retraction.** MC-16(i)/ML-L1 asserted the ideal target from a law on the observable $\mathbb T$. **Retracted**: $\mathbb T$ does not contain the mark, and observable laws admit distinct latent lifts (LC-3(iii)).

**Theorem LC-11 (existence, correctly typed). [proved]**
(i) *Ideal target — on the marked space.* For every **declared** marked law $\Pi^\bullet$ on the standard Borel $\mathbb T^\bullet$ (LC-1/LC-2), the ideal operator
$$M^\star_{\Pi^\bullet}\in\mathbb M,\qquad M^\star_{\Pi^\bullet}(c,Q,\gamma)=\Pi^\bullet\big(g_Q\text{-mark}\in\cdot\mid\kappa(O)=c\big)\ \text{(singleton class, top rung)}$$
exists by regular conditional probability, and satisfies the projective coherence of LC-5 (marginalization commutes with conditioning). Existence is relative to $\Pi^\bullet$ — never to $\Pi_{\mathrm{obs}}$.
(ii) *Estimator.* $A_\phi$ exists unconditionally as a total measurable map (LC-10) — a statement about a construction, needing no law at all. The two existence claims are different types and are kept apart. $\square$

## B. Identification

**Theorem LC-12 (identification — the lift class, with outer semantics declared as the target). [proved / declared]**
Given $\Pi_{\mathrm{obs}}$ (the $N\to\infty$ observable limit): the latent object is determined only up to the **lift class** $\mathcal L(\Pi_{\mathrm{obs}})=\{\Pi^\bullet:\mathrm{proj}_*\Pi^\bullet=\Pi_{\mathrm{obs}},\ \text{support coupling}\}$. Define the **identified operator**
$$M^\dagger_{\Pi_{\mathrm{obs}}}(c,Q,\gamma)\ =\ \text{the constraint class with eventwise endpoints}\ \big[\Pi_{\mathrm{obs}}(\text{record forces }E\mid c),\ \Pi_{\mathrm{obs}}(\text{record compatible with }E\mid c)\big],\ E\in\mathcal E_{Q,\gamma}.$$
Then: (i) **[proved]** for every lift $\Pi^\bullet\in\mathcal L$ and every index, $M^\star_{\Pi^\bullet}(\iota)\in M^\dagger(\iota)$ — the identified operator is a valid **outer** description of the whole lift class (pointwise: $l_i\le\mathbf 1\{\text{mark}\in E\}\le u_i$, integrate); (ii) **[proved]** eventwise the endpoints are attained by admissible lifts (put conditional mass at the extreme admissible marks — the LC-3(iii) construction), so each single evaluation-interval is sharp; (iii) **[declared, honest]** *joint* sharpness — whether every point of the polytope at every index is realized *simultaneously* by one lift respecting all projective-coherence constraints — is **not claimed**; the compilation target is $M^\dagger$ itself under outer semantics, which is all the decision layer ever consumes (one-sided validity, as everywhere in this program). The audit's sharpness gap is thereby closed by re-scoping the claim, not by an unproved theorem. $\square$

**Corollary LC-13 (point identification). [proved]** $M^\dagger(\iota)$ is a singleton iff the per-task censoring width for $\iota$'s events vanishes $\Pi_{\mathrm{obs}}$-a.s. — the frozen per-task partial identification, integrated; unchanged from Phase 9, now stated at the correct (marked/observable) types.

## C. Learning

**LC-14 (what the learning tier must now provide — contract for §5). [declared]**
The learning statement must control $A_\phi(H_N)$ **as an element of $(\mathbb M,d_{\mathbb M})$** against the target $M^\dagger_{\Pi_{\mathrm{obs}}}$: a uniform-over-$\mathcal I$ deviation bound and a consistency statement, under declared task-level sampling, complexity, and concentration assumptions — not a fixed finite family of event evaluations. Delivered as Theorem LC-15/16 (`operator_learnability_theorem.md`). Nothing from tiers A or B is used there except as the *definition of the target*; in particular, existence of $M^\star$ (tier A) is never invoked to assert convergence (tier C), and convergence is toward $M^\dagger$ (tier B's object), never toward any single lift's $M^\star$ — the identification width is not, and cannot be, closed by sampling. $\square$
