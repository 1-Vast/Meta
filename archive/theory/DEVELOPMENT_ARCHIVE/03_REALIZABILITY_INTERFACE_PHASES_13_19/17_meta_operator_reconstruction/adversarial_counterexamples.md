# Adversarial Counterexamples (§6)

> **Status:** Phase-17, 2026-08-03. Mandated self-refutation pass: each attack is constructed in earnest against the reconstruction; findings are repaired in place or absorbed as declared scope. Results **MR-15–MR-19**, tagged **[attack refuted] / [attack found — repaired] / [attack found — scoped]**.

---

## 1. Invalid parameters

**Attack.** Re-run the Phase-16 witnesses: endpoint order violation ($l=0.9,u=0.1$), coherence violation between a coarse event and its pullback, CDF non-monotonicity, empty decoded class.
**Outcome: [attack refuted].** All four are unconstructible: the offending coordinate combinations are excluded from $\Theta=[0,1]\times\mathbb B^m$ by the linear description of $\mathbb B$ (order, ties, monotonicity) and by the carried feasibility witness (nonemptiness); convex assembly cannot exit $\mathbb B$ (MR-9). Residual edge: $\mathbb B=\emptyset$ if the *declared* skeleton constraints are jointly inconsistent — an LP-checkable deployment error, flagged at declaration (MR-7), not a runtime state. A second edge probed: the $\lambda$-reparameterization in MR-12(i) must keep $(\lambda,\tilde b_j)$ feasible — verified: $\{(\lambda,\tilde b):\lambda\in[0,1],\ \tilde b\in\lambda\mathbb B\}$ is the convex hull of $\{0\}\times\{0\}$-section and $\{1\}\times\mathbb B$, a compact convex set; the assembly formula on it is affine. No invalid parameter exists. $\square$

## 2. Non-identifiability

**Attack (i): parameter gauge.** Distinct $\theta$ with identical operators (degenerate $\varphi$-weights; redundant anchors; $\lambda=0$ making all $b_j$ irrelevant).
**Outcome: [attack found — scoped, harmless by design].** $\theta$ is *deliberately* not an identified object: only $A_\theta$'s values matter (the frozen gauge philosophy, restated at MC-7/PM-9); the minimizing set of the convex program is a face of $\Theta$, any element equally good. No theorem in the package quantifies over "the" $\theta^\star$.
**Attack (ii): elicitation mismatch.** The risk minimizer targets calibrated quantile bands, not the outer identified band — a learner could masquerade preference as identification.
**Outcome: [attack found — repaired by typing].** MR-13(ii): learned bands are typed as decision information; certificates are canonical and $\theta$-invariant; the interface forbids scoring or emitting worst-case claims from the learned component. The mismatch is real, permanent, and *displayed*, not hidden. $\square$

## 3. Infinite-dimension leakage

**Attack (i): the sieve returns through the skeleton.** Demand accuracy $\to$ arbitrary precision; the skeleton must refine; $p$ grows.
**Outcome: [attack found — scoped, with the impossibility proved].** Exactly MR-4: no fixed $p$ covers all resolutions — proved, not hidden; the relaxation (fix the skeleton) is proved weakest (MR-5). Within a deployment there is no $\varepsilon$-indexed growth (MR-3(iii)); across deployments, refinement is a declared re-deployment. The verdict is downgraded accordingly (see `reconstruction_verdict.md`) — the honest classification, not a defect.
**Attack (ii): the $\top$-stratum exactness hole (Phase-16's finding).** Two histories, same statistic, different canonical margins.
**Outcome: [attack refuted — by removal].** The reconstruction never factors $b_{\mathrm{can}}$ through $z$: the canonical component is computed from $H$ directly (MR-3), and only the *learned* component reads the statistic. The false exact-factorization claim is retracted with its cause deleted.
**Attack (iii): ranking from continuous affinities.** Order laws derived from several continuous scalar values are not determined by scalar marginal bands (frozen DR-J2 lifts to laws).
**Outcome: [attack found — scoped].** The interface serves ranking through the Route-A order object with its own declared atlas (order events computed from per-task identified joint sets — the frozen $\Sigma$ machinery); deriving order bands from Route-B marginals is forbidden as a type error, and a joint continuous-vector object is **not claimed** — it would be a genuine theory extension, recorded as out of scope. $\square$

## 4. Hidden oracle assumptions

**Attack.** Audit every fixed component for uncomputable inputs: $b_{\mathrm{can}}$ (forced/compatible indicators, margins) — computable from observable records under the frozen constructions; $z(H)$ — counts and frequencies, observable; $\kappa$ — declared observable map; the loss — evaluated against per-task *identified* query information, never latent marks; $\mathbb B$ — explicit linear description from declared skeleton; $\varphi$ — declared fixed functions. The Phase-15 optimization oracle (global grid search against $A_\phi$ over all of $Z$) is **gone**: learning is expected-risk minimization over sampled tasks (MR-12), a convex program on observables; no component evaluates anything off the data or the declarations.
**Outcome: [attack refuted].** No oracle survives inspection. The only entities that are not data are declarations — and every declaration is echoed in the ledger by the standing H6 discipline. $\square$

## 5. Summary of the pass

| Attack | Result |
|---|---|
| invalid parameters (all four Phase-16 witnesses + two new edges) | refuted — unrepresentable |
| parameter gauge | found; harmless by design, scoped |
| elicitation masquerade | found; repaired by certificate/preference typing |
| accuracy sieve via skeleton | found; impossibility proved, relaxation declared — drives the verdict |
| $\top$-stratum exactness | refuted by construction removal |
| continuous-affinity joint ranking | found; scoped out with a typed prohibition |
| hidden oracles | refuted — none |

Two findings survive as *declared scope* (fixed skeleton; no joint continuous-vector ranking object), one as a *typed separation* (elicitation). None invalidates the construction; the first mandates the verdict's relaxation qualifier.
