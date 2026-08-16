# Final Handoff Object (Part 5)

> **Status:** Phase-19.5, 2026-08-03. The directly implementable object, assembled from the four closures. Every slot cites its theorem; every assumption is named; nothing is left to inference. **DT-9**.

---

## Input

$$\big(\,z_H,\ S,\ Q,\ \gamma\,\big)$$
- $z_H$: frozen deployment state — trained parameters $\hat\theta$, fiber counts $(N_c)$, per-context population bands with margins, assumption tags (DT-1; vacuous state = declared fallback).
- $S$: current support — finite **set** of $(x,\tilde y)$ pairs, $k\le5$, plus mandatory $\varepsilon$ and optional declared auxiliary label; consumed through exactly two computations: the identification object $I(S)$ (exact intervals / order sets — an explicit typed component) and the context $\kappa(S)$. Invariances: permutation, duplicate idempotence, gauge, channel typing (MI-3, DT-0.4/0.5).
- $(Q,\gamma)$: declared query and specification indices (finite deployment skeleton; declared value grid with mesh $h$ for the scalar head; declared order events for the ranking head).

## Output

Per index, the typed triple plus certificates:
$$\Big(\ K(b_\theta)\big|_{\mathrm{supp}\,I(S)}\ \text{(valid band class: scalar CDF-bands / order polytope)},\ \ \mathrm{confidence}(z_H),\ \ \mathrm{rung}(z_H)\ \Big)\ \ +\ \ \text{certificate row from }I(S)\ (\theta\text{-invariant})\ +\ \text{flags}.$$
Validity, coherence, nonemptiness, $W_1$-closedness hold for **every** parameter and every implementation state (DT-3). Two heads, no cross-arrow: scalar affinity (Route B) and ranking (Route A) — DT-8.

## Metric

Operator-value metric per head: Hausdorff-$W_1$ on scalar band classes (stability $\varepsilon D_V+2h$, proved), Hausdorff-TV/Hoffman on order polytopes (constant $\tfrac12\bar H$, proved); plus confidence ($|\cdot|$ on $[0,1]$) and rung (discrete) coordinates with declared weights. Risk differences are controlled by these metrics for the declared Lipschitz (scalar) and decomposable/LP (order) loss classes.

## Loss

- **Scalar head, point-supervised channel only:** interval score at declared level $\alpha$ — convex, Lipschitz, elicits central quantile bands; calibration diagnostics on this channel alone (DT-5/6).
- **Scalar head, censored channel:** forced/compatible estimation updates to the population bands (confidence-typed; never scored as observed outcomes).
- **Ranking head:** order-event scores against per-task **identified order sets** (Kendall-type decomposable losses: pairwise marginal-sufficient; exact-match/top-$k$: order-polytope LPs) — separate supervision, DT-8.
- Optimization in the perspective variables (convex; minimizer set convex; DT-0.3).

## Assumptions (complete list — nothing else is assumed anywhere)

1. **Declared skeleton** (contexts, atlases, value grid/mesh, horizon, partition) — fixed per deployment; refinement = re-deployment (scope, per DT-0.1's corrected status).
2. **Bounded affinity interval** $V$; **Lipschitz loss class** in the value (scalar head).
3. **Task sampling:** (IID) or (C-IID-$\kappa$) with fiber counts, transport radius if declared; missing-fiber term always charged. Undeclared shift ⇒ proved adversarial failure mode + fallback.
4. **Conditioning rung declarations** where population conditionals are consumed ($\kappa$-DESIGN, SUFF-$\kappa$, LIK — each optional, each priced, each echoed).
5. **Implementation obligations, exactly two:** **(C3)** uniform coefficient accuracy for the implementer's *specified* class on the compact statistic domain (the transfer to operator/risk error is proved with explicit constants); optimization tolerance to the (reparameterized-convex) optimum.

## What the theory guarantees vs what it does not

Guaranteed unconditionally by type: validity of every emission; honesty of every certificate; correct fallbacks; no false claim reachable by any parameter or implementation state. Guaranteed conditionally (tags above): ERM generalization; population coverage; calibration on the point channel. Not guaranteed, by proof of impossibility or by named obligation: latent calibration from censored data (DT-6); ranking derived from scalar affinities (DT-7); (C3) for unspecified classes.
