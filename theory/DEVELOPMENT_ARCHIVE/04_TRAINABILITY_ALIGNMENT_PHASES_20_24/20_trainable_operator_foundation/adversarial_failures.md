# Adversarial Failures (Deliverable 7)

> **Status:** Phase-20, 2026-08-03. Mandated self-refutation over five named attacks. Each is constructed in earnest against the now-defined $\mathcal H$, $\Omega$, objective, and bridge; findings are repaired or absorbed as declared scope. Results **TF-16–TF-20**, tagged **[refuted] / [found — repaired] / [found — scoped]**.

---

## 1. Latent task vector collapse

**Attack.** Two tasks the theory should distinguish are mapped to the same output because the realization $G$ collapses their statistics; or a learned latent under-separates tasks.
**Outcome: [refuted].** The only task-dependent inputs are $z(S,Q,\gamma)$, $b^{\mathrm{pop}}_{\kappa(S)}$, and $I(S)$ — all fixed, computable, and $\omega$-independent. Tasks that these *separate* cannot be collapsed by any $G$: $G$ acts after $z$, and $I(S)$ (support restriction) and the certificate row bypass $G$ entirely, so two tasks with different identified objects produce different outputs regardless of $G$. Tasks that $z/I$ do *not* separate are provably interchangeable for the whole valid family (family-relative sufficiency, MI-6 as corrected DT-0.2 — sufficiency stands even where minimality fails). There is no additional learned latent vector in the interface (outputs are typed bands, not embeddings — MC-7); no slot exists in which collapse could occur. $\square$

## 2. Finite-dimensional bottleneck

**Attack.** The fixed $D$ hides a capacity failure: some target is unrepresentable, invalidating the interface.
**Outcome: [found — scoped, with theorems].** True and quantified: TF-11 (fixed-$D$ floor $\varepsilon_0(D)\ge0$, generically $>0$) and TF-12 (all-resolution impossibility under uniform regularity). This is a *capacity* statement, priced in advice quality ($R$), and it touches nothing else: validity holds at every $D$ and every $\omega$ (TF-5.4), certificates are $\omega$-invariant. The bottleneck is displayed as a floor theorem, not concealed; the resolution ($D(\varepsilon)\to\infty$ across skeletons) is stated. No invalidity, only an honest accuracy/capacity trade the engineer chooses. $\square$

## 3. Hidden task identity

**Attack.** $\omega$ (rich enough) memorizes training tasks — a lookup keyed by task identity leaks into deployment.
**Outcome: [found — repaired by type + priced by risk].** (i) *Structural:* the deployment input contains no task ID (MC-7(ii)); $z$ is a declared statistic, $I(S)$ an identified object — neither carries an identity key. A memorizing $\omega$ can only overfit the *coefficient map on $Z$*, which at deployment is evaluated on the current task's statistic, not a stored key. (ii) *Priced:* any memorization that does not generalize is penalized by the generalization gap $2\Gamma_N$ (TF-13) — the covering bound over the compact $\Omega$ makes train/population divergence bounded, so a lookup that helps training but not population risk is provably suboptimal for large $N$. (iii) *Firewalled:* even a maximally overfit $\omega$ cannot emit an invalid object or a false certificate (TF-5.4, $\omega$-invariant certificate row). Identity leakage is thus unrepresentable in the output type, structurally absent from the input, and risk-penalized in the objective. $\square$

## 4. Leakage

**Attack battery.** Query-answered-from-population (marginal-to-conditional); loss reading latent marks; $\kappa$ reading unidentified quantities; history multiplicities entering identification; the imitation design $\mu_Z$ smuggling test information.
**Outcome: [refuted, all].** The first four are re-verified against the new surface exactly as in MI-audit §3 (rung-typed population consumption; loss on identified targets only; $\kappa$ a declared observable map; set/multiset channel typing) — the realization $G$ adds no new read, since it consumes only $z$ (already channel-typed). The fifth is new to this phase: $\mu_Z$ is a **declared design over the statistic domain $Z$**, not over test tasks; TF-14(ii)'s $\sup$-upgrade uses only a mass floor on $Z$'s mesh, which is a resolution declaration, not access to query answers. No leakage channel survives. $\square$

## 5. Impossible approximation

**Attack.** The approximation theorem secretly requires an oracle (the target $g^\star$) unavailable at training, so TF-10/TF-14 are vacuous.
**Outcome: [found — repaired by construction].** $g^\star$'s *canonical* instance is computable from observables (the frozen forced/compatible construction + population state $z_H$) — so the imitation objective TF-14 is evaluable without any oracle; that is precisely why imitation is the route that closes the operator-metric arrow. For the *risk-optimal* $g^\star$ (task-supervised route), no oracle is used at all — only sampled tasks (TF-13). The one genuine limit is TF-11's floor (a target may sit outside a fixed family's closure), which is a capacity statement, not an oracle assumption. Approximation is therefore possible and oracle-free per tolerance, with the honest fixed-capacity caveat. $\square$

## Summary

| Attack | Result |
|---|---|
| latent task vector collapse | refuted — no collapse slot; separated tasks bypass $G$ |
| finite-dimensional bottleneck | found; scoped as the TF-11/12 floor theorems, validity untouched |
| hidden task identity | found; unrepresentable in type, absent from input, risk-penalized |
| leakage (five channels) | refuted — all channels typed, $G$ adds no read |
| impossible approximation | refuted — imitation objective is oracle-free; only the honest capacity floor remains |

No attack invalidates the definitions; the two findings that survive are the (proved, displayed) capacity floor and the (structurally firewalled, risk-priced) memorization possibility — neither leaves $F_\omega$ undefined or dishonest.
