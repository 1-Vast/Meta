# Continuous Affinity Scope (Part 3)

> **Status:** Phase-19.5, 2026-08-03. Closes audit gap 3 (interval-censored calibration invalid). The audit's counterexample is adopted as the permanent witness; the scope is chosen to be the largest one the mathematics supports — no arbitrary censoring assumption is added. New results **DT-5–DT-6**, tagged **[proved] / [declared] / [rejected]**.

---

## The choice: **A, with a two-channel refinement** — point-supervised calibration; censored tasks feed estimation only

**Rejected options, with reasons.**
- **B (new elicitable loss for interval-censored supervision) [rejected].** A proper elicitation theorem under censoring requires a declared censoring mechanism (coarsening-at-random or stronger) — precisely the "arbitrary assumption" the mandate forbids; without one, the audit's counterexample is fatal and general: with every observable target the compatible interval $[0,1]$ and violation = distance-to-compatible-region, the degenerate zero-width band inside $[0,1]$ achieves zero loss — the objective elicits nothing about any latent law, and **no** loss can elicit latent quantiles from censored data without a mechanism assumption (the censored likelihood is not identified; this is the frozen partial-identification fact wearing a loss-function costume). **MI-12's calibration claim for censored supervision is retracted.**
- **C (distributional prediction, CRPS-type) [rejected].** Same obstruction: distributional scores are proper for *observed* outcomes; censored outcomes re-import the mechanism problem. Also heavier than the decisions require.

**Definition DT-5 (the adopted scope). [declared]**
Split the meta-training tasks by their own identification status at the queried functional:
- **Point-supervised channel** ($Q_T$ point-identified by the task's own data): these tasks — and only these — enter the **calibration-bearing loss**: the classical interval score at declared level $\alpha$, for which the central-quantile elicitation theorem is valid as stated (the audit's "pass for point supervision"). All calibration diagnostics (fiber-wise coverage vs $1-\alpha$) are computed on this channel.
- **Censored channel** ($Q_T$ interval/set-identified): these tasks contribute **only** to the estimation of population bands via their forced/compatible indicators — the DR-L/Manski route, whose outputs are *confidence-typed* interval statements (coverage of the identified population quantity), never calibration statements about a latent law. Their information is used; their supervision is never scored as if observed.

**Theorem DT-6 (the scope is honest and maximal). [proved]**
(i) *Validity:* on the point channel the elicitation theorem is classical; on the censored channel the interval bounds are the proved Phase-8.1/9 constructions; no statement crosses channels — the audit's counterexample is unconstructible because censored targets never reach the elicitable loss.
(ii) *Maximality without new assumptions:* any calibration claim from censored supervision requires a censoring-mechanism axiom (the counterexample shows loss-design alone cannot substitute), so no larger calibration scope exists within the program's assumption discipline; and discarding censored tasks entirely would waste identified information the estimation channel provably uses. The two-channel split is thus the unique scope that neither overclaims nor discards.
(iii) *Bounded scalar regression itself:* carried unchanged — declared value interval, $W_1$-closed CDF-band classes at declared mesh, Lipschitz losses, tagged (IID)/(C-IID) statistics — the audit's "conditionally supported" core, now with the calibration boundary drawn exactly. $\square$
