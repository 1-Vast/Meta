# Parameterized Meta-Operator Closure — Stopping Criterion

> **Status:** Phase-13 terminal decision, 2026-08-03. Sources: the three Phase-13 files (PM-1–PM-9) and the audit `../12_final_operator_audit/FINAL_VERDICT.md` (`OPERATOR_LEARNABILITY_INCOMPLETE`). Phases 0–11 unmodified; the operator not redesigned; no architectures. Retracted in this phase: OM-2's universal full-metric form (audit witness adopted); OM-5-pre's unconditional rung cancellation (restated conditional on the fiber event).

---

## The decision

$$\boxed{\textbf{PARAMETERIZED\_META\_OPERATOR\_CLOSED}}$$

## Repair audit

| Required repair | Delivered | Status |
|---|---|---|
| **1. Full metric transfer** | Weighted metric $d_{\mathbb M}=\sup_\iota(\alpha d_K+\beta d_C+\gamma d_R)$ (PM-1, complete, equivalent to the Phase-11 metric up to constants; Route-A scope for TV-compactness declared). Three separate lemmas: probability $d_K\le\tfrac12\bar H\,d_{\mathrm{desc}}$ (Hoffman, the valid part of OM-2 — PM-2a); confidence — universal Lipschitz form refuted by the audit's witness, valid under declared (CONF-ALIGN), and for the canonical pair exactly $d_C=\delta_N$, schedule-controlled (PM-2b); rung — explicit margin condition ($d_{\mathrm{desc}}<\Delta_r\Rightarrow r_1=r_2$ for description-assigned rungs; for the canonical pair $\Delta_r=\infty$ on the all-fibers-observed event, failure priced as an event, cancellation retracted — PM-2c). Combined transfer $d_{\mathbb M}\le\alpha C_K(\eta+\rho)+\beta\delta_N+\gamma\cdot\mathbf 1\{\mathrm{Miss}\}\to0$ with each coordinate carried by its own lemma (PM-3) | **met** |
| **2. Zero-fiber finite-sample repair** | Missing-fiber event defined over relevant contexts; $\Pr(\mathrm{Miss}_N)\le\sum_{c\in C^+}(1-\pi(c))^N\le|C^+|e^{-N\pi_{\min}}$ (PM-4); finite-$N$ theorem restated with the explicit probability term and an unconditional expectation bound using the trivial metric bound on bad events — **no uniform finite-$N$ claim without the event** (PM-5(i)); a.s. consistency re-derived through the event via SLLN on fiber frequencies (PM-5(ii)); per-emission "fiber unobserved" flag added to the ledger | **met** |
| **3. Parameterized meta-learning theorem** | Structural bridge: the canonical operator factors through a finite sufficient statistic $s(H)$ on a compact stratified domain, with explicit stratum-wise Lipschitz evaluation (PM-6) — approximation becomes a deterministic problem against a computable oracle. Abstract family $\{A_\theta\}_{\theta\in\Theta}$ reading $s(H)$ with shared canonical postprocessing (coherence inherited, membership in $\mathbb M$ proved — PM-7); declared capacity assumption (P-CAP) on a compact finite-dimensional domain; uniform approximation theorem $\inf_\theta\sup_H d_{\mathbb M}(A_\theta(H),A_\phi(H))\le\varepsilon'$ with explicit constant (PM-8); total-error theorem for a trained $\theta_N$: approximation $+$ optimization $+$ statistical $+$ missing-fiber terms, each separately attributed, $\to0$ under the schedules (PM-9). Four tiers separated — existence / identification (with $\theta$-gauge explicitly non-identified and irrelevant) / statistical learning / approximation learning — no tier used for another; notation aligned ($A_\phi$ = canonical operator; $M^\dagger$ = statistical target; the two convergences composed, never merged) | **met** |

## Residual open items (non-blocking)

Route-B stability class for continuous scalar outcomes (declared scope boundary; Route A covers the ranking/finite-outcome interface); data-driven lower bounds for $\pi_{\min}$ with tagged confidence (exact-sum form always available); sharp Hoffman constants; capacity-verification templates for concrete families (P-CAP is a checkable declaration per instantiation — deliberately outside this program's scope).

## Closing

The bridge the audits kept demanding now exists as three theorems with nothing hidden between them: descriptions reach the full operator metric coordinate by coordinate — Hoffman for the polytopes, the schedule for the confidence, an explicit margin-or-event analysis for the rungs; finite samples pay for unseen fibers in an explicit probability term rather than a silent exemption; and a parameterized family reaches the canonical operator through a finite compact statistic against a computable oracle, so that its total distance to the identified target splits into approximation, optimization, statistical, and coverage terms — each declared, each priced, each vanishing under its own schedule. The parameterized meta-operator is closed.

**Verdict: `PARAMETERIZED_META_OPERATOR_CLOSED`.**
