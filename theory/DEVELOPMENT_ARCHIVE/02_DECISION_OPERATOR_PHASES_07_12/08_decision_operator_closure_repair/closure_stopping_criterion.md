# Decision-Operator Closure — Stopping Criterion

> **Status:** Phase-8.2 terminal decision, 2026-08-03. Sources: the five closure files in this directory, the audit `../10_final_theory_audit/` (verdict `DECISION_OPERATOR_REPAIR_FAILED`), and the frozen Phases 0–7 (cited only, unmodified). The 8.1 verdict is superseded; the refuted 8.1 statements (DR-S4-R as stated; the tolerance-abstention clause; the subset-typed ranking class; the marginal-as-conditional step in $V_M$/DR-M1-R) are **retracted** and replaced, not patched in place.

---

## The decision

$$\boxed{\textbf{DECISION\_OPERATOR\_CLOSED}}$$

## Audit against the four stop conditions

| Stop condition | Delivered | Status |
|---|---|---|
| **1. Selection valid** | DC-S1: branch-switching discontinuity with the union-confinement hypothesis **in the statement**, connectedness proof, oscillation $\ge d$; DC-S2 (audit's $(a-t)^2$ witness, adopted) and DC-S3 (fat-bridge ramp) proving the hypothesis indispensable; DC-S4 positive side (Berge strict-quasiconvexity; Michael l.s.c.+convex) distinguishing disconnected branches / connected action space / continuous-selector existence; DC-S5: for finite action spaces confinement is automatic — the ranking discontinuity warning is unconditional. Carried: DR-S1–S3, DR-S5 (argmin set or declared $\tau$; no hidden measures; $\eta$-optimality) | **valid** |
| **2. Abstention loss-consistent** | $a_{\mathrm{abs}}\in\mathcal A$ with declared cost; $R(a_{\mathrm{abs}})=\sup_{v\in J}L(a_{\mathrm{abs}},v)$; DC-A1: abstain **iff** criterion-optimal — the only legitimate threshold rule is the argmin comparison itself (threshold = declared cost vs loss-typed strict robust value, derived); DC-A2: forbidden-threshold witnesses including the audit's $(1,\,c{=}2,\,T{=}\tfrac12)$ triple; DC-A4: exhaustive abstention/failure split — tolerance violations are infeasibility **reports**, never abstentions | **loss-consistent** |
| **3. Joint ranking object valid** | DC-R1: $\Omega_m=S_m$ with laws on the **full** space and queried events (orders, pairwise, top-$k$) as interval **constraints** — the subset-law pathology structurally impossible; pairwise probabilities typed as marginals; DC-R2: pairwise consistency = linear-ordering-polytope membership (dicycle-complete for $m\le5$, scoped for $m\ge6$; LP-checkable always); DC-R3: Kendall-type decomposable losses are marginal-sufficient (linearity); DC-R4: equal-marginal witness ($\tfrac12\delta_{123}+\tfrac12\delta_{321}$ vs uniform on $S_3$) proving non-decomposable losses need the law; DC-R5: simultaneous risk brackets and robust listwise decisions as LPs over the polytope, coverage inherited from the (correctly consumed) 8.1 event-interval theorem | **valid** |
| **4. Conditioning theorem valid** | Route B chosen and proved: DC-C1 states the gap exactly (marginal $\ne$ conditional, with witness); DC-C2: under declared ($\kappa$-DESIGN) + (C-IID-$\kappa$) + (SUFF-$\kappa$), the fiber estimate **is** the current-observation-conditioned class (conditioning statistic a function of $O$; support restriction likelihood-free by the frozen DE-H2); DC-C3: the conditioning ladder — nothing → support-only/minimax endpoint, C-IID-$\kappa$ → marginal-typed, +SUFF-$\kappa$ → conditional, +LIK → posterior — with the DE-H4 impossibility as the floor, and the rung-typing rule forbidding consumption above one's rung; DC-C4: $V_M$ re-typed with rung tags, making the audited error (marginal used as conditional) a type violation, not a runtime risk | **valid** |

**Scoping corrections absorbed:** policy-typed floors ($R_{\mathrm{rand}}$ for "no rule" claims); guarantee row = selected action's own outer risk, $G_{\mathrm{cert}}$ "achieved" only with attainment certificate (else infimum $+\eta$); archive feasibility channel formally order- and multiplicity-invariant ($V_I$(iii)).

**Closure:** DC-I1 proves the composed interface valid with no step resting on a refuted statement; DC-I2 verifies against the audit's own list that every remaining judgment is either a cited theorem or a declared, echoed input slot — a builder inherits no theorem-level choices.

## Residual open items (tightness/refinement only — non-blocking)

Facet description of $P^m_{LO}$ for $m\ge6$ (LP feasibility suffices operationally); sharper-than-union-bound simultaneous constants; falsification protocols for (SUFF-$\kappa$) beyond its declared status (it is testable in principle against held-out fibers, but the contract requires only that it be declared and echoed); tightness rates at order boundaries (validity unaffected).

---

**Verdict: `DECISION_OPERATOR_CLOSED`.**
