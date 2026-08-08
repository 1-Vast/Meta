# Decision-Operator Realizability — Stopping Criterion

> **Status:** Phase-8 terminal decision, 2026-08-03. Sources: the five Phase-8 files (DR-F, DR-S, DR-J, DR-L, DR-M numbers) and the contract `THEORY_TO_MODEL_INTERFACE.md`. Phases 0–7 were consumed as frozen citations only; no identification or decision theorem was modified, reinterpreted, or extended — Phase 8 only *typed* them into a realizable contract, with two explicit bookkeeping tightenings flagged as such (the DE-O2 tie-break $\mu_0$ reclassified as a declared $\tau$; the ledger extended by floor and tolerance rows).

---

## The decision

$$\boxed{\textbf{DECISION\_OPERATOR\_COMPLETE}}$$

## Audit against the stop condition

| Required | Delivered | Status |
|---|---|---|
| **Robust floor theorem** | $R_{\mathrm{set}}(J,\mathcal A,L)=\inf_a\sup_{v\in J}L(a,v)$ with: reduction to the Chebyshev radius / frozen minimax identity under absolute error (DR-F1(i)); sup-loss, ranking (pairwise and listwise, with $R_{\mathrm{rand}}$), and structured-action specializations (DR-F1(ii)–(iv)); the unbeatability theorem for arbitrary history-dependent randomized rules (DR-F2) with the exact characterization of the only floor-moving routes (DR-F3); and the outer-semantics approximation calculus — conservative validity, Lipschitz tightness, prohibition of inner envelopes (DR-F4) | **exists, proved** |
| **Honest selection theorem** | set-valued argmin as the canonical object (DR-S1); uniqueness under strict quasiconvexity — squared-loss prediction never needs a tie-break (DR-S2); the equivariance obstruction making Option A / declared-$\tau$ Option B an exhaustive forced dichotomy with no hidden measures (DR-S3); tie-boundary discontinuity warning (DR-S4); and the $\eta$-argmin theorem giving approximators a valid "$\eta$-optimal" claim (DR-S5) | **exists, proved** |
| **Joint-object learnability condition** | the quotient lattice with strictness proofs (DR-J1/J2); the listwise object $\Sigma(J)$ with the pairwise-over-admission witness and its floor cost $2$ vs $3$ (DR-J3); minimal decision-sufficient object for the ranking application (DR-J4); and the composed learnability condition — per-task identification gate $\wedge$ declared (fiber-relative) exchangeability gate $\wedge$ declared transport gate, with the explicit interval class, Hoeffding margins, and the three-width decomposition (DR-L3/L4) | **exists, proved conditional on declared axioms, with the gate-failure impossibilities frozen-cited** |
| **Meta-learning interface mathematically defined** | exactly three objects with domains, codomains, and validity predicates — $I_\theta:\mathcal W\times\mathcal O\times\mathcal X^m\to\mathcal K_m\times\mathrm{Flags}$ ($V_I$: outer envelope); $M_\phi:\mathcal W\times\text{tags}\to\mathfrak Q_g\times(0,1]$ ($V_M$: outer confidence class); $D_\psi:(\mathcal K_m\times\mathrm{Flags})\times(\mathfrak Q_g\times(0,1])\times\mathrm{Ctx}\to(2^{\mathcal A}\!\setminus\!\{\emptyset\}\cup\{\mathrm{abstain}\})\times\mathrm{Ledger}$ ($V_D$: H1–H6 + floor + honest selection) — with the honest-composition theorem, the end-to-end tightness calculus, and the proof that the three-way factorization is minimal (DR-M1/M2 + impossibility guard). No architectures, parameterizations, or training losses appear | **exists, defined** |

The compilation contract (`THEORY_TO_MODEL_INTERFACE.md`) instantiates all four in the mandated Input / Learned / Forbidden / Output / Failure form, every clause theorem-cited, with falsification checks extending the frozen protocol suite.

## Open items, each certified non-blocking

| Open item | Why it does not block |
|---|---|
| Tightness rates for discrete/order losses (DR-F4 remark: no Lipschitz transfer at order boundaries) | Validity is unconditional (outer semantics); only the *sharpness* of conservative floors near boundaries is unquantified. The contract's abstain/flag clause 6 covers the boundary region honestly. |
| Exact characterization of when $\Sigma(J)$ is computable from finitely many declared functional intervals (beyond the convex/linear regime) | $\Sigma^{\mathrm{pair}}$ is an always-available outer proxy (DR-J3(ii)); tightening is a quality objective by contract. |
| Optimal fiber-smoothing under a *declared* modulus on $C$ (rates interpolating $n_{c}\to n$) | The two endpoints are settled (exact fiber counts; impossibility without declared structure, CI-A5(iii)/DR-L2(ii)); the interpolation is a refinement inside an already-declared assumption. |
| Sharp constants in the listwise interval class (union-bound over $S_m$ vs queried subsets) | Conservative constants are valid; sharpening tightens, never corrects. |

None prevents constructing, certifying, auditing, or falsifying a realization; each is a tightness refinement inside clauses the contract already handles conservatively.

---

## Closing

The program now terminates in a typed contract: identification produces outer sets, population learning produces outer classes, decision produces $\eta$-optimal actions with a ledger — and the three arrows cannot be rewired without contradicting a frozen theorem. What a future system may learn is exactly two outer approximations; what it may never do is exactly six theorems; what it must say when it cannot act is exactly six clauses. The floor is loss-typed, the selection is honest, the learnability is gated, the interface is closed.

**Verdict: `DECISION_OPERATOR_COMPLETE`.**
