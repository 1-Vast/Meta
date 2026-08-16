# Conditional Population Repair (Part V)

> **Status:** Phase-9, 2026-08-03. Phases 0–7 frozen and cited. Repairs the audit's conditioning defects in DC-C2/DC-C3 (`../11_final_closure_audit/05_conditioning_audit.md`): (a) the noise-kernel quantifier — sufficiency assumed for one kernel was consumed for all; (b) the rung-4 overclaim — a declared likelihood maps a prior *class* to a posterior *class*, not to a single posterior; (c) the zero-fiber case. New results carry **ML-K** numbers, tagged **[proved] / [conditional] / [declared] / [impossible]**.

---

## 1. The declared objects, with the corrected quantifier

- **($\kappa$-DESIGN) [declared]** A measurable context map $\kappa$ from observable records to a finite $C_\kappa$, computable from any task's record (historical or current; never from unidentified quantities).
- **(C-IID-$\kappa$, query-indexed) [declared]** Conditional on $\kappa=c$, the pairs (record, member) — including the current task's — are IID under the actual data-generating law; the decision target is the query-indexed pushforward $g(f)$ at the current $(Q,g)$.
- **(SUFF-$\kappa$-$\lambda$) [declared, kernel-indexed — the repair of defect (a)]** Sufficiency is a property of the **joint law $\Pi\otimes\lambda$ of (member, record)**, hence indexed by the noise kernel $\lambda$:
$$g(f)\ \perp\ O\ \big|\ \kappa(O)\qquad\text{under }\Pi\otimes\lambda.$$
Two admissible declarations: **(SUFF-$\kappa$-$\lambda^*$)** — sufficiency for a *declared actual kernel* $\lambda^*$; or **(SUFF-$\kappa$-U)** — sufficiency uniformly over all $\lambda\in\Lambda_{\mathrm{adm}}$. The audit's witness is adopted as the reason the index is load-bearing: with constant $\kappa$, an uninformative admissible kernel satisfies sufficiency while a revealing admissible kernel violates it — so per-kernel sufficiency **does not transfer** across kernels, and the Phase-8.2 statement quantifying over all $\lambda$ from a single sufficiency assumption was false.

---

## 2. The corrected conditioning theorem

**Theorem ML-K1 (when population information can condition the current task). [conditional on the declared stack]**
Assume ($\kappa$-DESIGN) + (C-IID-$\kappa$, query-indexed) + $N_{\kappa(O_*)}\ge1$. Then:
(i) Under **(SUFF-$\kappa$-$\lambda^*$)** with declared actual kernel $\lambda^*$: for the actual law,
$$P\big(g(f_*)\in\cdot\ \big|\ O_*,Q\big)\ =\ P\big(g(f_*)\in\cdot\ \big|\ \kappa(O_*),Q\big)\ =\ \pi_g\big(\cdot\ \big|\ \kappa(O_*)\big),$$
and the fiber estimate (ML-L3, multiset-typed count $N_{\kappa(O_*)}$) is a valid $1-\delta$ confidence class **for the current-observation-conditioned object** — with the support restriction to the identified object applied on top, likelihood-free (DE-H2).
(ii) Under **(SUFF-$\kappa$-U)**: the same identity holds simultaneously for every admissible kernel, so the conclusion is kernel-robust.
(iii) Under neither: the only kernel-free conditional statement is the support restriction; the population class may be consumed **only at its marginal rung** (DC-C3 as corrected below). The conditioning is a **declared-assumption compilation target, not a derivation from identification** — the audit's formulation, adopted verbatim into the contract. $\square$

**Theorem ML-K2 (posterior classes — the repair of defect (b)). [proved]**
A declared likelihood $\lambda^*$ defines the posterior **map**, not a posterior: applied to the prior ambiguity class $\widehat{\mathcal Q}$ it yields the **posterior class**
$$\widehat{\mathcal Q}^{\,\mathrm{post}}\ =\ \big\{\,P_0(\cdot\mid O_*;\lambda^*)\ :\ P_0\in\widehat{\mathcal Q},\ P_0\text{-marginal likelihood}>0\,\big\},$$
which is a singleton iff $\widehat{\mathcal Q}$ is (up to likelihood-null degeneracies). Rung 4 of the conditioning ladder is therefore re-labeled **"posterior class"**; the 8.2 phrase "single posterior" is retracted and holds only in the doubly degenerate case (singleton prior class + declared kernel) — which is exactly Phase-7's DE-H5 statement, now enforced at the ladder's top rung too. $\square$

**The corrected ladder (replacing DC-C3's table). [proved/conditional as marked]**
| Declared | Valid conditioned object | Consumption type |
|---|---|---|
| nothing / $N_{\kappa(O_*)}=0$ | all laws on the identified support | rung 1: frozen minimax endpoint (vacuous fallback — defect (c) closed) |
| (C-IID-$\kappa$) | fiber class, **marginal-typed** | rung 2: only decisions declared insensitive to residual-record dependence |
| + (SUFF-$\kappa$-$\lambda^*$) or (SUFF-$\kappa$-U) | fiber class **as the $O_*$-conditional** (ML-K1) | rung 3: full conditional machinery |
| + declared $\lambda^*$ acting on the class | **posterior class** $\widehat{\mathcal Q}^{\,\mathrm{post}}$ (ML-K2) | rung 4: Bayes-class tier; singleton only if prior class singleton |

**Impossibility floor (carried). [impossible]** Below rung 3 the gap between the $\kappa$-conditional and the full-$O_*$-conditional is not estimable within the frozen noise model (DE-H4: one prior, many posteriors); no construction removes this without a sufficiency or likelihood declaration.

---

## 3. Falsifiability note (declared assumptions are still testable)

(SUFF-$\kappa$-$\lambda^*$) has observable consequences: under it, within a fiber the conditional record law factors from the decision target; systematic dependence between residual record features and *identified* historical targets inside a fiber refutes the declaration (a per-fiber independence audit, run on the subsample of historical tasks whose own data identify the target). The contract requires declaring it and echoes it; refutation downgrades the operator to rung 2 — a monotone, honest degradation (ML-C2).
