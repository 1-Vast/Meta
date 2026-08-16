# Final Composition Theorem (§6)

> **Status:** Phase-9 closure, 2026-08-03. The complete typed chain from historical tasks to decisions, with every domain and codomain. Supersedes prior interface statements wherever they conflict; all conflicts are audited defects repaired in this phase or in `../09_meta_learning_operator_formalization/`. New results **MC-20–MC-21**.

---

## The chain

$$\begin{array}{ccccc} H_N & \xrightarrow{\ A_\phi\ } & M\in\mathbb M & & \\[2pt] & & \big\downarrow{\scriptstyle\ \text{evaluate at }(\kappa(O_*),\,\gamma_*)} & & \\[2pt] O_* & \xrightarrow{\ I_\theta\ } & (\widehat J,\widetilde J,\mathrm{flags})\ \ +\ \ \Delta_{\mathrm{pop}} & \xrightarrow{\ \text{support restriction}\ } & J_Q\text{-object} & \xrightarrow{\ D_\psi\ } & \text{decision} \end{array}$$

**1. Historical tasks.** $H_N=(T_1,\dots,T_N)\in\mathbb T^N$ — ordered sequence; multiset quotient under declared task-exchangeability; the set quotient reserved to the feasibility channel (MC-2/3/4).

**2. Meta operator.** $A_\phi:\bigcup_{N\ge0}\mathbb T^N\to\mathbb M$, $\mathbb M=\{M:C_\kappa\times\Gamma\to\mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung}\}$ with rung-consistency, projective consistency, zero-fiber convention; validity $V_A$ (typed sample use, rung-typed coverage, ceiling, auditability) — MC-6/8; never an embedding, latent vector, or task ID (MC-7).

**3. Current observations + identification.** $I_\theta:\mathrm{set}(H_N)\times\mathcal O\times\mathcal X^m\to\mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathrm{Flags}$, $V_I:\ \widetilde J\subseteq J_{Q_*}(O_*)\subseteq\widehat J$ (order projection outer), declared-closure-class semantics; population evaluation $\Delta_{\mathrm{pop}}=[A_\phi(H_N)](\kappa(O_*),\gamma_*)\in\mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung}$, with $O_*$-dependence factoring through the declared $\kappa$ (MC-11's stack).

**4. Joint decision object.**
$$J_Q\text{-object}\ =\ \big(\widehat J\ \text{with}\ \widehat\Sigma,\ \ \widetilde J,\ \ \Delta_{\mathrm{pop}}\!\upharpoonright_{\mathrm{support}}\big)\ \in\ \mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\big(\mathfrak Q(\Omega_m)\times(0,1]\times\mathrm{Rung}\big),$$
support restriction likelihood-free (DE-H2), on the full outcome space $\Omega_m=S_m$ with events as constraints (DC-R1); emptiness → failure flag; material mass outside $g(\widehat J)$ → SUFF-$\kappa$ falsification signal (MC-13).

**5. Decision.**
$$D_\psi:\ \big(\mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathrm{Flags}\big)\times\big(\mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung}\big)\times\mathrm{Ctx}\times\mathcal T_{\mathrm{tie}}\ \longrightarrow\ \Big(\big(2^{\mathcal A}\!\setminus\!\{\emptyset\}\big)\cup\mathrm{FailureReports}\Big)\times\mathrm{Ledger},$$
under the repaired $V_D$: rung check (no marginal-to-conditional consumption, MC-14); guarantee row = selected action's own outer risk (attainment-typed $G_{\mathrm{cert}}$); policy-typed witness floors ($R_{\mathrm{rand}}$ for "no rule" claims); argmin-set or declared-$\tau$ selection, selector-continuity claims restricted to the proved regimes (ML-X1/X2); per-law LP robustness criterion (ML-X3); abstention iff criterion-optimal, failure per DC-A4; H1–H6.

---

## The theorems

**Theorem MC-20 (final composition — validity). [proved, given the cited components]**
Under $V_A$, $V_I$, the §4 conditioning stack, and $V_D$, every emission of
$$\mathbb D\ =\ D_\psi\ \circ\ \big(I_\theta\ \times\ \mathrm{eval}_{(\kappa(O_*),\gamma_*)}\!\circ A_\phi\big)$$
is valid at its stated type: floors by witness membership + monotonicity (policy-typed); guarantees by outerness; conditional rows by MC-11/12 at rung 3 (the conditional identity **proved** from the declared kernel-indexed sufficiency — no leap), posterior classes at rung 4, marginal-typed consumption only at rung 2, vacuous fallback at rung 1; ranking statements on the full-space polytope with LP-decidable robustness; selection and abstention by the corrected theorems. Degradation is total, monotone, loss-typed, and rung-typed (MC-10). No step cites a refuted statement: every defect of the three audit rounds (`../09_phase8_audit/`, `../10_final_theory_audit/`, `../11_final_closure_audit/`) is repaired at the type level or replaced by a corrected theorem, and the repair map is recorded in the respective phase logs. $\square$

**Theorem MC-21 (stop-condition audit). [proved]**
1. *Meta-learning object mathematically defined:* $\mathbb T$ (MC-1), $H_N$ sequence-typed with the channel quotients (MC-2/3/4), $\mathbb M$ and $A_\phi$ with domain, codomain, validity conditions, a constructive witness, and derived exclusions of embeddings/latents/IDs (MC-6/7/8).
2. *Finite-history learning separated from existence:* the three-tier theorem — unconditional existence (MC-16), partial identification with the exact iff (MC-17), estimation under declared task-IID/C-IID-$\kappa$ + finite complexity + concentration (MC-18) — with tier independence certified (MC-19) and the task-level/within-task firewall enforced (§5.0).
3. *Current-task conditioning valid:* Route B proved (MC-11: the conditional-independence identity for the declared kernel; MC-12: estimability; MC-13: support-consistency falsification; MC-14: posterior classes, ladder, and the impossibility floor making any marginal-to-conditional leap a type violation).
4. *Decision operator receives valid learned information:* MC-9 (decision sufficiency of the composed pair), MC-10 (totality), MC-20 (chain validity). $\square$
