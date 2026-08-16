# Final Compilation Interface (Part VI)

> **Status:** Phase-9, 2026-08-03. The complete mathematical chain from historical tasks to decisions, with every domain and codomain, superseding the interface statements of Phases 8–8.2 wherever they conflict (all conflicts are the audited defects, repaired in this phase's files). New results **ML-I1–I2**. No architectures; no implementation vocabulary.

---

## The chain

$$H_N\ \xrightarrow{\ A_\phi\ }\ M\ \xrightarrow{\ \text{evaluate at }(\kappa(O_*),\gamma_*)\ }\ \Delta_{\mathrm{pop}}\ \ ;\qquad O_*\ \xrightarrow{\ I_\theta\ }\ (\widehat J,\widetilde J,\mathrm{flags})\ \ ;\qquad \big(\widehat J,\ \Delta_{\mathrm{pop}},\ L,\ \tau\big)\ \xrightarrow{\ D_\psi\ }\ \text{decision}.$$

**1. Historical tasks.**
$$H_N=(T_1,\dots,T_N)\in\mathcal T^N\quad(\text{ordered sequence; multiset quotient under declared EXCH; never a set — ML-T3}).$$

**2. Meta-learning operator.**
$$A_\phi:\ \textstyle\bigcup_{N\ge0}\mathcal T^N\ \longrightarrow\ \mathcal M,\qquad \mathcal M=\Big\{M:\ C_\kappa\times\Gamma\to\mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung}\Big\}\ \text{(ML-A1, with its consistency conditions)}.$$
Validity $V_A$: multiset-typed sample use; rung-typed simultaneous coverage at fiber counts $N_c$; zero-fiber vacuous fallback; ceiling (no feasibility feedback); auditability. Output an operator object, never an embedding (ML-A2, derived from the frozen gauge/dimension theorems).

**3. Current-task adaptation — two typed arms (ML-C1).**
$$I_\theta:\ \mathrm{set}(H_N)\times\mathcal O\times\mathcal X^m\to\mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathrm{Flags},\qquad I_\theta(O_*)=(\widehat J,\widetilde J,\mathrm{flags}),\ \ V_I:\ \widetilde J\subseteq J_{Q_*}(O_*)\subseteq\widehat J;$$
$$M_\phi:\ \mathcal T^N\times\mathcal O\ \to\ \mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung},\qquad M_\phi(H_N,S_*)=[A_\phi(H_N)]\big(\kappa(O_*),\gamma_*\big)=\Delta_{\mathrm{pop}},$$
with the $O_*$-dependence factoring through the declared $\kappa$ (ML-K1's stack; corrected kernel-indexed sufficiency; posterior **classes** at rung 4 — ML-K2).

**4. Joint decision object.**
$$\big(\widehat J\ \text{with order projection}\ \widehat\Sigma,\ \widetilde J,\ \Delta_{\mathrm{pop}}\!\upharpoonright_{\mathrm{support}}\big)\in\ \mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathfrak Q(\Omega_m),$$
the support restriction ($P(\widehat\Sigma)=1$; laws on the full $\Omega_m$, events as constraints — DC-R1) applied likelihood-free (DE-H2); emptiness → failure flag, never renormalization.

**5. Decision.**
$$D_\psi:\ \big(\mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathrm{Flags}\big)\times\big(\mathfrak Q(\Omega)\times(0,1]\times\mathrm{Rung}\big)\times\mathrm{Ctx}\times\mathcal T_{\mathrm{tie}}\ \longrightarrow\ \Big(\big(2^{\mathcal A}\setminus\{\emptyset\}\big)\cup\mathrm{FailureReports}\Big)\times\mathrm{Ledger},$$
under the 8.2 $V_D$ with the Phase-9 corrections in force: robustness by the per-law/LP criterion (ML-X3, replacing the interval-separation "iff"); selector-continuity claims restricted to the proved regimes (ML-X1/X2); abstention by criterion-optimality, failure by DC-A4; policy-typed floors and attainment-typed guarantees (8.2 scoping, carried).

---

## Closure theorems

**Theorem ML-I1 (validity of the chain). [proved, given the cited components]**
Under $V_A$, $V_I$, the ML-K stack, and $V_D$: every emission of $D_\psi\circ(I_\theta\times M_\phi)\circ(\mathrm{id},A_\phi)$ is valid at its stated type. The proof is the 8.2 composition argument with each audited step replaced by its Phase-9 repair: sequence/multiset typing (ML-T3/T4) where the set type broke DR-L3-R's applicability; kernel-indexed sufficiency (ML-K1) where the quantifier was wrong; posterior classes (ML-K2) where a singleton was overclaimed; per-law robustness (ML-X3) where the interval test was mis-characterized; the proved selector regimes only (ML-X1/X2). No step cites a refuted statement. $\square$

**Theorem ML-I2 (stop-condition audit). [proved]**
1. *The meta-learning object exists mathematically:* ML-L1 ($M^*_\Pi$ via regular conditional laws; exchangeability suffices at the existence tier), and its estimable surrogate is a well-defined element of $\mathcal M$ for every finite $H_N$ including degenerate fibers (ML-A3, zero-fiber fallback).
2. *Finite-history learning is separated from existence:* the three-tier theorem (ML-L1/L2/L3 + ML-L4) — existence (no data), identification (infinite-history class, point iff a.s. zero censoring), estimation ($\eta_{N_c}$ rates under IID/C-IID-$\kappa$) — with no cross-tier type use.
3. *No task frequency information is destroyed:* $H_N$ is sequence/multiset-typed end to end; the set quotient is applied only inside the feasibility channel, where duplication-invariance is a frozen theorem, not a loss (ML-T3/T4) — the audited archive-type defect is structurally impossible.
4. *The decision operator receives valid learned information:* ML-C1 — the composed pair is exactly $D_\psi$'s typed domain, rung-tagged so that population information is consumed only at the rung its declared assumptions support, with monotone loss-typed degradation at every failure (ML-C2). $\square$
