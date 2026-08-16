# Final Composition (§6)

> **Status:** Phase-10, 2026-08-03. The complete typed chain, every arrow with valid domain and codomain, on the repaired spaces. New results **LC-18–LC-19**.

---

## The chain

$$H_N\ \xrightarrow{\ A_\phi\ }\ M\in\mathbb M\ \xrightarrow{\ \mathrm{ev}_{(\kappa(O_*),\,Q_*,\,\gamma_*)}\ }\ \Delta_{\mathrm{pop}}\qquad\Big/\qquad (O_*,S_*,Q_*)\ \xrightarrow{\ I_\theta\ }\ (\widehat J,\widetilde J,\mathrm{flags})\qquad\Longrightarrow\qquad J_Q\text{-object}\ \xrightarrow{\ D_\psi\ }\ \text{decision}$$

**1. Marked historical tasks.** The theory lives on $\mathbb T^\bullet=\mathbb T\times\mathbb Y$ with a declared marked law $\Pi^\bullet$ (LC-1/2); the **sample is observable**: $H_N=(T_1,\dots,T_N)\in\mathbb T^N$, the projections of $N$ marked draws — marks are never observed, and no arrow consumes them (they exist so that latent conditionals are well-typed; LC-3). $H_N$ is sequence/multiset-typed (MC-2/3/4 carried).

**2. Meta-learning operator.**
$$A_\phi:\ \bigcup_{N\ge0}\mathbb T^N\ \longrightarrow\ (\mathbb M,\ d_{\mathbb M},\ \text{evaluation }\sigma\text{-algebra}),$$
measurable and total (LC-10); validity $V_A$ + index-complete coverage (LC-9, LC-15); codomain the query-indexed operator space (LC-5/8) — never an embedding, latent vector, or task ID (MC-7 + LC-8 §3).

**3. Current support and observation.** $T_*=(O_*,S_*,Q_*,\gamma_*)$. Two typed arms:
$$I_\theta:\ \mathrm{set}(H_N)\times\mathcal O\times\mathcal X^{m_{Q_*}}\to\mathcal K_{m_{Q_*}}\times\mathcal K^{\mathrm{fin}}_{m_{Q_*}}\times\mathrm{Flags},\qquad I_\theta(O_*)=(\widehat J,\widetilde J,\mathrm{flags}),\quad V_I:\ \widetilde J\subseteq J_{Q_*}(O_*)\subseteq\widehat J;$$
$$\mathrm{ev}:\ \mathbb M\times C_\kappa\times\mathcal Q_0\times\Gamma_0\to\textstyle\bigsqcup_Q\big(\mathfrak Q(\Omega_Q)\times(0,1]\times\mathrm{Rung}\big),\qquad \Delta_{\mathrm{pop}}=[A_\phi(H_N)](\kappa(O_*),Q_*,\gamma_*)\ \in\ \mathfrak Q(\Omega_{Q_*})\times(0,1]\times\mathrm{Rung},$$
the query index explicit and preserved (LC-6/7); $O_*$ enters the population arm only through $\kappa$; conditioning valid at rung 3 under the query- and kernel-indexed sufficiency (MC-11, re-anchored on $\mathbb T^\bullet$ by LC-4(ii)).

**4. Identification/population object.**
$$J_Q\text{-object}\ =\ \big(\widehat J\ \text{with}\ \widehat\Sigma_{Q_*},\ \widetilde J,\ \Delta_{\mathrm{pop}}\!\upharpoonright_{\mathrm{support}}\big)\ \in\ \mathcal K_{m_{Q_*}}\times\mathcal K^{\mathrm{fin}}_{m_{Q_*}}\times\big(\mathfrak Q(\Omega_{Q_*})\times(0,1]\times\mathrm{Rung}\big),$$
support restriction likelihood-free (frozen DE-H2; automatic on the marked space by the support coupling, LC-4(iii)); all components typed at the **same** $Q_*$.

**5. Decision.**
$$D_\psi:\ \big(\mathcal K_{m}\times\mathcal K^{\mathrm{fin}}_{m}\times\mathrm{Flags}\big)\times\big(\mathfrak Q(\Omega_{Q})\times(0,1]\times\mathrm{Rung}\big)\times\mathrm{Ctx}\times\mathcal T_{\mathrm{tie}}\ \longrightarrow\ \Big(\big(2^{\mathcal A}\setminus\{\emptyset\}\big)\cup\mathrm{FailureReports}\Big)\times\mathrm{Ledger},$$
under the repaired $V_D$ (rung typing; policy-typed floors; attainment-typed guarantees; per-law LP robustness; criterion-optimal abstention vs failure; H1–H6) — all statements $Q_*$-typed.

---

**Theorem LC-18 (typed composition validity). [proved, given the cited components]**
Every arrow above has the stated domain and codomain; no arrow consumes a mark, erases $Q_*$, crosses the set/multiset channel typing, or consumes a population value above its rung. Under $V_I$, $V_A$ (+LC-15's assumptions where its rates are quoted), the conditioning stack, and $V_D$, every emission of the composed operator is valid at its stated type, with total, monotone, loss-typed degradation. The three audited defects cannot recur: the latent conditional is well-typed by construction (LC-1–3), the query index is load-bearing and preserved by typing (LC-6/7), and the learning claims quoted downstream are operator-level (LC-15) or explicitly tagged pointwise-fallback (LC-15′). $\square$

**Theorem LC-19 (stop-condition audit). [proved]**
1. *Latent probability space repaired:* marked space $\mathbb T^\bullet$ standard Borel; marked vs observable laws distinguished; unanchored latent conditionals excluded by typing; the false existence claim retracted and replaced (LC-1–4, LC-11).
2. *Query index preserved:* $M(\kappa,Q,\gamma)$ with query-typed values and projective coherence; necessity proved by the two-query witness; preservation enforced through every arrow (LC-5–7).
3. *Operator space defined:* $(\mathbb M,d_{\mathbb M})$ complete, evaluation topology separable-metrizable, evaluation σ-algebra, measurable total $A_\phi$; operator (not embedding/latent/ID) by carried and strengthened exclusion (LC-8–10).
4. *Finite-task learning theorem exists:* uniform-over-index deviation bound with VC-type complexity, a.s. operator-metric consistency toward the identified operator, assumption-free coverage fallback, and decision-layer transfer (LC-15/15′/16), with tiers A/B/C kept apart (LC-14). $\square$
