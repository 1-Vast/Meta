# Full Metric Transfer (§1)

> **Status:** Phase-13 (parameterized meta-operator closure), 2026-08-03. Phases 0–11 unmodified; the operator is not redesigned. Audit of record: `../12_final_operator_audit/FINAL_VERDICT.md` (`OPERATOR_LEARNABILITY_INCOMPLETE`). **Retracted:** OM-2's universal form "$d_{\mathbb M}\le\tfrac12\bar H\,d_{\mathrm{desc}}$" — the audit's witness (identical polytopes and rungs, confidences $0$ vs $1$: $d_{\mathrm{desc}}=0$, $d_{\mathbb M}=1$) is valid; the Hoffman argument transfers **only** the probability coordinate. New results carry **PM-** numbers, tagged **[proved] / [conditional] / [declared] / [retracted]**.

---

## 1. The weighted operator metric

**Definition PM-1.** With declared weights $\alpha,\beta,\gamma>0$ (echoed; default $\alpha=\beta=\gamma=1$), per index:
$$d_{\mathbb V}\big((K_1,c_1,r_1),(K_2,c_2,r_2)\big)\ =\ \alpha\,\underbrace{d_H^{TV}(K_1,K_2)}_{d_K}\ +\ \beta\,\underbrace{|c_1-c_2|}_{d_C}\ +\ \gamma\,\underbrace{\mathbf 1\{r_1\ne r_2\}}_{d_R},\qquad d_{\mathbb M}=\sup_{\iota\in\mathcal I}d_{\mathbb V}.$$
A genuine complete metric (each summand a metric on its complete factor; sum of metrics; bounded by $\alpha+\beta+\gamma$); equivalent to Phase-11's max-product metric within constant factors, so all Phase-11 completeness/measurability results carry. **Scope note (from the audit, adopted as declared scope):** the TV-compactness of $\Delta(\Omega_Q)$ — hence the hyperspace completeness — holds for the Route-A finite-outcome class; continuous scalar outcome spaces are *not* covered by Route A and require Route B's declared stability class.

## 2. The three coordinate lemmas — separately, as mandated

**Lemma PM-2a (probability coordinate). [proved]**
On the Route-A class (bounded $\bar n,\bar e$; uniform Hoffman constant $\bar H$; both polytopes nonempty; same constraint pattern):
$$d_K(K_1,K_2)\ \le\ C_K\,d_{\mathrm{desc}}(M_1,M_2),\qquad C_K=\tfrac12\bar H.$$
This is exactly the valid content of Phase-11's OM-2 (the Hoffman transfer), unchanged.

**Lemma PM-2b (confidence coordinate). [conditional on a declared alignment; universal form refuted]**
No universal bound $d_C\le L_C\,d_{\mathrm{desc}}$ exists — the audit's witness refutes it, since confidence is not a function of the description. The valid statements:
(i) **(CONF-ALIGN) [declared]** If both operators' confidence coordinates are assigned by one declared map $c=\psi(\text{description})$ with $\psi$ $L_C$-Lipschitz w.r.t. $d_{\mathrm{desc}}$, then $d_C\le L_C\,d_{\mathrm{desc}}$. **[proved given the declaration]**
(ii) **Canonical pair (the one the learnability theorem compares):** the estimator carries $c=1-\delta_N$ by the declared schedule, the target carries $c=1$ (OM-4); hence $d_C=\delta_N$ exactly — controlled by the **schedule** (S5: $\delta_N\downarrow0$), not by $d_{\mathrm{desc}}$. The transfer theorem below uses (ii); (i) is available for comparing two schedule-aligned estimators. $\square$

**Lemma PM-2c (rung coordinate — explicit margin condition). [proved]**
(i) *General form:* if both rung coordinates are assigned by one map $r=\rho(\text{description})$ that is constant on the $d_{\mathrm{desc}}$-ball of declared radius $\Delta_r$ around $M_2$'s description (the **margin condition**), then $d_{\mathrm{desc}}(M_1,M_2)<\Delta_r\ \Rightarrow\ r_1=r_2$, i.e. $d_R=0$.
(ii) *Canonical pair:* both estimator and target take $r=r^{\mathrm{decl}}(\iota)$ — a declared function of the index — **except** at zero fibers, where the estimator degrades to rung 1 while a positive-mass context's target may sit at rung 2/3 (the audit's finite-$N$ defect). So for the canonical pair the margin condition holds with $\Delta_r=\infty$ **on the all-relevant-fibers-observed event** $\{N_c\ge1\ \forall c:\Pi_{\mathrm{obs}}(c)>0\}$, and fails only on its complement — which is a *probability event*, priced in `zero_fiber_finite_sample.md`, never silently cancelled (Phase-11's OM-5-pre, which cancelled the rung coordinate at every index, is **corrected accordingly: [retracted]** in its unconditional form, restated as conditional on the fiber event). $\square$

## 3. The combined transfer

**Theorem PM-3 (full-metric transfer, correctly scoped). [proved from PM-2a–c]**
For the canonical estimator–target pair, on the event {all relevant fibers observed} $\cap$ {uniform endpoint deviation $\le\eta$}:
$$d_{\mathbb M}\big(\widehat M_N,\ M^\dagger\big)\ \le\ \alpha\,C_K\,(\eta+\rho)\ +\ \beta\,\delta_N\ +\ \gamma\cdot0\ \longrightarrow\ 0$$
as $\eta\to0$, $\delta_N\to0$ (schedule), $\rho=0$; off the fiber event the trivial bound $d_{\mathbb M}\le\alpha+\beta+\gamma$ applies and the event's probability is charged explicitly (PM-4). For general schedule-aligned pairs, replace $\beta\delta_N$ by $\beta L_C\,d_{\mathrm{desc}}$ (PM-2b(i)) and the fiber event by the margin condition (PM-2c(i)). **Hence $d_{\mathbb M}\to0$, with each coordinate's convergence supplied by its own lemma — no coordinate is carried by another's proof.** $\square$
