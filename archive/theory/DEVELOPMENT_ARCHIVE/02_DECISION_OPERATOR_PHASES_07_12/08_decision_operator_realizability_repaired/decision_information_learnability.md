# Learnability of Decision Information — REPAIRED (Part IV)

> **Status:** Phase-8.1, 2026-08-03. Supersedes `../08_decision_operator_realizability/decision_information_learnability.md`. Repairs audit targets **T5** (exchangeability does not concentrate) and **T6** (simultaneous coverage; over-broad "iff"). The typing rule, DR-L1, and the transport/ambiguity rungs are carried; the assumption ladder and the learning theorem are restated as **DR-L2-R / DR-L3-R / DR-L4-R**.

---

## 1. Carried (audit: pass or untouched)

- Typing rule: current-task information unconditional; population information conditional, confined to weighting the inside of the identified object (DE-H2/H3).
- **DR-L1 [carried]:** $\Delta_{\mathrm{population}}=(\widehat{\mathcal Q},\,1-\delta,\,\text{axiom tags})$, an outer confidence class of laws on the decision-relevant pushforward space.
- Transport rung (DE-T2 formulas; DE-T3 impossibility without any declared class) and ambiguity-set semantics (outer validity of $\widehat{\mathcal Q}$) — carried.

---

## 2. DR-L2-R: the assumption ladder, with concentration made explicit

**(i) Exchangeability (EXCH) — retained, correctly scoped. [proved insufficiency]**
(EXCH) makes the predictive problem well-posed (de Finetti: a mixture of IID; empirical frequencies converge to the **directing-measure-conditional** parameter, not the marginal). It does **not** yield concentration around population marginals: audit counterexample of record — all task indicators equal to one shared Bernoulli $Z(\tfrac12)$: the sequence is exchangeable, every empirical frequency is $0$ or $1$, and never approaches $\tfrac12$. Consequently **no rate in this program is claimed under (EXCH) alone**; under (EXCH), the honest statements are mixture-conditional and the confidence class must target the conditional law — or the rung must be upgraded.

**(ii) IID.** $f_1,\dots,f_n\stackrel{\text{iid}}\sim\pi$, and $f_\beta\sim\pi$ independently (or via a declared transport class to $\pi$). This is the rung under which Hoeffding/DKW rates are actually proved.

**(iii) Conditional IID (C-IID).** Tasks carry observed labels $c\in C$; conditional on $c$, tasks are IID from $\pi(\cdot\mid c)$, and $f_\beta\sim\pi(\cdot\mid c_\beta)$. All rates hold with $n\to n_{c_\beta}$ (fiber count). Cross-fiber borrowing still requires a declared modulus on $C$ (CI-A5(iii) lift — carried).

**(iv) Declared concentration condition.** Where neither (ii) nor (iii) is defensible, a separately declared concentration inequality for the actual dependence structure (e.g. mixing coefficients with their constants) may substitute; it must be declared with its constants — it is an assumption, not a derivation.

**(v) Complexity assumption.** Simultaneous statements over a family of targets require either a **finite declared index family** (union bound — used below) or a **declared uniform-convergence condition** (VC/covering bound with its constants) for infinite families. Nothing uniform is claimed without one of these.

---

## 3. DR-L3-R: the learning theorem, repaired

**Theorem DR-L3-R (simultaneous, correctly-assumed coverage). [conditional on the declared rungs]**
Fix the decision-relevant map $g$ and a **finite** queried order/event family $S$ (for listwise: $S\subseteq S_m$, $|S|\le m!$; for pairwise: the queried pairs). Assume:
**(a) per-task gate** — each historical task yields certified bounds $l_i(\sigma)\le\mathbf 1\{g(f_i)=\sigma\}\le u_i(\sigma)$ from its own identified object (frozen conditions; vacuous $[0,1]$ allowed);
**(b) (IID)** — or **(C-IID)** with $n\to n_{c_\beta}$;
**(c) declared transport** of TV-radius $\rho$ (possibly $0$).
Then, with $\displaystyle\eta_n=\sqrt{\frac{\ln(4|S|/\delta)}{2n}}$, with probability $\ge1-\delta$ **simultaneously for all $\sigma\in S$**:
$$p(\sigma)\ \in\ \Big[\ \tfrac1n\textstyle\sum_i l_i(\sigma)-\eta_n-\rho\ ,\ \ \tfrac1n\sum_i u_i(\sigma)+\eta_n+\rho\ \Big].$$
*Proof.* For each $\sigma$: $\mathbb E\,l\le p_h(\sigma)\le\mathbb E\,u$ (pointwise bounds), and each empirical mean is an IID average of $[0,1]$ variables — two one-sided Hoeffding events of probability $\le\delta/(2|S|)$ each; union over $2|S|$ events; then $|p(\sigma)-p_h(\sigma)|\le\rho$ by the transport declaration. $\square$
The class $\widehat{\mathcal Q}$ is the (convex, compact) polytope of laws on $S$ consistent with these intervals; its **simultaneous** coverage is exactly the displayed event. For scalar CDF targets, DKW already gives uniform-in-$t$ coverage at $\eta_n=\sqrt{\ln(2/\delta)/2n}$ without a union bound — the two routes are both available and both declared.

**Scoped necessity (replacing the over-broad "iff"). [proved, scoped]**
Within the **frozen distribution-free information model** (no declared likelihood for historical observations), gate (a) is also necessary: a task whose own data leave $g(f_i)$ entirely unrestricted contributes exactly the vacuous interval, and no distribution-free functional of its data can do better (any tighter contribution would distinguish members its data cannot distinguish — the frozen indistinguishability). **Outside** that model — under a *declared* per-task measurement/likelihood model — population aggregates can be identifiable without per-task identification (deconvolution-type inference); such a declaration is a new assumption with its own conditions and is merely permitted, not provided, by this contract. The gates are therefore: sufficient as stated; necessary exactly in the distribution-free model.

**DR-L4-R (three widths — carried with the corrected constants).** Width $=$ (systematic censoring width $\tfrac1n\sum(u_i-l_i)$-population analogue: reduced only by per-task coverage) $+$ ($2\eta_{n_{\mathrm{eff}}}$ with the **union-bound constant** $\ln(4|S|/\delta)$: reduced by fiber-relative task count) $+$ ($2\rho$: reduced only by tighter declared transport). Failure of any gate degrades to the vacuous class → the DE-T4 minimax endpoint: honest, conservative, never false. $\square$

---

## 4. Repaired summary

$$\boxed{\begin{array}{c}\text{Rates require (IID) or (C-IID) or a declared concentration condition — never bare exchangeability (shared-Bernoulli witness).}\\ \text{Simultaneous listwise coverage holds with the union-bound constant }\ln(4|S|/\delta)\text{ over the declared finite family (or a declared uniform bound).}\\ \text{Per-task identification is necessary only in the distribution-free model; a declared likelihood model may buy aggregates — as a declared assumption.}\end{array}}$$
