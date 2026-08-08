# Population Decision Information (Part III–IV)

> **Status:** Phase-7, 2026-08-03. Frozen corpus cited, not modified. New results carry **DE-H** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Question: can *previously observed members of the family* supply the decision primitive of Part II — and does the answer yield a single law, a set of laws, or neither?

---

## 1. The two channels — feasibility vs frequency

An archive of previous members $f_1,\dots,f_n$ can act on the problem through exactly two logically disjoint channels.

**Theorem DE-H1 (channel separation). [proved]**
(i) *Feasibility channel (identification-level, frozen).* The archive constrains the window system $\mathbb T$ on the covered set $U$ — under the frozen closure assumptions (F17; Phase-6 gluing/union for irregular/imperfect archives) — and thereby changes $I(O)$ itself. This channel uses the archive only as a **set** of witnessed traces; the multiset/frequency structure is discarded; its legitimacy conditions are exactly the frozen ones and are not reopened here.
(ii) *Frequency channel (decision-level, new).* The archive read as a **sample** — how often which members occur — carries information about a *population law* $\pi$ on $\mathcal F$ (or on any functional pushforward of it). This channel cannot change $I(O)$: membership of $f$ in $I(O)$ is a property of $(f,O,\varepsilon)$ quantifier-free in the archive frequencies. It can only weight the inside of $I(O)$.
*Proof of disjointness.* (i) is invariant under duplicating archive members (set semantics); (ii) is precisely what duplication changes. Any archive influence factors as (window constraint, empirical frequency), and the frozen identification map consumes only the first coordinate. $\square$

**Theorem DE-H2 (the ceiling: frequencies never shrink the identified set). [proved]**
Let $\pi$ be any law on $\mathcal F$ and $\lambda$ any noise law consistent with the frozen model ($\operatorname{supp}\lambda(\cdot\mid f)\subseteq\{\tilde y:\|\tilde y-f|_D\|_\infty\le\varepsilon\}$). Then every version of the posterior satisfies $\pi(f\in I(O)\mid \tilde y)=1$ for $\pi\lambda$-a.e. $\tilde y$; and conversely no posterior statement removes any member of $I(O)$ from admissibility: a $\pi$-null subset of $I(O)$ remains admissible — "improbable under an assumed $\pi$" is a graded preference, not an identification statement.
*Proof.* By the support condition, a member $f$ with $f\notin I(\tilde y)$ — i.e. $\|\tilde y-f|_D\|_\infty>\varepsilon$ — has the observed $\tilde y$ outside $\operatorname{supp}\lambda(\cdot\mid f)$, hence carries zero likelihood; disintegration of $\pi\lambda$ then forces every posterior version to put mass $1$ on $\{f:f\in I(\tilde y)\}$. The converse is definitional: $I$ is frozen and $\pi$-free. $\square$

**Corollary DE-H3 (the mandate's prohibition, as a theorem).** Any procedure that reports $I(O)\cap\operatorname{supp}\hat\pi_n$ (or a high-$\hat\pi_n$-probability subset) *as the identified set* asserts an identification claim with no observational basis in the current member — it is exactly the fabrication the information ceiling forbids. Shrinking $I(O)$ via history is legitimate only through channel (i), i.e. as new window/feasibility evidence under a declared closure class. **Frequency may move the *selection*; only evidence moves the *set*.** [proved — restatement of DE-H1/H2]

---

## 2. The likelihood gap: the frozen noise model does not license a posterior

The frozen noise model is adversarial-bounded, not stochastic: it specifies a support, not a law. This has a decision-level consequence usually skipped.

**Theorem DE-H4 (one prior, many posteriors). [proved]**
Fix a prior $\pi$. Let $\Lambda_{\mathrm{adm}}$ be the set of all noise kernels $\lambda$ consistent with the frozen support condition. Then:
(i) each $\lambda\in\Lambda_{\mathrm{adm}}$ yields a (generally different) posterior $\pi(\cdot\mid\tilde y;\lambda)$, all supported in $I(O)$;
(ii) the induced **posterior ambiguity set** $\mathcal P(\tilde y)=\{\pi(\cdot\mid\tilde y;\lambda):\lambda\in\Lambda_{\mathrm{adm}}\}$ can be maximally spread: witness $\mathcal F=\{f_1,f_2\}$, $\pi=(\tfrac12,\tfrac12)$, traces within $2\varepsilon$ so both are consistent with $\tilde y$; take $\lambda_1$ concentrating $f_1$'s noise near the observed $\tilde y$ and $f_2$'s noise elsewhere in its ball, $\lambda_2$ the reverse: the two posteriors approach $\delta_{f_1}$ and $\delta_{f_2}$. Within the frozen model, **the data cannot arbitrate**.
Hence a single posterior requires a declared stochastic noise law — an axiom beyond the frozen model, hereafter **(LIK)**. $\square$

**Declared axioms of the frequency channel** (each optional, each priced):
- **(EXCH)** $f_1,\dots,f_n,f_\beta$ are exchangeable draws from a population; working strengthening **(IID)**: i.i.d. from a law $\pi$.
- **(LIK)** a declared noise law $\lambda$ with the frozen support property (so DE-H2's ceiling survives by construction).

---

## 3. Single law vs set of laws

**Theorem DE-H5 (the minimal extension is naturally set-valued). [proved]**
Under the frequency channel, set-valuedness of the decision object arises from **two independent sources**, each irreducible without a further axiom:
(i) *finite history:* $n$ draws determine $\pi$ only up to a confidence/ambiguity class $\Pi_n$ (Part VIII gives rates; no finite $n$ collapses $\Pi_n$ to a point);
(ii) *likelihood ambiguity:* even with $\pi$ known exactly, DE-H4 yields the posterior set $\mathcal P(\tilde y)$ absent (LIK).
A **single law** therefore has the logical status of a *doubly degenerate case*: it requires (IID with $n\to\infty$, or a declared prior over $\pi$) **and** (LIK). A **set of admissible laws** is what the declared assumptions actually justify at every finite $n$; it degrades gracefully: as assumptions weaken, the set grows, and at the empty-assumption endpoint it becomes all laws supported on $I(O)$ — whose induced decision rule is exactly the frozen minimax rule (endpoint theorem DE-T4, `distribution_shift_and_robustness.md`). "Neither" (no law-like object at all) is the Level-0 choice-function tier of DE-P2(i) — well-defined but structureless and unlearnable. $\square$

**Logical comparison.**

| Object | Assumptions consumed | Status |
|---|---|---|
| single law $\mu$ on $J_Q(O)$ | (IID, $n=\infty$) + (LIK), or declared subjective law | strongest; exceeds what finite data justify — using it is an explicit axiom, not an inference |
| ambiguity set $\mathcal Q_n$ of laws | (EXCH/IID) + confidence level; (LIK) optional (else compose with $\mathcal P(\tilde y)$) | matches the data; the natural object of the extension |
| no law (choice function / pure minimax) | none | frozen fallback; certificate-bearing |

---

## 4. What changes under the three criteria

Given the identified set $J$ (bounded), a declared class $\mathcal Q$ of laws on $J$, and loss $L$:

- **Bayes risk** ($\mathcal Q=\{\mu\}$): $\rho_B(a)=\int L(a,v)\,d\mu$. Fully determined; consumes the maximal tier.
- **Worst-case expected risk** ($\Gamma$-minimax): $\rho_\Gamma(a)=\sup_{\mu\in\mathcal Q}\int L(a,v)\,d\mu$. Interpolates: $\mathcal Q$ singleton $\to$ Bayes; $\mathcal Q=\Delta(J)$ $\to$ frozen minimax exactly (DE-T4).
- **Minimax regret**: $\rho_R(a)=\sup_{v\in J}[L(a,v)-\inf_{b\in\mathcal A}L(b,v)]$; menu-dependent (DE-P4); consumes a declared reference.

**Theorem DE-H6 (comparison). [proved]**
(i) *Coincidence at zero radius:* if $J$ is a singleton, all three criteria (and every monotone completion) select the same actions.
(ii) *Prediction collapse:* if $\mathcal A\supseteq J$ and $L(v,v)=0\le L(a,v)$ (pure prediction), then $\inf_b L(b,v)=0$, so **minimax regret $\equiv$ minimax loss**. The two diverge only under action constraints or state-dependent loss floors.
(iii) *General divergence:* they can select three different actions on one identified set. Witness $J=\{0,1\}$, $\mathcal A=[0,1]$, $L(a,v)=|a-v|+c(v)$ with $c(0)=0,c(1)=5$ (a state-dependent floor):
 minimax loss minimizes $\max(a,\,|a-1|+5)$ $\Rightarrow a^*=1$ (value $5$);
 minimax regret subtracts the floor, minimizes $\max(a,1-a)$ $\Rightarrow a^*=\tfrac12$;
 Bayes with weighting $(w,1-w)$, squared distance term, selects $a^*=1-w$ — sweeping $(0,1)$ as $w$ varies.
(iv) No criterion is privileged by the mathematics: each is a completion choice in the sense of DE-P3(iii), and DE-S4 shows only the minimax tier is canonical. Selecting Bayes vs $\Gamma$-minimax vs regret is itself part of $\Delta$ and must be declared in the ledger. $\square$

---

## 5. The cross-member assumption, stated exactly

**If historical members are to support a decision object, the required assumption is:**
$$\textbf{(EXCH)}\ \ (f_1,\dots,f_n,f_\beta)\ \text{exchangeable}\quad[+\ \textbf{(LIK)}\ \text{for point posteriors};\ +\ \textbf{(COV)}\ \text{(Part VIII) for estimability}].$$
Nothing weaker ties the current member's weighting to historical frequencies at all (Part IX proves the corresponding impossibility); nothing stronger is needed for the set-of-laws formulation. Under (EXCH) alone the object learned is the population law up to the exchangeable/finite-sample ambiguity $\Pi_n$; under (EXCH)+(LIK) the update is honest Bayes conditioning, automatically supported in $I(O)$ (DE-H2) — **population information tilts the choice inside the identified set and can never leave it.**
