# Decision Learnability (Part VIII)

> **Status:** Phase-7, 2026-08-03. Frozen corpus cited, not modified. New results carry **DE-L** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Question: without reference to any computational architecture, can the decision object of Part II be **learned from multiple historical members** — and with what guarantees, under what ceiling?

---

## 1. The four objects (mandated definitions)

**DE-L1 (definitions).**
- **Population-level object.** A law $\pi$ on $\mathcal F$ — or, sufficient for a fixed decision problem, its pushforward $\mu_\psi=\pi\circ\psi^{-1}$ under the decision-relevant functional $\psi$ (e.g. $\psi=e_Q$ for joint prediction: $\mu_Q$ on $\mathbb R^m$; $\psi=\operatorname{sgn}(f(x_a)-f(x_b))$ for ranking: the single number $p_{ab}$). Under finite data the honest population object is an **ambiguity class** $\Pi_n$ (a confidence set of laws), per DE-H5.
- **Current-member identified object.** $I(O)$, equivalently its pushforward $J_Q(O)$ — frozen, untouched.
- **Update/conditioning operation.** Under (EXCH/IID)+(LIK): Bayes conditioning $\mu\mapsto\mu(\cdot\mid O)$, automatically supported in $J_Q(O)$ (DE-H2). Without (LIK): set-conditioning — the posterior ambiguity set $\{\mu(\cdot\mid O;\lambda):\lambda\in\Lambda_{\mathrm{adm}}\}$ (DE-H4), composed with $\Pi_n$. In all cases the update is **restriction-and-reweighting inside the identified set**; it never extends support.
- **Final decision functional.** $D(J_Q(O),\Delta)=\operatorname*{arg\,min}_{a\in\mathcal A}\ \sup_{\mu\in\mathcal Q(O)}\int L(a,v)\,d\mu$, where $\mathcal Q(O)$ is the conditioned ambiguity class ($\Gamma$-minimax; Bayes when $\mathcal Q$ is a singleton; frozen minimax when $\mathcal Q=\Delta(J_Q(O))$ — DE-T4), followed by undominated tie-breaking (DE-O2).

---

## 2. Reduction to the decision-sufficient functional

**Proposition DE-L2. [proved]** For a fixed context $(\mathcal A,L,\mathsf C)$ whose loss factors through $\psi$, the decision depends on the population only through $\mu_\psi$ conditioned on $O$; learning targets may be reduced accordingly (ranking: estimate one Bernoulli parameter, not a law on $\mathcal F$). This is the population-level analogue of DE-S1 and is what makes learnability a *finite-dimensional* estimation problem per decision, even over an infinite-dimensional family. $\square$

---

## 3. Learning theorems

Historical members $f_1,\dots,f_n$ under **(IID)** from $\pi$; each observed only on its own design (frozen archive stratum).

**Theorem DE-L3 (identified history — rates). [conditional on (IID), (COV)]**
Suppose (COV): each historical value $\psi_i=\psi(f_i)$ is identified by member $i$'s own data (frozen conditions: F17 coverage; Thm 7.1 for difference functionals).
(i) *Events (ranking):* $\hat p_n=\tfrac1n\#\{i:\psi_i\in E\}$ satisfies $|\hat p_n-p|\le\sqrt{\ln(2/\delta)/2n}$ with probability $\ge1-\delta$ (Hoeffding).
(ii) *Scalar laws:* $\sup_t|\hat F_n(t)-F(t)|\le\sqrt{\ln(2/\delta)/2n}$ w.p. $\ge1-\delta$ (Dvoretzky–Kiefer–Wolfowitz, Massart constant).
(iii) *Risk transfer:* if $v\mapsto L(a,v)$ has total variation $\le\bar V$ uniformly in $a$, then $\sup_a\big|\int L(a,\cdot)\,d\hat\mu_n-\int L(a,\cdot)\,d\mu\big|\le\bar V\cdot\sup_t|\hat F_n-F|$ (integration by parts) — so plug-in decisions have excess risk $\le2\bar V\sqrt{\ln(2/\delta)/2n}$ under the declared criterion.
Hence the decision object is learnable at rate $n^{-1/2}$, with explicit constants, from identified history. $\square$

**Theorem DE-L4 (censored history — the irreducible width). [conditional on (IID)]**
Without full (COV), member $i$ contributes only the interval $[l_i,u_i]\ni\psi$-event indicator (DE-R5 construction: $l_i=1$ iff member $i$'s own identified set forces the event, $u_i=1$ iff compatible). Then w.p. $\ge1-\delta$:
$$p\ \in\ \Big[\tfrac1n\textstyle\sum_i l_i-\eta_n,\ \ \tfrac1n\sum_i u_i+\eta_n\Big],\qquad \eta_n=\sqrt{\ln(4/\delta)/2n},$$
(two Hoeffding applications, since $\mathbb E l\le p\le\mathbb E u$ by pointwise $l_i\le\chi_i\le u_i$). The interval's *systematic* width $\mathbb E u-\mathbb E l$ equals the population probability that a member's own data leave the event undecided — a **coverage** quantity: $n$-irreducible, reduced only by richer per-member designs. Sampling error $\eta_n\downarrow0$; censoring width does not. This is the exact quantitative form of second-order partial identification (DE-R5). $\square$

*Failure without (IID):* `distribution_shift_and_robustness.md` (DE-T3) proves the corresponding impossibility; nothing in DE-L3/L4 survives an undeclared population change.

---

## 4. The information ceiling, preserved

**Theorem DE-L5 (no fabrication). [proved]**
(i) No decision layer can output an **unconditional** error guarantee below the frozen radius: if a procedure claimed worst-case error $<\tfrac12\omega_{x,D}(2\varepsilon)$ at some configuration, valid without population axioms, it would contradict Theorem 1's exact minimax lower bound (which quantifies over precisely the members of $I(O)$). Guarantees strictly inside the frozen radius exist only as **conditional** statements, tagged with (EXCH/IID)+(LIK)(+coverage/transport declarations), and fail exactly when those axioms fail.
(ii) The learned object influences only the selection *within* $J_Q(O)$ and the conditional risk accounting; by DE-H2 its posterior mass cannot leave the identified set, and by DE-H3 it cannot shrink it. **Historical population information may tilt a decision under ambiguity; it cannot manufacture task-specific information the current observations do not identify.** $\square$

**Learnability verdict.** Yes — the decision object (a monotone completion, generated concretely as a conditioned population law or ambiguity class of laws) is learnable from historical members, at $n^{-1/2}$ with explicit constants where their own values are identified, with an honest $n$-irreducible censoring width where they are not, under the declared axioms (IID/EXCH, LIK, COV) — and it is **not learnable at all without a cross-member axiom** (DE-T3).
