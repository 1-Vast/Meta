# Distribution Shift and Robustness (Part IX)

> **Status:** Phase-7, 2026-08-03. Frozen corpus cited, not modified. New results carry **DE-T** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Setting: historical members drawn from a population law $\pi_h$; the new member from $\pi_c$, possibly $\ne\pi_h$. Identification of the current member is untouched by any of this (it never used $\pi$); what is at risk is exactly the **decision object**.

---

## 1. The minimal cross-population assumption

**Theorem DE-T1 (exchangeability is the minimal identical-treatment axiom). [proved]**
(i) If $(f_1,\dots,f_n,f_\beta)$ are **exchangeable** — the weakest symmetry making "the new member is like the old ones" meaningful — then the predictive law of $f_\beta$ given the history is well-defined (for the infinite-exchangeable case, a mixture of i.i.d. laws by de Finetti; the population law itself becomes a parameter with its own posterior), and the Part VIII program applies with $\Pi_n$ enlarged by the mixing uncertainty.
(ii) Exchangeability of the *joint* $(n{+}1)$-tuple is genuinely required: symmetry among the historical members alone ($f_1,\dots,f_n$ exchangeable, $f_\beta$ arbitrary) yields **no** constraint linking $f_\beta$'s weighting to history — the adversarial construction of DE-T3 satisfies it. Assuming "identical distributions" is the strengthening (IID); it must be **declared**, never implied. $\square$

---

## 2. What can be weakened — declared transport classes

**Theorem DE-T2 (quantitative weakenings). [proved]**
Let the decision object be an event probability $p=\pi_c(E)$ (ranking; general bounded losses analogous), with $p_h=\pi_h(E)$ estimable per Part VIII.
(i) **Density-ratio (Λ-transport):** if $\pi_c\ll\pi_h$ with $d\pi_c/d\pi_h\le\Lambda$, then
$$p\;\le\;\Lambda\,p_h,\qquad 1-p\;\le\;\Lambda\,(1-p_h)\;\Rightarrow\; p\in[\,1-\Lambda(1-p_h),\ \Lambda p_h\,]\cap[0,1],$$
and for bounded loss $\mathbb E_{\pi_c}L\le\Lambda\,\mathbb E_{\pi_h}L$.
(ii) **Total-variation ball:** if $d_{TV}(\pi_c,\pi_h)\le\rho$ (sup over events), then $|p-p_h|\le\rho$ and $|\mathbb E_{\pi_c}L-\mathbb E_{\pi_h}L|\le 2\rho\,\|L\|_\infty$ (via the signed-measure decomposition; the factor $2$ with the sup-event convention).
(iii) **Robust ranking threshold:** under (ii) composed with the estimation interval of DE-L3/L4, the ordering is *decision-robust* (Tier 2 of DE-R6) iff
$$\big|\hat p_h-\tfrac12\big|\;>\;\rho+\eta_n+\tfrac12(\hat p^+-\hat p^-),$$
i.e. the historical margin must clear shift radius $+$ sampling error $+$ censoring half-width. Each term is separately declared and separately reducible (by tighter transport, more members, better coverage respectively). $\square$

These are the two canonical weakenings of (IID); moment-band and $f$-divergence classes behave analogously (convex ambiguity classes; DE-T4 covers well-posedness). What **cannot** be weakened away is having *some* declared class: that is the content of the next theorem.

---

## 3. Failure under undeclared shift — the impossibility

**Theorem DE-T3 (unrestricted shift renders history worthless — and committal rules harmful). [proved]**
Consider ranking with both signs admissible (DE-R1(iii)), $0$–$1$ loss, and no declared relation between $\pi_c$ and $\pi_h$ (i.e. $\pi_c$ may be any law supported in the admissible class).
(i) For **every** decision rule $d$ (a function of the history and $O$; randomization allowed), $\sup_{\pi_c}\ \mathbb E\,[\text{loss}]\ \ge\tfrac12$. Precisely: if $d$ is deterministic, the adversary sets $\pi_c=$ the point mass on an admissible member with the sign opposite to $d$'s output — worst-case loss $1$; if $d$ randomizes, then against one of the two admissible point masses $d$ outputs the wrong sign with probability $\ge\tfrac12$, forcing expected loss $\ge\tfrac12$.
(ii) The history-free randomized coin (or abstention, if priced $c\le\tfrac12$) achieves worst-case $\tfrac12$ (resp. $c$). Hence the **value of history under undeclared shift is exactly zero**, and any rule that *commits* to the historical majority sign has worst-case loss $1>\tfrac12$: history is then strictly harmful.
(iii) Consequently no learning theorem of Part VIII survives the deletion of its cross-population axiom; the failure mode is not degraded rates but **adversarial reversal**. $\square$

*Reading.* This is the population-level recurrence of the frozen theme (F4, Thm 9.2: archive data never transport information across an undeclared gap): **no volume of historical members transports decision information across an undeclared population change.** The impossibility recurs one level up, now about laws instead of values.

---

## 4. Robust ambiguity over population laws — possibility and the two endpoints

**Theorem DE-T4 (Γ-minimax well-posedness and the endpoint identities). [proved]**
Let $\mathcal Q$ be a nonempty, convex, weak-*-compact class of laws on the (bounded) identified pushforward $J=J_Q(O)$, $L$ bounded and lower semicontinuous in $a$ on a compact action set $\mathcal A$.
(i) *Existence:* $\rho_\Gamma(a)=\sup_{\mu\in\mathcal Q}\int L(a,\cdot)\,d\mu$ is l.s.c. (supremum of l.s.c. functionals), hence attains its minimum on $\mathcal A$: a robust action exists.
(ii) *Bayes endpoint:* $\mathcal Q=\{\mu\}$ gives Bayes risk.
(iii) *Frozen endpoint:* $\mathcal Q=\Delta(J)$ (all laws on $J$) gives
$$\sup_{\mu\in\Delta(J)}\int L(a,v)\,d\mu(v)\;=\;\sup_{v\in J}L(a,v)$$
(Dirac masses are extreme and achieve the supremum of a bounded measurable integrand), so the robust rule **equals the frozen minimax rule** — with scalar absolute loss, the Chebyshev center with its radius certificate.
(iv) *Monotone interpolation:* $\mathcal Q\subseteq\mathcal Q'$ $\Rightarrow$ $\rho_\Gamma\le\rho_{\Gamma'}$ pointwise; the whole decision extension is one **dial** from the single law (maximal assumptions) to the frozen operator (no assumptions), with the declared class $\mathcal Q$ as the dial position. Shift robustness is the special case where $\mathcal Q$ is the transport image of the historical estimate (DE-T2's classes conditioned on $O$). $\square$

**Answers to the mandated determinations.**
- *Minimum assumption:* joint exchangeability (DE-T1) — or, below it, any declared transport class (DE-T2); "identical distributions" is a declared strengthening.
- *What can be weakened:* (IID) → exchangeable → Λ-ratio or TV-ball transport, with explicit degradation formulas (DE-T2).
- *What failure occurs:* adversarial reversal — history worthless, committal rules strictly harmful (DE-T3); decision learning fails while identification stands untouched.
- *Robust ambiguity possible?* Yes — well-posed for convex compact classes, continuously interpolating between Bayes and the frozen minimax fallback (DE-T4); the fallback is itself the unique assumption-free position of the dial.
