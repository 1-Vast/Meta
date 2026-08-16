# Current-Task Adaptation (§3)

> **Status:** Phase-9 closure, 2026-08-03. New results carry **MC-** numbers, tagged **[proved] / [conditional]**.

---

## 1. The two arms

New task $T_*=(O_*,S_*,Q_*,\gamma_*)$, all components present.

**Identification arm (frozen theory, Phase-8 typing).**
$$I_\theta:\ \mathrm{set}(H_N)\times\mathcal O\times\mathcal X^m\ \longrightarrow\ \mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathrm{Flags},\qquad I_\theta(O_*)=(\widehat J,\ \widetilde J,\ \mathrm{flags}),$$
valid iff $\widetilde J\subseteq J_{Q_*}(O_*)\subseteq\widehat J$ (witnesses inner, envelope outer; order projection $\widehat\Sigma\supseteq\Sigma$) under the declared closure class. Input is the **set** quotient (MC-4).

**Population arm (evaluation of the meta-object).**
$$A_\phi(H_N)\in\mathbb M,\qquad \Delta_{\mathrm{pop}}\ :=\ \big[A_\phi(H_N)\big]\big(\kappa(O_*),\ \gamma_*\big)\ =\ \big(\widehat{\mathcal Q}_*,\ 1-\delta,\ r\big).$$
The current task touches population information **only** through the declared context statistic $\kappa(O_*)$ (a function of the observable record; $S_*\subseteq O_*$) — any finer dependence would smuggle current-task specifics into the frequency channel (§4). Meta-learning (one pass over $\mathrm{mult}(H_N)$) and adaptation-time conditioning (evaluation at the current context) are one factored operator.

## 2. Composition produces decision-sufficient information

**Theorem MC-9 (decision sufficiency of the composed pair). [proved, given the cited components]**
Under $V_I$, $V_A$, and the §4 conditioning stack, the pair
$$\Big(\ (\widehat J,\widetilde J,\mathrm{flags}),\ \ (\widehat{\mathcal Q}_*,1-\delta,r)\ \Big)$$
is **decision-sufficient** for the declared specification $\gamma_*$, in the exact sense of the corpus:
(a) the identification component determines the dominance order, the admissible action set, the loss-typed floor bracket $[R_{\mathrm{set/rand}}(\widetilde J),\ \sup_{v\in\widehat J}L(\cdot,v)]$, and the admissible-order object — everything the current task's observations determine about the decision (DE-S2, DR-F4-R, DR-J);
(b) the population component supplies the Phase-7 minimal decision primitive — a declared, rung-tagged generator of a monotone completion of that dominance order — at exactly the rung its assumptions support (§4);
(c) the two meet only at the likelihood-free support restriction ($P(\widehat\Sigma)=1$, mass inside $g(\widehat J)$), applied inside the decision operator;
(d) nothing else is needed: $D_\psi$'s typed domain is exactly this pair plus the declared context, and every judgment $D_\psi$ must make is either a cited theorem or a declared slot (8.2 DC-I2 with the Phase-9 corrections).
Hence $D_\psi\big(\widehat J,\ \Delta_{\mathrm{pop}}\!\upharpoonright_{\mathrm{support}},\ L,\ \tau\big)$ is well-typed and every emission is valid at its stated type. $\square$

**Proposition MC-10 (totality — adaptation never errors, it degrades). [proved]**
Every partial input yields a defined honest output: $N_{\kappa(O_*)}=0$ → rung-1 vacuous class → frozen minimax endpoint (DE-T4); off-coverage query → loss-typed guarantee values; empty section / empty restricted class / untypable single-valued demand → the DC-A4 failure reports (never priced as actions). The composition *uses* population information exactly when the declarations earn it, and *never requires* it. $\square$

**Separation, restated once.** The composed operator changes nothing upstream: $I(O_*)$, $J_{Q_*}(O_*)$, all frozen certificates are $\phi$-independent and multiplicity-independent (H1, DE-H2/H3, MC-4). Meta-learning tilts selections inside the identified set; identification bounds them. Composed, never merged.
