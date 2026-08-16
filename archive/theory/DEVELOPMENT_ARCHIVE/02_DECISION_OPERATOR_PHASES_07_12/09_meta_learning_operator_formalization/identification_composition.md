# Identification Composition (Part III)

> **Status:** Phase-9, 2026-08-03. Phases 0–7 frozen and cited. New results carry **ML-C** numbers, tagged **[proved] / [conditional]**. Connects the meta-learning operator to the frozen identification layer for the current task.

---

## 1. Current-task adaptation: the two arms

Current task $T_*=(O_*,S_*,Q_*,\gamma_*)$ (all components present; $S_*\subseteq O_*$).

**Identification arm (frozen + Phase-8 typing).**
$$I_\theta:\ \mathrm{set}(H_N)\times\mathcal O\times\mathcal X^m\ \longrightarrow\ \mathcal K_m\times\mathcal K^{\mathrm{fin}}_m\times\mathrm{Flags},\qquad I_\theta(O_*)\ :=\ I_\theta(\mathrm{set}(H_N),O_*,Q_*)=(\widehat J,\ \widetilde J,\ \mathrm{flags}),$$
with $V_I$: outer $\widehat J\supseteq J_{Q_*}(O_*)$ (order projection outer), certified witnesses $\widetilde J\subseteq J_{Q_*}(O_*)$, set-quotient input (ML-T4 — the 8.2 invariance axiom, now derived from the channel typing).

**Population arm (evaluation of the meta-object).**
$$M_\phi(H_N,S_*)\ :=\ \big[A_\phi(H_N)\big]\big(\kappa(O_*),\ \gamma_*\big)\ =\ \big(\widehat{\mathcal Q}_{*},\ 1-\delta,\ r\big),$$
i.e. the transferable object produced from the **multiset** of history, evaluated at the current context. Two typing facts make this well-posed and honest:
(i) the current task enters the population arm **only** through the declared context statistic $\kappa(O_*)$ — computable from the observable record ($S_*$ is part of it); any finer dependence would smuggle current-task information into the frequency channel (DC-C/ML-K);
(ii) the mandated signature $M_\phi(H_N,S_*)$ is exactly this evaluation: meta-learning ($A_\phi$, done once over history) and adaptation-time conditioning (evaluation at $\kappa(O_*)$) are one factored operator, so nothing about the current member beyond its context and — separately, in $D_\psi$ — its identified support ever touches population information.

---

## 2. The composition theorem

**Theorem ML-C1 (the composition produces the decision-relevant object). [proved, given the cited components]**
Under $V_I$, $V_A$, and the conditioning stack of `conditional_population_repair.md`:
$$\big(I_\theta(O_*),\ M_\phi(H_N,S_*)\big)\ =\ \Big(\underbrace{(\widehat J,\widetilde J,\mathrm{flags})}_{\text{identified joint object, outer + witnesses}},\ \underbrace{(\widehat{\mathcal Q}_*,1-\delta,r)}_{\text{rung-tagged population class}}\Big)$$
is precisely the typed domain of the decision operator $D_\psi$ (8.2 `final_interface.md`, with the Phase-9 corrections), and the pair is **decision-relevant** in the exact sense of the program: (a) the first component determines the dominance order, the floors/guarantee brackets, and the admissible-order object (DR-F4-R, DR-J); (b) the second determines the declared completion of the dominance order — the Phase-7 minimal decision primitive — at the rung its assumptions support (DC-C3 corrected in ML-K); (c) support restriction ($P(\widehat\Sigma)=1$, $P$ supported in $g(\widehat J)$) is applied inside $D_\psi$, is likelihood-free (DE-H2), and is the only point where the two arms meet. Consequently
$$D_\psi\big(\widehat J,\ \widehat{\mathcal Q}_*\!\upharpoonright_{\mathrm{support}},\ L,\ \tau\big)$$
is well-typed, and every emission is valid under DC-I1's (corrected) composition argument. $\square$

**Proposition ML-C2 (degeneracy behavior — the composition is total). [proved]**
Every partial or failed input yields a defined, honest output rather than an error: empty fiber ($N_{\kappa(O_*)}=0$) → rung-1 vacuous class → $D_\psi$ at the frozen minimax endpoint (DE-T4); off-coverage query → loss-typed guarantee values (8.2); empty section / empty conditioned class / untypable demand → the DC-A4 failure reports. The chain therefore never *requires* population information — it uses it exactly when the declarations support it, which is stop-condition 4's content. $\square$

**Remark (separation preserved).** The composition changes nothing upstream: $I(O_*)$, $J_{Q_*}(O_*)$, and all frozen certificates are $\phi$-independent and $H_N$-multiplicity-independent (H1 + DE-H2/H3 + ML-T4). Meta-learning tilts selections; identification bounds them. The two arms are composed, never merged.
