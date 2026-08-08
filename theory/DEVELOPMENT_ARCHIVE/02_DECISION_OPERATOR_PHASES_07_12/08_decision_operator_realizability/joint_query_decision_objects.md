# Joint Query Decision Objects (Part III)

> **Status:** Phase-8, 2026-08-03. Phases 0–7 frozen and cited (especially DE-J, DE-R). New results carry **DR-J** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. The future application requires ranking; the realizable carrier is therefore the joint object $J_Q(O)$ and its decision-relevant quotients — never scalar sections alone. This file fixes, for each decision type, the **minimal decision-sufficient object** a system must represent.

---

## 1. The quotient lattice of $J_Q(O)$

For $Q=\{x_1,\dots,x_m\}$ and $J=J_Q(O)\subseteq\mathbb R^m$, each decision type consumes a **pushforward** $g(J)$ under its decision-relevant map $g$:

| Decision | map $g$ | pushforward object |
|---|---|---|
| point prediction, sup-loss / separable Bayes | coordinates $v\mapsto v_j$ | marginal sections / marginal laws |
| pairwise ranking $(x_a,x_b)$ | $v\mapsto v_a-v_b$ | difference sets $\Delta_{ab}$ (DE-J6) |
| listwise ranking | $v\mapsto\operatorname{ord}(v)$ (order type in $S_m$, ties resolved by declared convention) | admissible-order set $\Sigma(J)=\operatorname{ord}(J)\subseteq S_m$ |
| general joint loss | $g=\mathrm{id}$ | $J$ itself (DE-J7) |

**Proposition DR-J1 (sufficiency per loss class). [proved]** If $L(a,v)=\ell(a,g(v))$, then $R_{\mathrm{set}}$, all criterion risks, and the argmin sets depend on $J$ only through $g(J)$ (and on the weighting only through its $g$-pushforward): the pushforward is decision-sufficient. *Proof:* every $\sup_{v\in J}$ and $\int\cdot\,d\mu$ factors through $g$. $\square$

**Theorem DR-J2 (strictness of the lattice — coarser objects are insufficient). [proved]**
None of the quotients is recoverable from the coarser ones:
(i) marginals $\nrightarrow$ differences: DE-J6's equal-marginal triple (diagonal / anti-diagonal / square), pairwise minimax error $0$ vs $1$ vs $1$;
(ii) differences $\nrightarrow$ order set: see DR-J3(ii) below — the pairwise object over-admits orders;
(iii) order set $\nrightarrow$ $J$: $\Sigma$ discards all magnitudes (any order-preserving deformation of $J$ has the same $\Sigma$ but different prediction risks). Each decision tier needs exactly its own quotient; a system serving several tiers must carry the finest one demanded. $\square$

---

## 2. The three mandated decisions

**Point prediction. [frozen + DR-F]** Sup-loss minimax: coordinatewise Chebyshev centers, marginal sections suffice (DE-J2); separable Bayes: marginal laws (DE-J3/J4); non-separable or joint-minimax: $J$ itself (DE-J4(ii)); floors from DR-F1(i)–(ii).

**Pairwise ranking. [frozen + DR-F]** The complete object is the family $\{\Delta_{ab}\}_{a<b}$ with the trichotomy DE-R1, the preference number $p_{ab}$ where declared (DE-R4), and floor $R_{\mathrm{set}}\in\{0,1\}$ per pair ($R_{\mathrm{rand}}=\tfrac12$ at ties).

**Listwise ranking — the new object. [proved]**

**Theorem DR-J3.**
(i) *Sufficiency:* for any listwise loss depending on $v$ through its order (Kendall discordance, Spearman footrule on ranks, top-$k$ membership), $\Sigma(J)$ — with the weighting's pushforward onto $S_m$ where declared — is decision-sufficient (DR-J1), and the listwise floor is
$$R_{\mathrm{set}}=\min_{\sigma\in\mathcal A}\ \max_{\pi\in\Sigma(J)}\ \ell(\sigma,\pi).$$
(ii) *The pairwise object strictly over-admits:* the sign data $\{\operatorname{sgn}\Delta_{ab}\}$ determine only the **pairwise-compatible** order set $\Sigma^{\mathrm{pair}}(J)=\{\sigma\in S_m:\ \forall a<b,\ \sigma\text{'s }(a,b)\text{-comparison is realized by some member}\}\ \supseteq\ \Sigma(J)$, and the inclusion is strict in general because pairwise signs are realized by *different* members while an order must be realized by *one*. Witness ($m=3$): $J=\{(0,1,2),\,(2,1,0)\}$ — every pairwise difference set contains both signs, so $\Sigma^{\mathrm{pair}}=S_3$ (all $6$ orders), but $\Sigma(J)=\{321\text{-order},\,123\text{-order}\}$: exactly the two full reversals.
(iii) *The over-admission costs real floor:* Kendall discordance on the witness — against $\Sigma(J)$ (two reversals, Kendall distance $3$ apart), the triangle inequality forces $\max\ge\lceil 3/2\rceil=2$, attained (e.g. the order $2\!\succ\!1\!\succ\!3$: discordance $2$ vs one truth, $1$ vs the other): $R_{\mathrm{set}}=2$. Against $\Sigma^{\mathrm{pair}}=S_3$, the adversary plays the reversal of any $\sigma$: $R_{\mathrm{set}}^{\mathrm{pair}}=3$. A system carrying only pairwise information reports a floor of $3$ where the true floor is $2$ — conservative here, but the dual failure (Σ under-approximated, floor falsely low) is forbidden by DR-F4(iii); and the **selected actions differ** ($\Sigma^{\mathrm{pair}}$ makes all orders minimax-equivalent; $\Sigma(J)$ singles out the middle orders). $\square$

**Corollary DR-J4 (minimal decision-sufficient object for the ranking application). [proved]**
For an application whose decisions are pairwise and listwise rankings over declared query sets $Q$: the minimal object is the **admissible-order set $\Sigma(J_Q(O))$ per queried $Q$** (which determines every $\operatorname{sgn}\Delta_{ab}$ by projection onto transpositions, but not conversely), together with — exactly when preferences beyond minimax/abstention are demanded — the declared population weighting's pushforward onto $S_m$ (the $p(\sigma)$ data of DR-L3). Magnitude information ($J$ beyond $\Sigma$) is required only if graded losses (margin-weighted discordance) are declared. $\square$

---

## 3. Realizable representation constraints

**Outer semantics for orders. [proved]** The system's representation $\widehat\Sigma$ must satisfy $\widehat\Sigma\supseteq\Sigma(J_Q(O))$ under the declared closure class (an order may be spuriously admitted, never spuriously excluded): excluded-but-admissible orders are false ranking certificates (Tier-1 claims without identification — the DE-R6 tiers collapse). $\Sigma^{\mathrm{pair}}$ is one always-available outer approximation (DR-J3(ii)); tightening it toward $\Sigma$ is a *tightness* objective, never a validity requirement. Conversely every Tier-1 ranking claim emitted must be certified by $\widehat\Sigma$'s **projection being a singleton on the relevant comparison** — checkable, since outer $\widehat\Sigma$ singleton $\Rightarrow$ true singleton.

**Size discipline. [proved]** $|\Sigma(J)|\le|S_m|=m!$ but the object is generated by $J$'s intersection pattern with the $\binom m2$ hyperplanes $\{v_a=v_b\}$; for convex $J$ (linear regime with interval noise), $\Sigma(J)$ is the set of orders of points of a convex set — representable by its $\binom m2$ sign-interval data **plus** the joint realizability predicate; the witness of DR-J3(ii) shows the predicate is not redundant. A representation carrying only pairwise marginals silently substitutes $\Sigma^{\mathrm{pair}}$ — permitted as an outer proxy, flagged as non-tight. $\square$
