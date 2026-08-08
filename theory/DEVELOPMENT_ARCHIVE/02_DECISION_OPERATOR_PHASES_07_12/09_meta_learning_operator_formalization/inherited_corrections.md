# Inherited Corrections (Selector Converse; Ranking Robustness)

> **Status:** Phase-9, 2026-08-03. The `../11_final_closure_audit/` refuted two further 8.2 statements outside the mandated Parts I–VI; since the Phase-9 chain cites the selection and ranking layers, they are corrected here so the final interface rests on no refuted claim. New results carry **ML-X** numbers, tagged **[proved]**. DC-S1 (the confinement jump theorem) and DC-R1–R4 were audited **valid** and are untouched.

---

## 1. ML-X1: the bridge converse of DC-S4(iii) is retracted

**The defect.** DC-S4(iii) claimed continuous selectors exist "whenever any argmin bridge connects the branches." Audit witness, adopted: $\mathcal A=[0,1]$, $\rho_t(a)=(2t-1)a$. Argmin: $\{1\}$ for $t<\tfrac12$, all of $[0,1]$ at $t=\tfrac12$, $\{0\}$ for $t>\tfrac12$. A bridge exists (at the single instant $t=\tfrac12$), yet every selector equals $1$ on $[0,\tfrac12)$ and $0$ on $(\tfrac12,1]$ — discontinuous regardless of its value at $\tfrac12$. An instantaneous bridge is not a path.

**Corrected statement (ML-X1). [proved]** The proved sufficiency results are exactly DC-S4(i)–(ii): (Berge) strict quasiconvexity per $t$ ⇒ unique continuous selector; (Michael) argmin correspondence l.s.c. with nonempty closed convex values ⇒ a continuous selector exists. Beyond these, existence of a continuous selector is equivalent — tautologically — to the argmin correspondence admitting a continuous selection, and **no bridge-existence shortcut is valid**: the witness shows the bridge must support a continuous transit over an interval, which is what l.s.c. encodes and what an instantaneous fat argmin fails to provide (the correspondence above is not l.s.c. at $\tfrac12$ from either side's perspective: nearby argmins are the *endpoints* $\{1\}$ and $\{0\}$, not nearby interior points). An exact characterization between the two proved sufficiency regimes is left open and is not consumed anywhere in the chain — the contract only ever *uses* (i), (ii), and the negative results. $\square$

## 2. ML-X2: DC-S5's "automatic confinement" rephrased

**The defect.** DC-S5 justified the discrete-action discontinuity via "(H2) holds automatically for any two actions" — false as stated when $\mathcal A$ has further actions the argmin can visit (confinement to a chosen *pair* is not automatic).

**Corrected statement (ML-X2). [proved]** For finite $\mathcal A$ with the discrete metric, no confinement framing is needed: **any map from a connected parameter interval to a discrete space that takes two distinct values is discontinuous** (its fibers would partition the interval into $\ge2$ disjoint nonempty relatively-closed... directly: a continuous map from a connected space to a discrete space is constant). Hence every selector whose output ordering changes along a continuous problem path jumps — the ranking realizability warning stands, now by the correct one-line argument, with no pairwise-confinement claim. $\square$

## 3. ML-X3: DC-R5's robustness "iff" corrected

**The defect.** DC-R5 characterized Tier-2 listwise robustness by interval separation ("one action's upper risk end below all others' lower ends"). Audit witness, adopted: risks $r_0(P_p)=1-p$, $r_1(P_p)=1.1-p$ over the class indexed by $p\in[0,1]$: action $0$ is strictly optimal under **every** law in the class, yet its risk interval $[0,1]$ is not separated below $[0.1,1.1]$. Separation is sufficient, not necessary — risks co-vary across the class.

**Corrected statement (ML-X3). [proved]**
*Definition (robustness, per-law):* action $\sigma^*$ is **decision-robust over $\widehat{\mathcal Q}$** iff $\sigma^*\in\arg\min_\sigma\mathbb E_P\,\ell(\sigma,\cdot)$ for every $P\in\widehat{\mathcal Q}$; equivalently
$$\sup_{P\in\widehat{\mathcal Q}}\ \Big[\mathbb E_P\,\ell(\sigma^*,\cdot)\ -\ \min_{\sigma}\mathbb E_P\,\ell(\sigma,\cdot)\Big]\ =\ 0.$$
*Decidability:* $\mathbb E_P\ell(\sigma^*)-\min_\sigma\mathbb E_P\ell(\sigma)=\max_\sigma\big[\mathbb E_P\ell(\sigma^*)-\mathbb E_P\ell(\sigma)\big]$ is a maximum of affine functions of $P$, hence convex in $P$; its supremum over the polytope $\widehat{\mathcal Q}$ is attained at a vertex, so robustness is finitely decidable (vertex check), and equivalently by LPs: $\sigma^*$ is robust iff $\max_{P\in\widehat{\mathcal Q}}\big[\mathbb E_P\ell(\sigma^*)-\mathbb E_P\ell(\sigma)\big]\le0$ for every competing $\sigma$ — each a linear program.
*Relation to the old test:* interval separation $\Rightarrow$ robustness (immediate), never conversely (witness above). DC-R5's confidence brackets and $\Gamma$-minimax LPs are unaffected; only the robustness *characterization* is replaced by the per-law/LP criterion, and DE-R6's Tier 2 is henceforth read with this definition. $\square$
