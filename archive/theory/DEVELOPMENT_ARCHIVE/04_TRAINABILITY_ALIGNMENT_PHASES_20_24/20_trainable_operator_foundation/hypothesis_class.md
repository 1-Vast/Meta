# Hypothesis Class (Deliverables 2–3)

> **Status:** Phase-20, 2026-08-03. Defines $\mathcal H=\{F_\omega\mid\omega\in\Omega\}$ and proves its admissible properties. New results **TF-4–TF-6**, tagged **[proved] / [declared]**.

---

## 1. The hypothesis class

**Definition TF-4 (operator hypothesis class).**
For $\omega=(\xi,b_1,\dots,b_m)\in\Omega$ and inputs $(z_H,S,Q,\gamma)$ (deployment state, support, query, specification — Phase-19.5 typing), define the **trainable operator**
$$F_\omega(S,Q,\gamma;z_H)\ =\ \Big(\ K\big(\ \mathsf{asm}(\,F_\xi(z),\,(b_j),\,b^{\mathrm{pop}}_{\kappa(S)}\,)\ \big)\big|_{\mathrm{supp}\,I(S)},\ \ \mathrm{confidence}(z_H),\ \ \mathrm{rung}(z_H)\ \Big)\ +\ \text{certificate row from }I(S),$$
where $z=z(S,Q,\gamma)$ is the statistic, $F_\xi=\pi_C\circ G(\xi,\cdot)$ (TF-3), $\mathsf{asm}$ the fixed convex assembly, $b^{\mathrm{pop}}$ from $z_H$ (frozen population band), and the side channels / certificate row canonical and $\omega$-invariant. The hypothesis class is
$$\mathcal H\ =\ \{\,F_\omega\ :\ \omega\in\Omega\,\},$$
indexed by the compact finite-dimensional $\Omega$. This is the object the audit found missing; it is now a set of fully-specified maps, each of type (input) $\to$ (valid operator value $+$ certificate).

## 2. Admissible properties — all five, proved

**Theorem TF-5 (admissibility of $\mathcal H$). [proved]**
Every $F_\omega\in\mathcal H$ satisfies:
1. **Measurability.** $(\omega,\text{input})\mapsto F_\omega(\text{input})$ is jointly measurable into $(\mathbb M,\text{evaluation σ-algebra})$: $z(\cdot)$ and $b^{\mathrm{pop}},I(S)$ are measurable in the input; $F_\xi$ is jointly continuous (TF-3); $\mathsf{asm}$ and $K(\cdot)$ are continuous into the evaluation topology; side channels/certificates are measurable functions of $z_H,I(S)$. Composition of (jointly) measurable maps.
2. **Continuity.** For fixed input, $\omega\mapsto F_\omega$ is continuous in the operator metric $d_{\mathbb M}$: $F_\xi$ is $L_\xi$-Lipschitz in $\xi$ (TF-3), $\mathsf{asm}$ is Lipschitz in $(\text{coeff},b_j)$, and the band-to-class map $K$ is Lipschitz in $d_{\mathbb M}$ by the stability theorems (Hoffman $\tfrac12\bar H$ on Route A; $D_V$ with mesh floor on Route B). Explicit constant in TF-6.
3. **Permutation invariance of support.** $F_\omega(\pi S,\dots)=F_\omega(S,\dots)$ for every permutation $\pi$ of the observation list: the input enters only through $z$, $b^{\mathrm{pop}}_{\kappa(S)}$, and $I(S)$, each a symmetric function of the support set (MI-7); the query/specification indices are held fixed (not permuted). Duplicate idempotence likewise (MI-2).
4. **Compatibility with set-valued output.** Every output is a valid $\mathbb M$-element — a nonempty compact convex (Route B: $W_1$-closed) constraint class $K(\cdot)$, coherent by construction, for **every** $\omega\in\Omega$ (MR-9 via $\mathsf{asm}$ landing in $\mathbb B$): validity is a property of the type, not of training.
5. **Operator metric consistency.** The class is a subset of the metric space $(\mathbb M,d_{\mathbb M})$ (Phase-11/PM-1 weighted metric); the map $\omega\mapsto F_\omega(\text{input})$ is Lipschitz into it (property 2), so $\mathcal H$ is a Lipschitz image of the compact $\Omega$ — compact in $d_{\mathbb M}$ per fixed input. $\square$

**Proposition TF-6 (the class Lipschitz constant). [proved]**
There is an explicit $L_{\mathcal H}$ with, for every input,
$$d_{\mathbb M}\big(F_\omega(\text{input}),\,F_{\omega'}(\text{input})\big)\ \le\ L_{\mathcal H}\,\|\omega-\omega'\|,\qquad L_{\mathcal H}=C_{\mathrm{stab}}\big(L_\xi+1\big),$$
$C_{\mathrm{stab}}=\tfrac12\bar H$ (Route A) or $D_V$ (Route B): the coefficient part contributes $C_{\mathrm{stab}}L_\xi\|\xi-\xi'\|$ (TF-3 + assembly + stability), the anchor part $C_{\mathrm{stab}}\|(b_j)-(b_j')\|$ (assembly weights sum to $\le1$). This constant is what the covering-number generalization bound (`generalization_bridge.md`) consumes. $\square$
