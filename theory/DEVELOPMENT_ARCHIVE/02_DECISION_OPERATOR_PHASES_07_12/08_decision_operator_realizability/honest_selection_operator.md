# Honest Selection Operator (Part II)

> **Status:** Phase-8, 2026-08-03. Phases 0–7 frozen and cited. New results carry **DR-S** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. Problem: the criterion argmin is in general **set-valued**; a trainable system must output something definite. This file derives the minimal valid object and characterizes exactly when a single-valued selector exists — with **no hidden measures** (the Phase-7 tie-break device DE-O2 used an auxiliary full-support $\mu_0$; Phase 8 surfaces that as a *declared* object or removes it).

---

## 1. The set-valued truth

For criterion risk $\rho(a)=\rho(a;J,\Delta,L)$ (Bayes, $\Gamma$-minimax over declared $\mathcal Q$, or the floor criterion $\sup_{v\in J}L(a,v)$):

$$\mathcal A^*(J,\Delta,L)\;=\;\operatorname*{arg\,min}_{a\in\mathcal A}\ \rho(a).$$

**Proposition DR-S1 (well-posedness). [proved]** With $\mathcal A$ compact and $\rho$ l.s.c. (DE-T4(i)), $\mathcal A^*$ is nonempty and closed. It is the **maximal honest answer**: every element is criterion-optimal, no element is privileged by anything so far declared. $\square$

---

## 2. Option A — the action set is the canonical object

**Definition.** Option A returns $\mathcal A^*(J,\Delta,L)$ itself (a set-valued operator). It consumes nothing beyond $(J,\Delta,L)$, hides nothing, and — by DR-S3 below — is the **only** canonical possibility when ties carry a symmetry.

## 3. Option B — declared tie-break $\tau$

**Definition.** A **tie-break** is a declared total preorder $\tau$ on $\mathcal A$ (or a declared selection function), fixed *before* seeing the current task's data, expressible in terms of declared structure only (loss, action-space structure, declared coordinates or reference points), and echoed in the output per H6. Then
$$D(J,\Delta,L,\tau)\;=\;\tau\text{-minimal element(s) of }\mathcal A^*(J,\Delta,L),$$
single-valued whenever $\tau$ is a total order with unique minima on closed sets (e.g. lexicographic on declared coordinates; minimal declared norm).
**Admissibility rule (no hidden measures):** $\tau$ is admissible iff it references no unidentified structure and no undeclared measure. The Phase-7 device "$\mu_0$-minimal element of the argmin" (DE-O2(ii)) is hereby reclassified: it is a **valid Option-B tie-break iff $\mu_0$ is declared in $\Delta$**; as an implicit implementation choice it is forbidden. This is a Phase-8 tightening of the interface, not a change to any Phase-7 theorem (DE-O2's existence claim is unaffected; only the selection step's bookkeeping moves into $\tau$).

---

## 4. When does a single-valued selector exist?

**Theorem DR-S2 (uniqueness from strict convexity — canonical selector exists). [proved]**
If $\mathcal A\subseteq\mathbb R^d$ is convex compact and $\rho$ is strictly quasiconvex, $\mathcal A^*$ is a singleton and Option A is already single-valued. Concretely: point prediction under squared loss with any criterion of the form $\rho(a)=\sup_{\mu\in\mathcal Q}\int\|a-v\|^2d\mu$ — each integrand is $\|a\|^2-2\langle a,m_\mu\rangle+c_\mu$; the common quadratic part makes the supremum $=\|a\|^2+\sup_\mu(\text{affine})$, strictly convex. **Squared-loss prediction never needs a tie-break.** Absolute loss (median-type intervals) and all discrete-action problems do not qualify. $\square$

**Theorem DR-S3 (equivariance obstruction — when no canonical selector can exist). [proved]**
Let a group $G$ act on the problem (on $\mathcal A$ and $\mathbb R^m$) leaving $J$, $L$, and $\Delta$ invariant. Then $\mathcal A^*$ is $G$-invariant as a set; a $G$-equivariant single-valued selector exists **iff** $\mathcal A^*$ contains a $G$-fixed point (equivalently, iff the action on $\mathcal A^*$ has a fixed point; for free actions on non-singleton $\mathcal A^*$: never). Witness: sign-symmetric pairwise ranking (DE-R3(iv)) — $G=\mathbb Z_2$ swaps the two orderings, fixes the problem, acts freely on the two-element $\mathcal A^*$: **no equivariant single-valued selector exists.** Any single-valued output therefore breaks the symmetry, and by the Phase-7 canonicity theorems (DE-S4/S5) the broken symmetry is exactly an undeclared asymmetric input — a hidden measure. Conclusion: **the dichotomy Option A / Option B is exhaustive and forced**: return the set, or declare $\tau$. A set-theoretic selector always exists (choice; measurable selection under standard conditions, Kuratowski–Ryll-Nardzewski), so the obstruction is canonicity, not existence — which is precisely why hiding it in an implementation is fabrication rather than necessity. $\square$

**Theorem DR-S4 (discontinuity at tie boundaries — the realizability warning). [proved]**
Along any continuous path of problems crossing a tie (e.g. an identified difference set $\Delta_{ab}$ sliding from $\subseteq(0,\infty)$ through the symmetric position), every single-valued selector has a jump discontinuity at the crossing (its value moves between separated components of $\mathcal A$, e.g. the two orderings). Hence any **continuous** approximator of a single-valued selector carries irreducible localized error near tie boundaries — the decision-layer recurrence of frozen MP-4 (optimal-operator discontinuity at section-topology transitions). $\square$

---

## 5. The realizable object: $\eta$-argmin with declared tolerance

**Definition.** $\mathcal A^*_\eta(J,\Delta,L)=\{a\in\mathcal A:\rho(a)\le\inf\rho+\eta\}$ for declared tolerance $\eta\ge0$.

**Theorem DR-S5 (honest selection for approximators). [proved]**
Let a system hold an approximate risk $\hat\rho$ with $\sup_{a}|\hat\rho(a)-\rho(a)|\le\eta/2$ (composed from outer-$\widehat J$ and class-estimation errors per DR-F4 / DE-U7). Then:
(i) every minimizer of $\hat\rho$ lies in $\mathcal A^*_\eta(\rho)$ — the system can validly claim "**$\eta$-optimal**", never "optimal";
(ii) $\mathcal A^*_\eta$ is outer-semicontinuous in $(\rho,\eta)$ and contains $\mathcal A^*$: reporting $\mathcal A^*_\eta$ (Option A at tolerance $\eta$) is conservative — the honest set-valued output degrades gracefully where the exact selector jumps (DR-S4);
(iii) with a declared $\tau$, the composed output "$\tau$-minimal element of the $\hat\rho$-argmin, claimed $\eta$-optimal" satisfies all Phase-7 honesty axioms; the claim's tolerance $\eta$ joins the ledger.
*Proof.* (i) $\rho(\hat a)\le\hat\rho(\hat a)+\eta/2\le\hat\rho(a^*)+\eta/2\le\rho(a^*)+\eta$. (ii) definitional. (iii) H1–H6 are untouched by the selection step (DE-O3); the new emitted quantity ($\eta$) is declared. $\square$

---

## 6. Summary — the honest selection theorem

$$\boxed{\begin{array}{c}\text{The minimal valid selection object is the set }\mathcal A^*\ (\text{realizably: }\mathcal A^*_\eta\text{ with declared }\eta).\ \text{A single-valued selector exists canonically iff}\\ \mathcal A^*\text{ is a singleton (guaranteed under strict quasiconvexity, e.g. squared loss); in the presence of a free problem symmetry it cannot exist}\\ \text{canonically (DR-S3), and the only honest alternative is a declared tie-break }\tau\text{ — never a hidden measure. Exact selectors are}\\ \text{discontinuous at ties (DR-S4); approximators must output }\eta\text{-argmin sets or accept flagged localized error there.}\end{array}}$$
