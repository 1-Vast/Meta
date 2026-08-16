# Preliminaries: All Base Objects (Item 3 + shared symbols)

> **Status:** Phase-25.1 (final self-contained theory closure), 2026-08-03. Phases 0–24 are the provenance; this package is **self-contained** — every symbol used by any retained theorem is defined here or in the two companion files, with no out-of-folder reference required to audit. No theory is redesigned, no claim strengthened, no continuum/ranking/varying-$z_H$ content restored. Definitions carry **FC-** numbers. Where a definition merely fixes the audited notation issue (Item 3), it is marked **[notation]**; all others are **[definition]** transcribed from the passed earlier results.

---

## 1. The fixed deployment (frozen throughout)

**FC-1 [definition].** The deployment is the fixed tuple $\mathcal D=(z_H^0,B(\cdot),\Delta_m,\mu,h)$: frozen state $z_H^0$; band rule $B(\cdot)$ (FC-4); coefficient simplex $\Delta_m$ (FC-3); ridge modulus $\mu>0$; output value mesh $h>0$. All quantities below are relative to this one $\mathcal D$; no coordinate of $\mathcal D$ ever varies.

## 2. Spaces

**FC-2 [definition].**
- **Statistic domain $Z$:** a compact metric space (finite union of compact cubes under the skeleton); $d_Z$ its metric. The statistic is $z=z(S,Q,\gamma)\in Z$ for a support $S$, query $Q$, specification $\gamma$.
- **Context map $\kappa:Z\to C_\kappa$**, $C_\kappa$ **finite** (declared); measurable.
- **Value space $V=[a_{\min},a_{\max}]\subset\mathbb R$**, compact; $D_V^{\mathrm{val}}:=a_{\max}-a_{\min}$ (the value-space diameter).
- **Law space $(\Delta(V),W_1)$:** Borel probability measures on $V$ with the Wasserstein-1 metric $W_1(P,P')=\int_V|F_P-F_{P'}|\,dv$; compact.

**FC-3 [definition] Coefficient space.** $C=\Delta_m=\{p\in\mathbb R^{m+1}:p_k\ge0,\ \sum_kp_k=1\}$, the $m$-simplex; norm $\|\cdot\|$ = Euclidean on $\mathbb R^{m+1}$; compact convex.

## 3. The band rule $B(z)$ (Item 3 — notation fixed)

**FC-4 [notation].** The assembly is a **fixed, deployment-determined matrix-valued rule**
$$B(z)=\big[\ \beta_0(z)\ \big|\ \beta_1\ \big|\ \cdots\ \big|\ \beta_m\ \big]\in\mathbb R^{(\dim\mathbb B)\times(m+1)},\qquad \beta_0(z)=b^{\mathrm{pop}}_{\kappa(z)},\ \ \beta_1,\dots,\beta_m\ \text{fixed anchors},$$
so the audited inconsistency ("$B$ called a fixed matrix while $\beta_0=b^{\mathrm{pop}}_{\kappa(z)}$") is resolved: $B$ is **not** one constant matrix; it is the rule $z\mapsto B(z)$, constant on each of the finitely many context cells $\{\kappa=c\}$. The assembly map is **pointwise linear in $p$**: $p\mapsto B(z)p=\sum_{k=0}^mp_k\beta_k(z)\in\mathbb B$ (a valid band vector, since $\{\beta_k(z)\}\subset\mathbb B$ and $\mathbb B$ is convex, $p\in\Delta_m$). Define the **assembly-norm constant**
$$\kappa_B\ :=\ \sup_{z\in Z}\ \|B(z)\|_{\mathrm{op}}\ <\ \infty$$
(finite: $\kappa(z)$ takes finitely many values, each $B$ a fixed matrix with bounded columns). Used in $D_V$ (FC-9).

## 4. Band space, class map, operator values

**FC-5 [definition].**
- **Band space $\mathbb B$:** the valid-description polytope — CDF-band vectors on the fixed grid of $V$ at mesh $h$, satisfying $0\le l_j\le u_j\le1$ and band monotonicity, with the closed/open convention (lower bounds on closed intervals, upper on open) making the induced class $W_1$-closed. $\mathbb B$ is compact convex; norm $\|\cdot\|_{\mathbb B}$ = sup over endpoint coordinates.
- **Class map $K:\mathbb B\to$ (nonempty compact convex subsets of $(\Delta(V),W_1)$):** $K(\beta)=\{P\in\Delta(V):\ P([a_{\min},t_j])\ge l_j,\ P([a_{\min},t_j))\le u_j\ \forall j\}$; nonempty and $W_1$-closed (FC-5 convention).
- **Operator value space $\mathbb M$:** triples (probability object $\in\{K(\beta)\}$, confidence $\in[0,1]$, rung $\in\{1,2,3,4\}$), the probability object further restricted to $\mathrm{supp}\,I(S)$; plus the $\omega$-invariant certificate row. The retained theorems compare two operator values differing **only** in the probability object (confidence, rung, certificate are computed from $z_H^0,I(S)$ alone, identical for the two maps compared), so those coordinates cancel in every distance below.

## 5. Loss, conditional risk, operative risk, target

**FC-6 [definition].** **Loss** $L:\mathbb B\times V\to[0,\infty)$, the declared interval/band score: **convex and $L_{\mathrm{Lip}}$-Lipschitz in the band argument, bounded by $\bar L$**, uniformly in the target value (A-LOSS). **Base conditional risk** $L_0(z,\beta)=\mathbb E[L(\beta,A_T)\mid\zeta=z]$ (A-STAT gives the conditional expectation; A-CONT selects a continuous, everywhere-defined version). **Operative regularized risk**
$$J_\mu(z,p)\ =\ L_0\big(z,B(z)p\big)+\tfrac\mu2\|p\|^2,$$
convex in $p$ (convex $L_0(z,\cdot)$ ∘ linear $B(z)$, plus ridge), **$\ge\mu$-strongly convex** on $\Delta_m$. **Target**
$$g^\star_\mu(z)=\operatorname*{arg\,min}_{p\in\Delta_m}J_\mu(z,p)\ \ \text{— unique, everywhere-defined (strong convexity + compactness + A-CONT).}$$
$L_p:=$ the Lipschitz constant of $p\mapsto J_\mu(z,p)$ on $\Delta_m$, uniform in $z$ (finite: $J_\mu(z,\cdot)$ convex and finite on the compact $\Delta_m$; used in the consistency decomposition).

All symbols appearing in the retained calibration and consistency theorems that are *not* statistical/learning objects are now defined; the remaining ones ($R_\mu,\mathcal E_\mu,d_{\mathbb M},D_V,\Phi$ and the meta-learning contract) are defined in the two companion files.
