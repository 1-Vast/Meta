# Selector Theorem Repair (Closure Target 1)

> **Status:** Phase-8.2, 2026-08-03. Replaces DR-S4-R (`../08_decision_operator_realizability_repaired/honest_selection_operator.md` §2), which the audit correctly refuted: the union-of-branches hypothesis appeared only in the proof. Everything else in that file (DR-S1–S3, DR-S5) stands. New results: **DC-S1–S5**, tagged **[proved]**.

Setting: $\mathcal A$ a metric space; a continuous path of decision problems $\{P_t\}_{t\in[0,1]}$ with criterion risk $\rho:[0,1]\times\mathcal A\to\mathbb R$ jointly continuous; argmin correspondence $\mathcal A^*(t)=\arg\min_a\rho(t,a)$, assumed nonempty (compactness + l.s.c., DR-S1). A **selector** is any map $s:[0,1]\to\mathcal A$ with $s(t)\in\mathcal A^*(t)$ for all $t$.

---

## 1. The corrected discontinuity theorem

**Theorem DC-S1 (branch-switching discontinuity — all hypotheses in the statement). [proved]**
Assume:
(H1) $A_0,A_1\subseteq\mathcal A$ are closed with $d:=\operatorname{dist}(A_0,A_1)>0$;
(H2) **union confinement:** $\mathcal A^*(t)\subseteq A_0\cup A_1$ for **every** $t\in[0,1]$;
(H3) endpoint assignment: $\mathcal A^*(0)\subseteq A_0$ and $\mathcal A^*(1)\subseteq A_1$.
Then every selector $s$ is discontinuous, and at some $t^*$ its oscillation is $\ge d$.
*Proof.* Let $T_j=\{t:s(t)\in A_j\}$, $j=0,1$. By (H2) they cover $[0,1]$; by (H1) they are disjoint; by (H3) $0\in T_0$, $1\in T_1$. If $s$ were continuous, each $T_j=s^{-1}(A_j)$ would be closed, splitting $[0,1]$ into two disjoint nonempty closed sets — contradicting connectedness. So $s$ is discontinuous; moreover, taking $t^*=\sup T_0\cap[0,\inf T_1]$-type boundary point (any $t^*$ in $\overline{T_0}\cap\overline{T_1}$, nonempty since the two sets cover the connected $[0,1]$ and are not both closed-and-open), every neighborhood of $t^*$ contains points of $T_0$ and of $T_1$, so the oscillation of $s$ at $t^*$ is $\ge\operatorname{dist}(A_0,A_1)=d$. $\square$

Hypothesis (H2) is the previously missing condition; it is what "branch-switching" *means*: the argmin never leaves the two separated branches, so no bridge exists.

---

## 2. Why weaker assumptions fail — two counterexamples

**DC-S2 (endpoint separation alone is insufficient — the audit's witness, adopted). [proved]**
$\mathcal A=[0,1]$, $\rho_t(a)=(a-t)^2$, $A_0=\{0\}$, $A_1=\{1\}$. (H1) and (H3) hold ($\mathcal A^*(0)=\{0\}$, $\mathcal A^*(1)=\{1\}$); (H2) fails ($\mathcal A^*(t)=\{t\}\not\subseteq\{0,1\}$ for $t\in(0,1)$). The unique selector $s(t)=t$ is continuous. Endpoint minimizers in separated sets therefore force nothing: the argmin *travels* between the branches. $\square$

**DC-S3 (disconnected endpoint branches + connected $\mathcal A$ still insufficient — the fat-bridge witness). [proved]**
$\mathcal A=[0,1]$, $\rho_t(a)=\max(0,1-3t)\,a+\max(0,3t-2)\,(1-a)$ (jointly continuous). Then $\mathcal A^*(t)=\{0\}$ for $t<\tfrac13$, $=[0,1]$ for $t\in[\tfrac13,\tfrac23]$, $=\{1\}$ for $t>\tfrac23$. Endpoint argmins are the separated $\{0\},\{1\}$; the feasible action space is connected; and the continuous selector $s(t)=0$ on $[0,\tfrac13]$, $s(t)=3t-1$ on $[\tfrac13,\tfrac23]$, $s(t)=1$ on $[\tfrac23,1]$ stays in $\mathcal A^*(t)$ throughout. A temporarily *fat connected* argmin is a bridge just as a traveling singleton is. $\square$

Together: neither endpoint separation, nor connectedness of $\mathcal A$, nor even momentary disconnection of $\mathcal A^*(t)$ forces a jump. Exactly (H2) does.

---

## 3. The positive side — when continuous selectors exist

**Proposition DC-S4 (existence of continuous selectors). [proved / classical]**
(i) *Berge route:* if $\rho(t,\cdot)$ is strictly quasiconvex on convex compact $\mathcal A\subseteq\mathbb R^d$ for each $t$, then $\mathcal A^*(t)$ is a singleton and $t\mapsto\mathcal A^*(t)$ is continuous (Berge maximum theorem: u.s.c. + single-valued): the unique selector is continuous. (Squared-loss prediction, DR-S2.)
(ii) *Michael route:* if $t\mapsto\mathcal A^*(t)$ is lower semicontinuous with nonempty closed convex values in a Banach space, a continuous selector exists (Michael's selection theorem). This is sufficient, not necessary — DC-S3's correspondence is not l.s.c. at $t=\tfrac13$ yet admits a continuous selector.
(iii) The trichotomy demanded by the mandate: **disconnected argmin branches** obstruct continuity only under (H2)-confinement (DC-S1); a **connected feasible action space** neither creates nor removes the obstruction (DC-S3); **continuous selectors exist** in the convex regimes (i)–(ii) and whenever any argmin bridge connects the branches (DC-S2/S3). $\square$

**Lemma DC-S5 (discrete actions: the ranking warning survives unconditionally). [proved]**
If $\mathcal A$ is finite (orderings, discrete structured actions) with the discrete metric, then (H2) holds automatically for any two actions ($\mathcal A^*(t)\subseteq\mathcal A$, every pair of distinct actions is separated), and DC-S1 specializes to: **any selector whose output differs at two parameter values is discontinuous** (a nonconstant map from a connected interval to a discrete space). Hence for ranking decisions the Phase-8 realizability warning stands with no extra hypothesis: whenever the optimal ordering changes along a continuous problem path, every single-valued selector jumps, and continuous approximators near the switch must output $\mathcal A^*_\eta$-sets or the declared abstention action (priced per `abstention_semantics_repair.md`). $\square$

---

## 4. Summary

$$\boxed{\begin{array}{c}\text{Jump theorem: separated branches }+\ \textbf{union confinement of every intermediate argmin}\ +\text{ endpoint assignment }\Rightarrow\text{ oscillation }\ge d.\\ \text{Both witnesses (traveling singleton; fat bridge) show the confinement hypothesis is indispensable. Continuous selectors: Berge (strict quasiconvexity),}\\ \text{Michael (l.s.c. + closed convex values). Finite action spaces satisfy confinement trivially — the ranking discontinuity warning is unconditional.}\end{array}}$$
