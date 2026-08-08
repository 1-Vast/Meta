# Learnability of Decision Information (Part IV)

> **Status:** Phase-8, 2026-08-03. Phases 0–7 frozen and cited (DE-H, DE-L, DE-T; CI-A for fibers). New results carry **DR-L** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. No architectures. Question: what, exactly, may a system learn from historical tasks, under which declared assumptions, so that the composed decision operator remains honest?

---

## 1. The separation, restated as a typing rule

- **Current-task information:** $I(O)$, $J_Q(O)$, and their quotients ($\Delta_{ab}$, $\Sigma$) — produced by the frozen identification theory from the current member's observations only. Type: *unconditional*.
- **Population information:** $\Delta_{\mathrm{population}}$ — everything imported from historical tasks. Type: *conditional on declared cross-task axioms*, and confined by DE-H2/H3 to **weighting the inside** of the current identified object.

**Definition DR-L1 ($\Delta_{\mathrm{population}}$).** The population decision object is a triple
$$\Delta_{\mathrm{population}}=\big(\ \widehat{\mathcal Q}\ ,\ 1-\delta\ ,\ \text{axiom tags}\ \big),$$
where $\widehat{\mathcal Q}$ is a convex, weak-*-compact **ambiguity class of laws on the decision-relevant pushforward space** ($g$-space: $\mathbb R$ for a difference, $[0,1]$ for a sign frequency, $\Delta(S_m)$-parameters for listwise), $1-\delta$ its confidence, and the tags the exact axioms consumed. The single-law case is the singleton class (DE-H5's degenerate endpoint). Nothing else about the population is a legitimate learning target: learning *feasibility* (windows) is the frozen archive channel (DE-H1(i)); learning member-level structure is a declared closure class, not an estimate.

---

## 2. The assumption ladder (each rung declared, each priced)

**DR-L2 (the four rungs). [proved / conditional as marked]**
(i) **Exchangeability (EXCH/IID).** Joint exchangeability of $(f_1,\dots,f_n,f_\beta)$: the minimal identical-treatment axiom (DE-T1); buys the Phase-7 rates — Hoeffding/DKW at $n^{-1/2}$ with explicit constants (DE-L3) for pushforward functionals identified per task.
(ii) **Conditional exchangeability (C-EXCH). [new formalization]** Tasks carry an observed auxiliary label $c\in C$ (the Phase-6 datum); assume exchangeability only **within fibers**: $(f_i)_{c_i=c}$ exchangeable with $f_\beta$ given $c_\beta=c$. Then all rung-(i) guarantees hold with $n\to n_{c_\beta}=\#\{i:c_i=c_\beta\}$ — the effective sample is the fiber count. **Borrowing across fibers is impossible without declared structure on $C$:** any smoothing $c\mapsto\pi(\cdot\mid c)$ requires a declared modulus of continuity, and no such modulus is derivable from data — this is exactly CI-A5(iii) lifted to the population level. [conditional on (C-EXCH); the impossibility part proved]
(iii) **Transport (TRANS).** A declared class linking historical and current populations — $\Lambda$-ratio band or TV-ball $\rho$ — with the explicit degradation formulas and the robust-decision thresholds of DE-T2. Without any declared class, DE-T3's adversarial reversal applies: history is worthless and committal rules harmful. [frozen citations; nothing new needed]
(iv) **Ambiguity-set semantics (AMB).** Whatever is learned enters the decision operator only as the class $\widehat{\mathcal Q}$ with **outer validity**: $\Pr(\text{true conditioned law}\in\widehat{\mathcal Q})\ge1-\delta$ under the declared rungs. An under-covering class is the population-level analogue of an inner envelope — false conditional certificates (DR-F4(iii) one level up). [proved by the same argument]

---

## 3. The joint-object learnability condition

**Theorem DR-L3 (learnability of the decision information — the composed condition). [conditional on the declared rungs]**
Fix the decision-relevant map $g$ (Part III). $\Delta_{\mathrm{population}}$ for $g$ is learnable from $n$ historical tasks with a valid confidence class iff the following three gates hold, and the resulting class has the stated size:
- **(a) Per-task gate (frozen identification, applied to history):** for each historical task $i$, $g(f_i)$ is identified — or interval/set-identified — by task $i$'s **own** data under the frozen conditions (F17 coverage; Thm 7.1 for differences; the order-type analogue: task $i$'s own $\Sigma^{(i)}$). Tasks contribute what their designs identify, nothing more.
- **(b) Cross-task gate:** a declared rung (i) or (ii) with effective count $n_{\mathrm{eff}}$ ($=n$ or $n_{c_\beta}$).
- **(c) Shift gate:** a declared rung (iii) class of radius $\rho$ (possibly $\rho=0$ under (IID)).
Then, for event/frequency targets (pairwise sign; each listwise order $\sigma$), with probability $\ge1-\delta$:
$$p(\sigma)\ \in\ \Big[\ \underbrace{\tfrac1{n_{\mathrm{eff}}}\#\{i:\Sigma^{(i)}=\{\sigma\}\}}_{\text{forced}}-\eta_{n_{\mathrm{eff}}}-\rho\ ,\ \ \underbrace{\tfrac1{n_{\mathrm{eff}}}\#\{i:\sigma\in\Sigma^{(i)}\}}_{\text{compatible}}+\eta_{n_{\mathrm{eff}}}+\rho\ \Big],\qquad \eta_n=\sqrt{\ln(4/\delta)/2n},$$
— the listwise generalization of DE-R5/DE-L4 (pointwise bounds $\mathbf 1\{\Sigma^{(i)}=\{\sigma\}\}\le\mathbf 1\{\operatorname{ord}(f_i)=\sigma\}\le\mathbf 1\{\sigma\in\Sigma^{(i)}\}$, two Hoeffding applications, TV-transport shift). The class $\widehat{\mathcal Q}$ is the polytope of laws on $S_m$ consistent with these $m!$ (or the queried subset of) intervals — convex and compact, as required by DE-T4/DR-S1. $\square$

**Corollary DR-L4 (the three widths, and what reduces each). [proved]**
The interval width decomposes as (systematic censoring width: population fraction of order-undecided tasks — reduced only by richer historical designs, never by $n$) $+$ ($2\eta_{n_{\mathrm{eff}}}$: reduced by more tasks *in the fiber*) $+$ ($2\rho$: reduced only by tighter declared transport). This is DE-U8 specialized to the ranking application, with the (C-EXCH) refinement: task count outside the fiber reduces **nothing**. $\square$

**Impossibility restated for the record. [proved, frozen]** Delete gate (b): DE-T3 — no learning. Delete gate (a): historical tasks contribute vacuous intervals $[0,1]$ — the class is all laws, $\Gamma$-minimax collapses to the frozen floor (DE-T4(iii)): *graceful, honest degradation to minimax is the failure mode, by construction.*

---

## 4. Summary

$$\boxed{\begin{array}{c}\Delta_{\mathrm{population}}\ \text{is an outer-valid confidence class of laws on the decision-relevant pushforward, learnable iff}\\ \text{(per-task identification)}\ \wedge\ \text{(declared exchangeability, fiber-relative)}\ \wedge\ \text{(declared transport)};\ \text{its width has three separately-priced terms;}\\ \text{and every gate failure degrades the operator to a strictly more conservative, still-valid position of the DE-T4 dial — never to a false claim.}\end{array}}$$
