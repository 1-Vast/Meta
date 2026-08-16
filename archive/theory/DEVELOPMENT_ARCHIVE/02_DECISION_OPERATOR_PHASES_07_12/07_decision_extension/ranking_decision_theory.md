# Ranking Decision Theory (Part VI)

> **Status:** Phase-7, 2026-08-03. Frozen corpus cited, not modified. New results carry **DE-R** numbers, tagged **[proved] / [conditional] / [impossible] / [open]**. The action is an **ordering** of two queries $x_a,x_b$: $\mathcal A\subseteq\{a\!\succ\!b,\ b\!\succ\!a,\ \text{abstain}\}$. The state functional is $\delta=f(x_a)-f(x_b)$ with identified set
> $$\Delta_{ab}(O)=\{f(x_a)-f(x_b):f\in I(O)\}.$$

---

## 1. Identifiability of the sign — the trichotomy

**Definition.** The **sign** is $\operatorname{sgn}\delta\in\{+,0,-\}$; the ordering decision is correct iff it matches $\operatorname{sgn}\delta$ (ties handled explicitly below).

**Proposition DE-R1 (trichotomy). [proved]**
Exactly one of:
(i) $\Delta_{ab}(O)\subseteq(0,\infty)$: the ordering $a\succ b$ holds for **every** admissible member — the decision is *identified*; no decision primitive is consumed; the frozen certificate applies (correctness is worst-case guaranteed).
(ii) $\Delta_{ab}(O)\subseteq(-\infty,0)$: symmetrically, $b\succ a$ identified.
(iii) $\Delta_{ab}(O)$ meets both $[0,\infty)$ and $(-\infty,0]$ nontrivially (or $=\{0\}$, the identified tie): the sign is **not identified**; any strict ordering output is a decision under residual ambiguity. $\square$

When is case (i)/(ii) checkable? In the frozen linear regime, $\delta$ is a linear functional; frozen Thm 7.1 gives the exact test ($\phi(x_a)-\phi(x_b)\in\operatorname{row}(G)$ identifies $\delta$ itself; with noise, $\Delta_{ab}$ is an interval of half-width $\varepsilon\Lambda_*^{(ab)}$, the min-dual-norm constant for the difference functional), including C9's phenomenon: the **difference** — hence possibly the sign — can be identified when neither value is.

**Proposition DE-R2 (marginals cannot decide rankability). [proved]**
Overlap of the marginal intervals $J_{x_a},J_{x_b}$ neither implies nor refutes sign identifiability. Witnesses: the shifted diagonal $\{(t,\,t-\tfrac14):t\in[\tfrac14,1]\}$ has heavily overlapping marginals $[\tfrac14,1]$ and $[0,\tfrac34]$ yet $\Delta_{ab}=\{\tfrac14\}$ — sign identified; while $J^{\mathrm{sq}}=[0,1]^2$ has overlapping marginals and $\Delta_{ab}=[-1,1]$ — sign unidentified. Rankability is a property of the **joint** pushforward (DE-J6). $\square$

---

## 2. Deciding when both signs remain admissible

Suppose case (iii), with both strict signs realized by admissible members.

**Proposition DE-R3 (what the set alone yields). [proved]**
Under $0$–$1$ loss on $\{a\!\succ\!b,\ b\!\succ\!a\}$:
(i) both strict actions have worst-case loss $1$ — minimax is **indifferent**; the identified set produces a tie, not a selection;
(ii) minimax regret is the same tie (prediction collapse, DE-H6(ii) with the $0$–$1$ table);
(iii) if abstention with constant loss $c\in(0,1)$ is admitted, minimax uniquely selects **abstain** — the only set-canonical strict preference available in case (iii);
(iv) by DE-S4 (applied to the sign-symmetric situation: the map $\delta\mapsto-\delta$ swaps the actions and, when it fixes $\Delta_{ab}$, forces any intrinsic equivariant rule into indifference or abstention), **no canonical rule strictly prefers one ordering** while both remain admissible. $\square$

**Theorem DE-R4 (the minimal object for a strict preference). [proved]**
A strict selection between the two orderings that satisfies R0–R2 (DE-P2) is rationalized by a monotone preorder on the two loss profiles; under $0$–$1$ loss the profiles are the indicators of the events $E_+=\{\delta>0\}$, $E_-=\{\delta<0\}$ (as subsets of $\Delta_{ab}(O)$, with $\{\delta=0\}$ split by the declared tie convention), and the preorder reduces to a **comparative probability on the sign events** — a relation "$E_+\ \trianglerighteq\ E_-$" — i.e. one bit of comparison; quantifying its regret requires one number
$$p_{ab}\;=\;w(E_+)\in[0,1]$$
for a declared weighting $w$: expected $0$–$1$ loss of $a\!\succ\!b$ is $1-p_{ab}$ (tie convention absorbed), so the Bayes action is $a\!\succ\!b$ iff $p_{ab}>\tfrac12$. **This single number (or an ambiguity interval for it) is the entire decision object of pairwise ranking under $0$–$1$ loss.** For graded ranking losses $\ell(|\delta|)$ on errors, the object grows to the law of $\delta$ — still one-dimensional, still a pushforward of the weighting on $I(O)$. $\square$

---

## 3. Using history: legitimacy and second-order partial identification

**The population assumption, exactly (mandated).** A strict preference justified by historical members requires:
$$\textbf{(EXCH)}\ (f_1,\dots,f_n,f_\beta)\ \text{exchangeable};\qquad \textbf{(LIK)}\ \text{declared noise law (for conditioning on }O\text{)};\qquad \textbf{(COV)}\ \text{each historical sign estimable from that member's own data}.$$
Under (EXCH)+(LIK), $p_{ab}=\Pr(f_\beta(x_a)>f_\beta(x_b)\mid O)$ is well-defined and — by DE-H2 — computed from a posterior supported inside $I(O)$: the preference tilts within the admissible set and never overrides case (i)/(ii) identification. Without (EXCH), Part IX's impossibility applies verbatim: historical sign frequencies carry no information about the current member's sign.

**Theorem DE-R5 (second-order partial identification of $p_{ab}$). [proved]**
(COV) is itself gated by the frozen identification theory: historical member $f_i$'s sign at $(x_a,x_b)$ is known only if $f_i$'s own design identifies it — e.g. $x_a,x_b$ (or the difference functional) covered per Thm 7.1/F17; off its covered set, F18's dichotomy applies to $f_i$ too. Let each historical member contribute an *interval* $[\,l_i,u_i\,]\ni\chi_i:=\mathbf 1\{\delta_i>0\}$ with $l_i=\mathbf 1\{\Delta^{(i)}_{ab}\subseteq(0,\infty)\}$, $u_i=\mathbf 1\{\Delta^{(i)}_{ab}\cap(0,\infty)\ne\emptyset\}$. Then the population frequency is only **interval-identified**:
$$p_{ab}\in\Big[\mathbb E\,l,\ \mathbb E\,u\Big],\qquad \widehat{[\,\cdot\,]}_n=\Big[\tfrac1n\textstyle\sum l_i,\ \tfrac1n\sum u_i\Big],$$
with interval width $=$ the population fraction of members whose sign their own data do not identify — **partial identification recurs at the decision-object level**: the object that was to resolve the current member's ambiguity is itself only set-identified, by the same theory, one level up. More historical members shrink the *sampling* error of the endpoints (Part VIII), never this width; only better per-member coverage does. $\square$

**Theorem DE-R6 (three-tier ranking hierarchy). [proved]**
Every ranking output belongs to exactly one tier, each consuming strictly more declared content, each to be labeled in the ledger:
- **Tier 1 — sign-identified** (case (i)/(ii)): correctness certified worst-case; no $\Delta$ consumed.
- **Tier 2 — decision-robust:** sign not identified, but every weighting in the declared ambiguity class $\mathcal Q$ (e.g. the estimated interval $[\hat p^-,\hat p^+]$ of DE-R5 with its confidence margin) yields the same Bayes ordering — i.e. $\tfrac12\notin[\hat p^-,\hat p^+]$. Consumes (EXCH)+(LIK)+(COV)+class declaration; robust to everything inside $\mathcal Q$.
- **Tier 3 — tie-broken:** $\tfrac12$ inside the interval; a strict output exists only under a single declared law (or declared tie-break); the choice is an axiom, not an inference.
Strictness: Tier 1 $\subsetneq$ Tier 2 $\subsetneq$ Tier 3 witnesses — diagonal-shift family (Tier 1); $J^{\mathrm{sq}}$ with history concentrated on $\delta>0$ members (Tier 2, not 1); $J^{\mathrm{anti}}$ with symmetric history (Tier 3, not 2). $\square$

---

## 4. Summary

$$\boxed{\begin{array}{c}\text{Ranking splits exactly into: identifiability of the sign — a joint-pushforward property, settled by the frozen theory (Thm 7.1/F18);}\\ \text{and preference between admissible signs — requiring precisely a comparative probability }p_{ab}\text{, legitimate only under (EXCH)+(LIK)+(COV),}\\ \text{itself generally only interval-identified from partially identified historical members.}\end{array}}$$
