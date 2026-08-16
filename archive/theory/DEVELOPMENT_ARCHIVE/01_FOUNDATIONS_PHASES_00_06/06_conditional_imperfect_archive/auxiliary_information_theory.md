# Auxiliary-Information Theory (Part I)

> **Status:** Phase-6 (final theoretical phase), 2026-08-03. Phases 0–5 are frozen and cited, not modified. Sources: the frozen corpus (Theorem 1, CP-2, F17, DM-3, DM-7 — see `../RESEARCH_HANDOFF.md`). New results carry **CI-A** numbers. Every result is tagged **[proved]**, **[conditional on stated assumptions]**, **[impossible]**, or **[open]**. All CI-A claims were adversarially refereed; the corrections (fiber formulation is primitive; augmentation is class-conditional; realizability flag = section-emptiness) are incorporated.

**Additional datum.** Each member $f_a$ carries an observed auxiliary quantity $c_a\in C$, with $C$ an arbitrary set — **no relation between $c$ and $f$ is assumed**. For the new member, $c_b$ is known exactly; at most $k\le5$ evaluations are available.

---

## 1. The primitive object: the conditional (fiber) family

Do not posit a conditional distribution or a regression of $f$ on $c$. The weakest object that lets $c$ act is the **fiber**:

$$\mathcal F_c \;=\; \{\,f : (c,f)\in\mathcal F^+\,\},\qquad \mathcal F^+=\{(c_a,f_a)\}\subseteq C\times\mathbb R^{\mathcal X}.$$

Conditioning on $c_b$ replaces $\mathcal F$ by $\mathcal F_{c_b}$. This requires nothing of $C$: it is set-membership.

**Definition (conditional trace structure).** $T_{D,x}(c)=\{(f|_D,f(x)):f\in\mathcal F_c\}$; section $S_\varepsilon(\tilde y\mid c)$; **conditional modulus of ambiguity** $\omega_{x,D}(t\mid c)=\sup\{|f(x)-g(x)|:f,g\in\mathcal F_c,\ \|f|_D-g|_D\|_\infty\le t\}$.

**Proposition CI-A1 (conditional minimax identity). [proved]**
If $c_b$ is realizable ($\mathcal F_{c_b}\ne\emptyset$), then given $c_b$ and $\varepsilon$-noisy data the minimax error of estimating $f_b(x)$ is exactly $\tfrac12\,\omega_{x,D}(2\varepsilon\mid c_b)$, attained by the center of $S_\varepsilon(\tilde y\mid c_b)$.
*Proof.* Theorem 1 applied verbatim to the nonempty family $\mathcal F_{c_b}$; $c_b$ fixed-and-known makes the conditioning legitimate. $\square$

**The fiber view is primitive** (refereed): every base theorem that holds for *arbitrary* families — Theorem 1, CP-2 minimal sufficiency, the useful-iff theorem below — transfers verbatim to $\mathcal F_{c_b}$.

---

## 2. The augmentation reduction — and its exact scope

When $C\subseteq\mathbb R^m$ there is a second, complementary picture: adjoin $m$ virtual coordinates to the domain and set $f^+(\text{virtual}_j)=c_j$. Then **auxiliary information is $m$ exact evaluations at fixed virtual locations**, at heterogeneous noise (exact at virtual coordinates, level $\varepsilon$ at real ones).

**Proposition CI-A2 (heterogeneous-noise Theorem 1). [proved]**
Under per-coordinate consistency $|\tilde y_i-f(x_i)|\le\varepsilon_i$, the minimax error is $\tfrac12\sup\{|f(x)-g(x)|:|f(x_i)-g(x_i)|\le2\varepsilon_i\ \forall i\}$ — same proof (per-coordinate midpoint lower bound; per-coordinate triangle upper bound), with the convention $\varepsilon_i=0$ meaning hard equality. State the modulus via the constraint set, not a quotient metric.

**Scope of the augmentation reduction (refereed correction — the original "every base theorem applies on $X^+$" was overbroad):**
- **General-family theorems** (Theorem 1, CP-2, useful-iff CI-A3) transfer verbatim to $X^+$. **[proved]**
- **Class-conditional theorems** (F17 rank/archive, DM-7 counting) transfer **iff the augmentation preserves the class** — i.e. **iff $c_a=L(f_a)$ for a linear functional $L$**, so that $V^+=\{(f,Lf)\}$ is again a dimension-$d$ subspace with $G$ extended by virtual rows. For a nonlinear label the augmented family leaves the exactly-$d$ linear class (the synergy witness §4 augments to a line-through-origin $\cup$ affine-line, not a subspace). **[conditional: linear label]**
- **Budget theorems** (DM-3 truncation): virtual coordinates are cost-free/always-observed — effective budget $k+m$; a verbatim budget-$k$ application on $X^+$ is a false restriction. **[proved, with the accounting fix]**
- Verbatim heterogeneous Theorem 1 on $X^+$ (exact at virtual coords) yields the **global** minimax $\sup_c\tfrac12\omega(2\varepsilon\mid c)$, **not** the conditional identity CI-A1 — because it forces $c_f=c_g$ but leaves $c$ free. The conditional-on-$c_b$ identity is genuinely a fiber statement. This is why §1, not §2, is primitive.

---

## 3. When is $c$ useful? The useful-iff theorem

**Definition.** $c_b$ **provides no information at $(D,x)$** iff for every $\varepsilon\ge0$ and every data vector the conditional section equals the unconditional section (realizable sets included).

**Theorem CI-A3 (useful-iff). [proved]**
$c_b$ provides no information at $(D,x)$ **iff** the conditional joint window equals the unconditional one: $T_{D\cup\{x\}}(c_b)=T_{D\cup\{x\}}$.
Equivalently: **auxiliary information is useful at $(D,x)$ iff it changes the joint window.**
*Proof.* ($\Leftarrow$) the window determines all sections and, via $\operatorname{proj}_D$, all realizable sets. ($\Rightarrow$) the $\varepsilon=0$ sections reconstruct the window as $T_S=\bigcup_u\{u\}\times S_0(u)$. The $\varepsilon=0$ clause is **load-bearing** (refereed): witness $\mathcal F=\{f_q:f_q(x_1)=q,f_q(x)=0,q\in\mathbb R\}$, $c=\mathbf 1_{\mathbb Q}(q)$, $c_b=1$ gives fiber window $\mathbb Q\times\{0\}\ne\mathbb R\times\{0\}$ yet identical $\varepsilon>0$ sections everywhere (density) — an $\varepsilon>0$-only definition breaks the iff. $\square$
*Coverage of pure misspecification-detection:* if $c_b$ shrinks the realizable set without shrinking any nonempty section, some nonempty section becomes empty (projection), so window equality fails and the iff counts detection as usefulness. **[proved]**

**Monotonicity. [proved]** $\omega_{x,D}(t\mid c)\le\omega_{x,D}(t)$ always (sup over a sub-collection of pairs): conditioning never increases ambiguity.

---

## 4. The five mandated determinations

**(1) $c$ provides no information — condition. [proved]** Exactly CI-A3: $T_{D\cup\{x\}}(c_b)=T_{D\cup\{x\}}$. Sufficient concrete case: the label is independent of the joint window (e.g. $c$ a pure external tag with $\mathcal F_c$ ranging over the same window for all $c$).

**(2) $c$ strictly reduces ambiguity — condition. [proved]** $\omega_{x,D}(0\mid c_b)<\omega_{x,D}(0)$, equivalently the fiber window is a proper sub-graph-collapse of the full window at $(D,x)$.

**(3) Synergy — $c$ + evaluations identify what neither identifies alone. [proved]**
Witness (refereed, hypotheses load-bearing): members indexed by $\theta\in\mathbb R$, binary $c\in\{0,1\}$; $f_{\theta,0}=\theta g_0$, $f_{\theta,1}=\theta g_1$ with $g_0(x_1)=g_1(x_1)\ne0$, $g_0(x)\ne g_1(x)$, both nonzero at $x$. Then (i) $c$ alone ($k=0$): value at $x$ ranges over $\mathbb R$ — nothing; (ii) one evaluation at $x_1$ alone: determines $\theta$ but leaves the two candidate values $\{\theta g_0(x),\theta g_1(x)\}$ ($\theta\ne0$); (iii) $c$ **and** the evaluation: exact. Concrete: $\mathcal X=\{x_1,x\}$, $g_0=(1,1)$, $g_1=(1,2)$. Edge case $\theta=0$ collapses candidates to $\{0\}$ — correctly excluded.

**(4) Minimax quantification of the reduction. [proved]** Yes: the reduction is *exactly* the drop in the conditional modulus. Ambiguity radius falls from $\tfrac12\omega_{x,D}(2\varepsilon)$ to $\tfrac12\omega_{x,D}(2\varepsilon\mid c_b)$; the certified information supplied by $c_b$ at $(D,x,\varepsilon)$ is
$$\Gamma_c \;=\; \tfrac12\big(\omega_{x,D}(2\varepsilon)-\omega_{x,D}(2\varepsilon\mid c_b)\big)\;\ge\;0,$$
zero iff $c_b$ is useless in the worst case at that configuration (the conditional analogue of the Phase-3 adaptation gain).

**(5) Substituting an incorrect $c'$. [proved, with the flag fix]**
Cross-fiber modulus $\omega_{x,D}(t\mid c',c)=\sup\{|f(x)-g(x)|:f\in\mathcal F_c,\,g\in\mathcal F_{c'},\ \|f|_D-g|_D\|_\infty\le t\}$.
- (i) **Dichotomy.** Conditioning on $c'$ when the truth has $c_b$ either **fires the realizability flag** — defined as **emptiness of $S_\varepsilon(\tilde y\mid c')$** (refereed: the operational check, robust to non-closed $T_D(c')$ where a distance-threshold would silently pass an empty section) — or errs silently. In the silent case the true value's distance from the reported interval $\operatorname{hull}(S_\varepsilon(\tilde y\mid c'))$ is $\le\omega_{x,D}(2\varepsilon\mid c',c_b)$, and in fact the entire reported interval lies within $\omega_{x,D}(2\varepsilon\mid c',c_b)$ of $f_b(x)$ (any $\varepsilon$-consistent $g\in\mathcal F_{c'}$ has trace within $2\varepsilon$ of the truth's trace).
- (ii) **Harmless-iff. [proved]** Substitution is harmless for all data iff $T_{D\cup\{x\}}(c')=T_{D\cup\{x\}}(c_b)$ — the CI-A3 iff between fibers (the proof never used containment, so it holds between arbitrary fibers).
- (iii) **No metric bound without declared structure. [impossible]** There is no general bound relating the silent error to any "distance" between $c'$ and $c_b$: witness $C=\mathbb R$, $\mathcal F_c=\{f:f(x_1)=0,\ f(x)=h(c)\}$ with $h$ arbitrary — data near $0$ is silent under every fiber and the silent error is exactly $|h(c')-h(c_b)|$, defeating any $B(|c'-c_b|,\varepsilon)$. A metric bound requires declared structure on $C$ (e.g. a modulus of continuity for $c\mapsto T_{D\cup\{x\}}(c)$ in the operational/Hausdorff topology). The set-theoretic theory treats $c$ as a pure label.

---

## 5. The conditional useful-iff, in one line

$$\boxed{\ \text{Auxiliary information is useful at }(D,x)\ \iff\ \text{it changes the joint window }T_{D\cup\{x\}}\ }$$
— proved, sharp, and assumption-free on $C$. Its magnitude is $\Gamma_c$; its failure modes are exactly the two of (5). What it can never do without added structure on $C$ is grade the harm of a wrong label by a distance on $C$.
