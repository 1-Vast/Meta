# Loss-Typed Robust Information Floor — REPAIRED (Part I)

> **Status:** Phase-8.1, 2026-08-03. Supersedes `../08_decision_operator_realizability/loss_typed_information_floor.md`. Repairs audit targets **T1** (outer-envelope certificate semantics) and **T3** (abstention loss-typing); all other content carried over. Results: **DR-F0–F3** carried (with F1(iii) repaired as **DR-F1-R**), **DR-F4** replaced by **DR-F4-R**, new elementary device **DR-F5-R**.

---

## 1. The object (unchanged)

For nonempty $J\subseteq\mathbb R^m$, action set $\mathcal A$, loss $L:\mathcal A\times\mathbb R^m\to[0,\infty]$:
$$R_{\mathrm{set}}(J,\mathcal A,L)=\inf_{a\in\mathcal A}\sup_{v\in J}L(a,v),\qquad R_{\mathrm{rand}}(J,\mathcal A,L)=\inf_{\xi\in\Delta(\mathcal A)}\sup_{v\in J}\int L\,d\xi\ \le\ R_{\mathrm{set}}.$$

**DR-F1-R (specializations, abstention repaired). [proved]**
(i) Scalar absolute error: $R_{\mathrm{set}}=\tfrac12(\sup J-\inf J)$ for bounded $J$ — the frozen radius recovered. (ii) Sup-loss: coordinatewise radii (frozen). (iii) *Pairwise ranking, loss-typed:* with $\mathcal A=\{a\!\succ\!b,\ b\!\succ\!a\}$ and $0$–$1$ loss, both signs admissible: $R_{\mathrm{set}}=1$, $R_{\mathrm{rand}}=\tfrac12$. **If $\mathcal A$ additionally contains a declared abstain action with declared constant loss $c$:** $R_{\mathrm{set}}=\min(1,c)$, $R_{\mathrm{rand}}=\min(\tfrac12,c)$. There is no loss-independent abstention value: every abstention statement in this program is a statement about the declared pair $(\mathcal A,L)$ — abstention is an *action*, and its cost enters the game like any other loss value. (iv) Structured actions: unchanged. $\square$

---

## 2. The floor theorem (unchanged — it always concerned the true set)

**DR-F2 [proved, carried verbatim].** For any rule $\Phi(O,H,\Delta)$, deterministic or randomized: $\sup_{f\in I(O)}\mathbb E\,L(a,e_Q(f))\ \ge\ R_{\mathrm{rand}}(J_Q(O),\mathcal A,L)$, and $\ge R_{\mathrm{set}}(J_Q(O),\mathcal A,L)$ for deterministic $\Phi$. **DR-F3 [carried]:** the floor moves only via new current-member evidence or declared member-level structure.

The floor is a property of the **exact** identified pushforward $J_Q(O)$. Nothing in DR-F2 refers to an approximation — the repair below is about what an *approximating system* may claim.

---

## 3. DR-F4-R: the three-type certificate discipline (replaces DR-F4)

Let $\widetilde J\subseteq J_Q(O)\subseteq\widehat J$. Monotonicity of $R_{\mathrm{set}}$ in its first argument (sup over a larger set, then inf — unchanged) gives the **bracket**
$$R_{\mathrm{set}}(\widetilde J)\ \le\ \underbrace{R_{\mathrm{set}}(J_Q(O))}_{\text{true floor, generally uncomputable}}\ \le\ R_{\mathrm{set}}(\widehat J).$$

**(a) Exact type.** $R_{\mathrm{set}}(J_Q(O))$ — the true minimax value: simultaneously the information floor (DR-F2) and the best achievable worst-case guarantee. The only object entitled to the word "floor".

**(b) Outer type — certified guarantee, NOT a floor. [proved]**
For outer $\widehat J$, define the **certified guarantee value** $G_{\mathrm{cert}}(\widehat J,\mathcal A,L)=R_{\mathrm{set}}(\widehat J,\mathcal A,L)$, attained by $\hat a\in\arg\min_a\sup_{v\in\widehat J}L(a,v)$. Then the *true* worst-case risk of $\hat a$ satisfies $\sup_{v\in J_Q(O)}L(\hat a,v)\le G_{\mathrm{cert}}$: playing $\hat a$ is guaranteed to do at least this well. $G_{\mathrm{cert}}$ is an **upper bound on the achievable value and on the true floor** — an achievability certificate. It is **false to present it as a lower information floor**: with $J=\{0\}$, $\widehat J=\{0,100\}$, absolute loss, $G_{\mathrm{cert}}=50$ while the action $a=0$ achieves true worst-case $0$ — the audit's counterexample, now the type-error witness of record. *Phase-8's DR-F4(i) interpretive sentence is retracted; its monotonicity computation stands.*

**(c) Inner type — certified floor lower bound. [proved]**
For inner $\widetilde J$ (a **witness set**: DR-F5-R), $R_{\mathrm{set}}(\widetilde J)\le R_{\mathrm{set}}(J_Q(O))$, and by DR-F2 applied to the witnesses (each witness is an admissible member, so the adversary may play it): **no rule can guarantee below $R_{\mathrm{set}}(\widetilde J)$** — a valid, computable floor statement. An inner set is *not* a feasibility certificate (it under-represents ambiguity and may never be reported as $\widehat J$); it is exclusively a floor-bounding device. The two certificate roles are disjoint and must never be swapped — that swap was exactly the Phase-8 failure.

**(d) Tightness (carried, re-typed).** If $d_H(\widehat J,J)\le h$ and $L$ is $\ell_v$-Lipschitz in $v$: $G_{\mathrm{cert}}\le R_{\mathrm{set}}(J)+\ell_v h$ — the *guarantee* is near-tight; and the reported bracket width $G_{\mathrm{cert}}-R_{\mathrm{set}}(\widetilde J)$ is the end-to-end tightness deficit, auditable per emission. Discrete-loss caveat unchanged (no Lipschitz transfer at order boundaries; validity of both bracket ends is unaffected).

**DR-F5-R (witness certificates — the elementary inner device). [proved]**
Any finite set of explicitly constructed members $f_1,\dots,f_r$ verified admissible under the declared closure class (e.g. solutions of the class-consistency system with the observed data, feasible by construction) yields $\widetilde J=\{e_Q(f_j)\}_{j\le r}\subseteq J_Q(O)$ and hence the valid floor bound $R_{\mathrm{set}}(\widetilde J)$. This uses only monotonicity and membership verification — no new theory. Two witnesses realizing both signs of a difference already certify the ranking floor $\min(1,c)$ (deterministic) at that pair. $\square$

---

## 4. Repaired summary

$$\boxed{\begin{array}{c}R_{\mathrm{set}}(J_Q(O),\mathcal A,L)\ \text{is the information floor (DR-F2) — exact-set only.}\\ \text{An approximating system emits the bracket: certified floor bound }R_{\mathrm{set}}(\widetilde J)\ \text{from inner witnesses, certified guarantee }G_{\mathrm{cert}}(\widehat J)\ \text{from the outer envelope.}\\ \text{Outer values are guarantees, never floors; inner sets are floor devices, never feasibility reports. All abstention values are computed from the declared }(\mathcal A,L).\end{array}}$$
