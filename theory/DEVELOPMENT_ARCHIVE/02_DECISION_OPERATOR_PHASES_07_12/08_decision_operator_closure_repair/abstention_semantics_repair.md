# Abstention Semantics Repair (Closure Target 2)

> **Status:** Phase-8.2, 2026-08-03. Replaces the abstention clauses of the 8.1 contract (Failure 2–3), which the audit refuted with the triple (strict robust value $1$, abstain cost $c=2$, tolerance $T=\tfrac12$): the old clause mandated abstention with loss $2$ over a strict action with loss $1$. New results **DC-A1–A4**, tagged **[proved]**. Principle: **abstention is an action inside the game; failure is a report about the game.** The two never substitute for each other.

---

## 1. Abstention typed as an action

Declare $a_{\mathrm{abs}}\in\mathcal A$ with declared loss $L(a_{\mathrm{abs}},v)$ — constant $c$ or value-dependent $c(v)$; write $\mathcal A_{\mathrm{strict}}=\mathcal A\setminus\{a_{\mathrm{abs}}\}$. Under the declared criterion risk $\rho$ (worst-case form: $\rho(a)=\sup_{v\in J}L(a,v)$, outer-envelope surrogate $\hat\rho(a)=\sup_{v\in\widehat J}L(a,v)$):

$$R(a_{\mathrm{abs}})\;=\;\rho(a_{\mathrm{abs}})\;=\;\sup_{v\in J}L(a_{\mathrm{abs}},v)\quad(=c\ \text{for constant cost}).$$

Abstention enters the argmin like every other action. Nothing else about it is special.

## 2. When abstention is selected — the derived threshold theorem

**Theorem DC-A1 (criterion-optimal abstention; the only valid "threshold rule"). [proved]**
Under the minimax criterion with constant declared cost $c$:
$$a_{\mathrm{abs}}\in\mathcal A^*\iff c\ \le\ R_{\mathrm{set}}(J,\mathcal A_{\mathrm{strict}},L),\qquad \{a_{\mathrm{abs}}\}=\mathcal A^*\iff c\ <\ R_{\mathrm{set}}(J,\mathcal A_{\mathrm{strict}},L).$$
*Proof.* $\rho(a_{\mathrm{abs}})=c$; the best strict value is $R_{\mathrm{set}}(J,\mathcal A_{\mathrm{strict}},L)$; compare. $\square$
So a threshold rule is legitimate **exactly when** the "uncertainty" it thresholds is the *loss-typed strict-action robust value* and the threshold is the *declared abstention cost* — i.e. when it is nothing but the argmin comparison. Realizable form: with the outer surrogate, "abstain iff $c\le\hat R_{\mathrm{strict}}:=R_{\mathrm{set}}(\widehat J,\mathcal A_{\mathrm{strict}},L)$" selects abstention only when it is optimal for the certified problem; since $\hat R_{\mathrm{strict}}\ge R_{\mathrm{strict}}$, this errs (if at all) toward abstaining — conservative in the guarantee sense, and every emitted number remains a valid outer risk.

**Corollary DC-A2 (forbidden threshold rules — with witnesses). [proved]**
(i) *Wrong uncertainty quantity:* "abstain iff identification radius $>$ threshold" is not optimal in general — witness: scaled diagonal $J=\{(t,t):t\in[0,M]\}$, $0$–$1$ ranking loss, any $c>0$: value radius $M/2$ arbitrarily large while the ordering (tie) is identified, strict robust value $0<c$ — the radius rule abstains where acting is free. (ii) *Wrong threshold:* any threshold $T\ne c$ compared against the strict robust value — the audit's triple: strict value $1$, $c=2$, rule thresholded at $T=\tfrac12$ abstains and incurs $2>1$; conversely $c=\tfrac14$, $T=\tfrac12$, strict value $\tfrac13$: rule acts ($\tfrac13<T$) though abstention at $\tfrac14$ is strictly better. A rule of the shape "uncertainty $>$ threshold $\Rightarrow$ abstain" is therefore **forbidden unless it provably coincides with DC-A1's argmin comparison** — which pins both the quantity and the threshold. $\square$

**Remark (Bayes/Γ-minimax criteria).** DC-A1 generalizes verbatim: $a_{\mathrm{abs}}\in\mathcal A^*$ iff its (class-)expected declared loss is minimal; for value-dependent $c(v)$, $\rho(a_{\mathrm{abs}})=\sup_{\mu\in\mathcal Q}\int c\,d\mu$. In every case abstention is *selected by the criterion*, never imposed on it.

## 3. When failure must be reported instead

**Definition DC-A3 (in-game vs out-of-game).** *Abstention* presupposes a well-posed decision problem: nonempty certified $\widehat J$, well-defined $\rho$, declared $(\mathcal A,L)$ including $a_{\mathrm{abs}}$'s cost. *Failure* is the report that the problem itself is broken or the demanded emission is untypable. Failure has no loss value and is not an element of $\mathcal A$; pricing it would smuggle epistemic breakdown into the game.

**Theorem DC-A4 (the exhaustive split). [proved]**
Every non-strict outcome is exactly one of:
- **Abstention (in-game):** $a_{\mathrm{abs}}$ is criterion-optimal (DC-A1). Emitted as the selected action with its declared cost in the guarantee row.
- **Failure (out-of-game), exactly the cases where DC-A1's comparison cannot legitimately run or its result cannot be emitted:**
  1. empty section / unrealizable support (misspecification): $\rho$ is a sup over the empty set — undefined; frozen flag;
  2. invalid certificates (envelope-containment or class-coverage self-check fails): $\hat\rho$ is not a certified upper bound — comparisons using it are void;
  3. empty conditioned/transported class: the declared population assumptions are jointly inconsistent with the data — surface, do not renormalize;
  4. untypable demand: single-valued strict output demanded at a symmetric tie with no declared $\tau$ (DR-S3) — emit the set or fail, never a hidden pick;
  5. **tolerance infeasibility:** a declared tolerance $T$ with $\min_{a\in\mathcal A}\hat\rho(a)>T$ (abstention included). No action — abstention included — meets the requirement; the honest emission is the infeasibility report, optionally accompanied by the argmin labeled "tolerance unmet". Converting this into abstention is exactly the refuted 8.1 clause: it can select an action (cost $c$) strictly worse than the best strict action, and it misrepresents "cannot meet the demand" as "chose not to act". $\square$

Tolerance is thereby re-typed: $T$ is a **feasibility requirement on the emission**, checked *after* the criterion has selected; it never participates in selection. (If a user wants tolerance-driven abstention, the legitimate encoding is to *declare* $c\le T$ — then DC-A1 delivers it as genuine optimality.)

---

$$\boxed{\ \text{Abstain iff }a_{\mathrm{abs}}\text{ wins the declared criterion (DC-A1) — the threshold is the declared cost against the loss-typed strict robust value, derived not posited.}\atop\text{Fail iff the problem or the demanded emission is broken (DC-A4 cases 1–5); tolerance violations are infeasibility reports, never abstentions.}\ }$$
