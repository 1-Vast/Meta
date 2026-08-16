# Target Definition (Route B: Risk-Optimal Meta-Learning Operator)

> **Status:** Phase-21 (target alignment repair), 2026-08-03. Phases 0–19.5 unmodified; Phase-20's alternation between two targets is **abandoned, not patched**. Audit of record: `../21_trainable_operator_audit/` (`TRAINABLE_OPERATOR_FOUNDATION_INVALID`; sole obstruction: two targets, never proved identical). **Route B is chosen; Route A is not used.** New results carry **PT-** numbers, tagged **[proved] / [declared]**. No architecture, no application, no implementation.

---

## 0. The single decision, stated once

**The learner has exactly one target: the risk-optimal coefficient map $g^\star:Z\to C$ defined below.** The canonical observable operator $A^\star$ of the earlier phases is **not** a target of this phase, is **not** used as an imitation label, and appears nowhere in the objective. The only residue of prior "canonical" objects is the fixed population band $b^{\mathrm{pop}}$, which enters solely as one *constant input* to a fixed affine assembly — never as a thing the learner is asked to match. This eliminates target-switching by construction: there is one object, and it is defined by the task risk.

## 1. Fixed decoder structure (so the target is a map $Z\to C$)

Under the declared skeleton (SKEL): $Z$ is a compact metric statistic domain; the **coefficient space**
$$C=[0,1]\times\Delta_{m-1}\quad\text{(mixing weight }\lambda\ \times\ \text{simplex weights }w\text{ over }m\text{ declared anchors)}$$
is compact convex. The **anchor bands** $b_1,\dots,b_m\in\mathbb B$ are **declared constants** (SKEL), chosen so their convex hull together with $b^{\mathrm{pop}}$ is all of $\mathbb B$ (e.g. the vertices of the polytope $\mathbb B$, finite by SKEL). The **assembly** is affine in the coefficient:
$$\mathsf{asm}(c;z)=(1-\lambda)\,b^{\mathrm{pop}}_{\kappa(z)}+\lambda\textstyle\sum_j w_j b_j,\qquad c=(\lambda,w),$$
continuous, and surjective onto $\mathbb B$ (its image is $\operatorname{conv}(\{b^{\mathrm{pop}}\}\cup\mathbb B)=\mathbb B$). The full operator value is $K(\mathsf{asm}(c;z))$ restricted to $\mathrm{supp}\,I(S)$, with $\omega$-invariant confidence/rung/certificate channels. **Because anchors are fixed, the entire learnable object is the coefficient map $F:Z\to C$** — matching the mandate's type exactly, with no expressiveness lost ($\mathsf{asm}(C;z)=\mathbb B$).

## 2. The operative local risk and the target

Let $T\sim P(T)$, $\zeta=z(S_T,Q_T,\gamma)\in Z$ the induced random statistic, $A_T$ the task's **identified** target information (observable; point channel per Phase-19.5). Define the **local risk** at a statistic $z$ and coefficient $c$:
$$\ell_0(z,c)=\mathbb E\big[L\big(\mathsf{asm}(c;z),\,A_T\big)\ \big|\ \zeta=z\big],\qquad \ell(z,c)=\ell_0(z,c)+\tfrac{\mu}{2}\|c\|^2,$$
where $L$ is the declared convex band-score loss and $\tfrac\mu2\|c\|^2$ is the **declared strongly-convex term** (`assumptions.md` A-SC: either the base score is already $\mu$-strongly convex in $c$ and $\mu$ is its modulus, or a ridge with declared $\mu>0$ is adjoined — the target is defined with respect to this operative risk, and that choice is owned, not hidden).

**Definition PT-1 (the target).**
$$\boxed{\ g^\star(z)\ =\ \operatorname*{arg\,min}_{c\in C}\ \ell(z,c)\ }\qquad(z\in Z).$$

**Theorem PT-2 (the target is a well-defined function). [proved]**
Under (A-STAT), (A-LOSS), (A-SC): for each $z$, $\ell(z,\cdot)$ is continuous and $\mu$-strongly convex on the compact convex $C$, hence attains a **unique** minimizer; $g^\star:Z\to C$ is therefore a single-valued function (not a correspondence). Its continuity — the property Phase 20 assumed without proof — is *derived* in `approximation_theorem.md` (PT-6), not posited here.
*Proof.* $c\mapsto L(\mathsf{asm}(c;z),A_T)$ is convex ($L$ convex in the band argument, $\mathsf{asm}$ affine in $c$); expectation preserves convexity, so $\ell_0(z,\cdot)$ is convex and finite (bounded loss); adding $\tfrac\mu2\|c\|^2$ gives $\mu$-strong convexity; continuity in $c$ from $L$-Lipschitz-in-band + affine $\mathsf{asm}$. Strongly convex + continuous on compact convex ⇒ unique minimizer (existence by Weierstrass, uniqueness by strict convexity). $\square$

## 3. Why this is the risk-optimal operator, not a surrogate

**Proposition PT-3 (Bayes optimality over all measurable maps). [proved]**
The population risk of any measurable coefficient map $F:Z\to C$ is $R(F)=\mathbb E_\zeta[\ell(\zeta,F(\zeta))]$ (pointwise separability, proved in `calibration_theorem.md` PT-7), so
$$R(F)\ \ge\ \mathbb E_\zeta\big[\min_{c}\ell(\zeta,c)\big]\ =\ \mathbb E_\zeta[\ell(\zeta,g^\star(\zeta))]\ =\ R(g^\star),$$
with equality iff $F=g^\star$ $\mu_\zeta$-a.e. **$g^\star$ is exactly the minimizer of the declared population risk over all measurable maps** — the risk-optimal meta-learning operator, defined by the task objective and nothing else. $\square$

**Ownership statement (honesty on the regularizer).** When (A-SC) is secured by an adjoined ridge, $g^\star=g^\star_\mu$ is the minimizer of the *regularized* risk; it is not claimed equal to the unregularized minimizer $g^\star_0$. The bias $\|g^\star_\mu-g^\star_0\|$ and the trade $\mu\downarrow$ (bias$\downarrow$, calibration constant$\uparrow$) are stated and bounded in `failure_modes.md` (PT-12). The target is whichever the deployment declares $\mu$ for; there is still exactly one target per declared objective.
