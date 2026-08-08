# Approximation Limits

> **Status:** Phase-5 derivation, 2026-08-02. Question: when can a finite-parameter approximator converge to $\mathbb A$? All results refereed; the instability analysis (CR-7) is the *corrected* version — the referee found a genuine quantifier error in the original and constructed the witness that makes the repaired theorem unconditional.

**Declared correctness metric (repaired per plan audit).** Sup over realizable covered configurations on compact subsets of the **relative interiors of all cells** — not only full-dimensional ones. (DM-6's regularity is per-cell, so this is available; without the extension the metric never tests the center where validity regions are null — stratum (i) with ambient $d>k$ — and the positive theorem would be silently vacuous there.) Envelope slack is measured one-sidedly. Worst-case only; per-query (sup-loss) only.

---

## 1. The master error decomposition (organizing display)

$$\text{prediction error}\ \le\ \underbrace{\tfrac12\,\omega\big(2\varepsilon+2h\big)+h}_{\substack{\text{representation error }h\\ \text{inside the modulus; nonlinear; can be }\infty}}\ +\ \underbrace{\eta_s}_{\substack{\text{selection slack}\\ \text{constant }1\text{ (tight)}}}\ +\ \underbrace{\big[\tfrac12 J\big]_{\text{transitions}}}_{\substack{\text{half-jump, continuous models}\\ \text{(sharp; localized)}}}\ \ (+\ \text{collar flips: CR-4}).$$

Representation error is governed by the $(m,L)$-entropy tradeoff (CR-9): $N(\mathcal C,2h)\lesssim(Lr/h)^m$ — impossibility statements bound the *pair* (dimension, Lipschitz budget), never dimension alone (a finite-parameter family with unbounded decoder constants evades bare-$m$ floors; refereed correction).

---

## 2. The positive theorem (CR-6)

**Theorem CR-6 (certified differentiable convergence).** Hypotheses: declared class, tame by route (A) or (B) (CR-1); **stability** on the compact configuration set $K$ (modulus locally bounded — the local, data-realizable form); exact archive. Then finite-parameter approximators converge to $\mathbb A$ on $K$ in the declared metric, differentiably and with valid certificates:
1. each branch is $C^r$ on its cell (DM-6(b)); on compact $K$ inside a cell's relative interior, polynomials approximate the branch in the $C^1$ norm — by the **simultaneous Weierstrass/Bernstein approximation theorem** (extend $C^1$-ly, mollify, Bernstein-approximate on an enclosing cube; derivatives converge too — this, not Stone–Weierstrass, is the correct attribution; refereed);
2. **outer rounding:** report $(\widehat{\text{lower}}-\eta,\ \widehat{\text{upper}}+\eta)$ with $\eta\ge$ the certified sup-error on $K$: enclosures stay outer — valid, conservative — at exactly $+\eta$ radius cost; $\eta\to0$ along the approximation sequence;
3. **margin guarding (CR-4):** selector flips are confined to the definable collar $\{\mu<2h\}$, flagged one-sidedly; the collar shrinks to the transition set as $h\to0$;
4. **flag realizability:** for effectively definable (e.g. semialgebraic) classes the flag sets are finite Boolean combinations of sign conditions on finitely many $C^r$ functions; one-sided approximation inner-approximates validity sets and outer-approximates their complements, confining misclassification to a definable collar. *(Effectivity scoping: cell data are computable for semialgebraic classes (CAD); nothing in the corpus makes DM-6 cells effective for general o-minimal structures — an explicit assumption for any engineering phase.)*

So convergence — including of first derivatives — holds cell-wise with certificates intact; this is the strongest convergence the impossibility results below permit.

---

## 3. Impossibility results

**I-1 (jumps; established).** No sequence of end-to-end continuous models converges uniformly on any neighborhood of a genuine jump (MP-4); soft or continuously randomized selection does not escape (CR-3). The lawful alternatives are exactly: discrete flags, or certificate inflation by $\sim J/2$ on the band.

**I-2 (instability trichotomy — CR-7; corrected and completed by the referee).** Define the *data-localized modulus* $\omega_{x,D,\tilde y}(t)$ (the modulus restricted to members $t$-consistent with the given data). The refereed facts:

(a) **The original global claim is false.** Global instability ($\omega_{x,D}(t)=\infty$ for all $t>0$) does *not* void certificates at fixed interior data: for the tanh family, a system with representation error $h$ can report the hull of the $(\varepsilon+h)$-local section inflated by $h$ — finite and valid against the whole $h$-ball of candidate classes. The composed bound's infinity is a vacuous upper bound there; the forced radius is governed by the **local** modulus.

(b) **What global instability does force (uniform form).** For every $h>0$ there exist realizable data — in an $\sim(\varepsilon+h)$-band at the validity boundary — where the $h$-ball-forced radius is $+\infty$, while exact knowledge attains radius $0$ at $\varepsilon=0$: a **void band that never disappears** under any representation error, shrinking only as $h\to0$.

(c) **Pointwise impossibility (local form; unconditional via the referee's witness).** If instability is *local at a realizable configuration* — $\omega_{x,D,\tilde y}(0)=0$ but $\omega_{x,D,\tilde y}(t)=\infty$ for all $t>0$ — then any representation error $h>0$ forces radius $+\infty$ *at that configuration*: zero tolerance, pointwise. Witness class satisfying simultaneously (i) no finite-dimensional Lipschitz-decoder representation and (ii) local instability at realizable data (domain $\{x_1\}\cup[0,1]$): $\mathcal F=\{f:\ f|_{[0,1]}=\theta+g,\ g\ 1\text{-Lipschitz},\ g(0)=0,\ f(x_1)=-1/\theta,\ \theta\ge1\}\cup\{f_0\equiv0\}$. The Lipschitz factor gives infinite entropy (i); at $D=\{x_1\}$, $x=0$, $\tilde y=0$ (realizable via $f_0$) the exact section is $\{0\}$ but every $t$-section contains $\{\theta\ge1/t\}$ (ii). **Hence: no finite-parameter Lipschitz-budget system outputs valid nontrivial certificates at that configuration — unconditionally.**

(d) **Exactly representable unstable classes are not impossible — just honest.** A finite-dimensional unstable class (tanh) admits exact representation ($h=0$); its honest certificates are finite at interior data and blow up only toward the validity boundary. Instability voids *tolerance to representation error*, not representability.

**I-3 (entropy).** $(m,L)$-floors (CR-9): for infinite-dimensional classes, achievable $h$ is bounded below by entropy numbers at any fixed Lipschitz budget — convergence rates in $m$ inherit the class's entropy exponents.

**I-4 (sample/task floors; scoped).** Exact-archive scope; common-core pattern: clipped excess $\sum_j(k_j-d)_+\ge d(N-d)=\dim\operatorname{Gr}(d,N)$, necessary and matched (DM-7) — with DM-7(e)'s own counterexample barring any general-class extension of the per-task surcharge; beyond-core patterns OPEN; per-task support floors $k\ge d$ ($2d{+}1$ nonlinear-generic; at $k\le5$: nonlinear classes admitted only at $d\le2$ generically or per-configuration verified. Archive noise: OPEN — no bound may be assumed.

---

## 4. Convergence verdict

A finite-parameter approximator converges to $\mathbb A$ **iff** all of: declared tame class; stability at the configurations claimed (else the I-2 void band/pointwise voids); representation within the $(m,L)$-entropy budget; archive meeting the DM-7 floors; discontinuities handled by flags or priced inflation; metric read on relative interiors of cells with one-sided slack. Under exactly these conditions CR-6 delivers convergence with derivatives and valid certificates; deny any one and a named impossibility result (I-1…I-4) blocks it.
