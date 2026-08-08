# Target Definition (Item 1, 2, 4)

> **Status:** Phase-23.1 (target identity repair), 2026-08-03. Phases 0–21 unmodified; the operator is not redesigned and the scope is not changed. Audit of record: `../23_final_theory_freeze_audit/FINAL_AUDIT.md` (`THEORY_STILL_INVALID`; minimal obstruction: the ridge in barycentric coordinates changes the regularized minimizer, so the package held two targets). This phase fixes **one** target and rewrites every theorem to it. New results carry **TI-** numbers, tagged **[proved] / [declared] / [retracted]**.

---

## 0. Notation (fixed once, used everywhere below)

At a fixed deployment (Item 5, carried): $z_H=z_H^0$ constant; band matrix $B=[\beta_0\cdots\beta_m]$ with $\beta_0=b^{\mathrm{pop}}_{\kappa(z)}$ and fixed anchors $\beta_1,\dots,\beta_m$; ridge modulus $\mu>0$ **declared and fixed**; Route-B output value mesh $h$ **declared and fixed**; coefficient space $C=\Delta_m$. Base conditional risk of a band vector $\beta\in\mathbb B$:
$$L_0(z,\beta)\ =\ \mathbb E\big[\,L(\beta,A_T)\ \big|\ \zeta=z\,\big],\qquad \zeta=z(S_T,Q_T,\gamma).$$

## 1. The single learning target

**Definition TI-1 (the one target). [declared]**
$$\boxed{\ g^\star_\mu(z)\ =\ \operatorname*{arg\,min}_{p\in\Delta_m}\ \Big[\ \underbrace{L_0(z,Bp)}_{=\,\ell(z,Bp)}\ +\ \tfrac{\mu}{2}\|p\|^2\ \Big]\ =:\ \operatorname*{arg\,min}_{p\in\Delta_m}\ J_\mu(z,p).\ }$$
This is the *only* target in Phase 23.1. Every theorem below (Bayes optimality, approximation, calibration, consistency) refers to $g^\star_\mu$ and to nothing else. The regularized operative risk is $J_\mu(z,p)=L_0(z,Bp)+\tfrac\mu2\|p\|^2$; its population form is $R_\mu(F)=\mathbb E_\zeta[J_\mu(\zeta,F(\zeta))]$.

**Theorem TI-2 (well-defined, unique, everywhere-defined). [proved]**
Under (A-STAT, A-LOSS, A-SC via the linear assembly, A-CONT): $J_\mu(z,\cdot)$ is continuous and **at least $\mu$-strongly convex** in $p$ on the compact convex $\Delta_m$ (RP-2 corrected: "$\ge\mu$", not "exactly $\mu$", to cover a strongly-convex base — the audit's §3 note), so its minimizer exists and is unique for every $z$; A-CONT selects a continuous, hence measurable, everywhere-defined version. $g^\star_\mu:Z\to\Delta_m$ is a single-valued function on all of $Z$. $\square$

## 2. Retraction of the equivalence claim (Item 2)

**Retraction TI-3. [retracted]**
Every statement asserting that the Phase-22.1 barycentric target equals the Phase-21 target — in particular "the target is unchanged and only its coordinates changed" and "Bayes optimality transfers with the same risk" — is **withdrawn**. They are false when the ridge is active, by the audit's witness:
> one anchor, $L_0\equiv0$, $w=1$: old ridge $\tfrac\mu2(\lambda^2+1)$ minimized at $\lambda=0$; new ridge $\tfrac\mu2((1-\lambda)^2+\lambda^2)$ minimized at $\lambda=\tfrac12$. The assembled-band images coincide, but the regularized minimizers differ.

Consequently: the Euclidean ridge is **not** invariant under the nonlinear, non-bijective-at-$\lambda{=}0$ map $(\lambda,w)\mapsto p=(1-\lambda,\lambda w)$; the old regularized target $g^{\mathrm{old}}_\mu$ is **not** carried into this phase, is **not** claimed equal to $g^\star_\mu$, and is **not** referenced by any theorem here. There is one target, defined in barycentric coordinates by TI-1, and the ridge that defines it is $\tfrac\mu2\|p\|^2$ — period.

## 3. What $g^\star_\mu$ is, and is not (Item 4)

**Declaration TI-4. [declared]**
$g^\star_\mu$ is the **regularized risk-optimal target**: the minimizer of the *regularized* population risk $R_\mu$ over measurable maps $Z\to\Delta_m$ (proved in `bayes_optimality.md`). It is **not** claimed to equal the unregularized Bayes target $g^\star_0=\arg\min_p L_0(z,Bp)$; the two differ by the ridge bias, which is owned, not hidden. The deployment declares $\mu$ (and $h$, $z_H^0$, $B$); given that declaration, $g^\star_\mu$ is one fixed function. Choosing a different $\mu$ defines a different deployment with its own single target — never two targets within one deployment. The unregularized $g^\star_0$ is mentioned only to state, explicitly, that it is *not* the target.
