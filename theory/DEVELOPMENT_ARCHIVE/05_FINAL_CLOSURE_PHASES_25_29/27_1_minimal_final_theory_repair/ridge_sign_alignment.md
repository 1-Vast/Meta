# Ridge Sign Alignment (Item 1)

> **Status:** Phase-27.1 (minimal final theory repair), 2026-08-03. Phases 0–24 provenance; the Phase-26.1 aligned theorems (AL-1–AL-12) are retained and only the three audited points are repaired. Audit of record: `../27_final_theory_audit/FINAL_AUDIT.md` (`FINAL_THEORY_NOT_READY`; three blockers). No theory redesign, no operator change, no scope enlargement. New results carry **MR-** numbers, tagged **[declared] / [proved]**.

---

## 1. The retained target uses the positive ridge — affirmed and disambiguated

**Declaration MR-1 (positive ridge is the intended and retained target). [declared]**
The single learning target is, and remains,
$$\boxed{\ g^\star_\mu(z)\ =\ \operatorname*{arg\,min}_{p\in\Delta_m}\Big[\ L_0(z,B(z)p)\ +\ \tfrac{\mu}{2}\|p\|^2\ \Big],\qquad \mu>0.\ }$$
This positive-ridge form is what the package defines (AL-4), what the empirical risk $\widehat R_{\mu,N}$ minimizes (AL-9), and what the calibration proof (AL-5) uses through the positive-$\mu$ strong-convexity lower bound. The present repair mandate confirms the positive ridge is intended; the audit's "sign conflict" was a discrepancy with a *prior mandate's* wording ($-\tfrac\mu2\|p\|^2$), not an internal inconsistency of the package. There is one target, and it carries the $+\tfrac\mu2\|p\|^2$ ridge.

## 2. The negative-ridge form is explicitly rejected

**Proposition MR-2 (why $-\tfrac\mu2\|p\|^2$ cannot be the retained target). [proved]**
Consider the negative-ridge functional $\tilde J(z,p)=L_0(z,B(z)p)-\tfrac\mu2\|p\|^2$. Then:
(i) *Strong convexity is destroyed.* $L_0(z,B(z)\cdot)$ is convex (convex loss ∘ linear assembly) but need not be strongly convex; subtracting $\tfrac\mu2\|p\|^2$ subtracts $\mu$ from every eigenvalue of the Hessian (where defined), so $\tilde J(z,\cdot)$ is at most $(-\mu)$-**weakly concave** in the worst case and generally **not convex**. Concretely, if $L_0(z,B(z)\cdot)$ is affine in $p$ on some direction (e.g. the zero-loss / linear-score cases already used as witnesses in this program), $\tilde J$ is *strictly concave* there.
(ii) *Existence/uniqueness fail.* A non-convex $\tilde J$ on the simplex may attain its minimum only at vertices, attain it non-uniquely, or have multiple local minima; the unique-everywhere-defined-minimizer guarantee (AL-4/AL-8) does not hold.
(iii) *Every downstream proof breaks.* The continuity modulus (AL-8), the calibration lower bound $\tfrac\mu2\|F-g^\star_\mu\|^2\le\mathcal E_\mu$ (AL-5), and hence the consistency chain (AL-12) all rely on the **positive**-$\mu$ strong-convexity inequality $J_\mu(z,p)-J_\mu(z,g^\star_\mu(z))\ge\tfrac\mu2\|p-g^\star_\mu(z)\|^2$. With the negative ridge the inequality reverses sign and none of these transfer.
**Therefore the negative-quadratic functional does not define the retained target: it destroys the strong-convexity guarantee on which existence, uniqueness, continuity, and calibration all depend.** The retained theory uses $+\tfrac\mu2\|p\|^2$ exclusively; $-\tfrac\mu2\|p\|^2$ is not an alternative form of the same theory but a different, ill-posed problem, and is excluded. $\square$

**Consequence.** With MR-1/MR-2, target identity holds unambiguously: one positive-ridge $g^\star_\mu$, and no reading of the theory admits the negative-ridge object. Blocker 1 is resolved by affirmation-plus-explicit-exclusion, exactly as the mandate directs — no target is redesigned, only disambiguated.
