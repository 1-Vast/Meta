# Output Object Alignment and Calibration (Items 1, 4)

> **Status:** Phase-26.1 (final theorem alignment repair), 2026-08-03. Phases 0–24 provenance; this package supersedes `25_1_final_theory_closure` where they conflict, and is self-contained. Audit of record: `../26_final_theory_freeze_audit/FINAL_AUDIT.md` (`FINAL_THEORY_NOT_READY`; theorem-alignment defects only — scope, target, strong convexity, approximation idea accepted). The operator is **not** redesigned and the scope is **not** enlarged. New results carry **AL-** numbers; every symbol used is defined here or in `hypothesis_class_and_statistics.md`. Tagged **[definition] / [proved] / [declared]**.

---

## Item 1 — Route A: the declared operator output has no support restriction

**Decision AL-1 (Route A, preferred). [declared]**
The support restriction is **removed from the declared operator output**. The operator is
$$\boxed{\ \mathsf A(F,z)\ =\ K\big(B(z)\,F(z)\big)\ }\qquad(\text{a nonempty compact convex, }W_1\text{-closed subset of }(\Delta(V),W_1)),$$
and **all** metrics, calibration, and consistency statements below refer to this object. The audited mismatch — FC-13 declaring $K(B(z)F(z))\cap I(S)$ while FC-9/FC-11 calibrated the unrestricted $K(B(z)F(z))$ — is removed at the source: there is now one output object, the unrestricted class, and it is the one that is calibrated. The symbols $I(S)$, and any intersection with an identified set, **do not appear** in the declared output or in any retained theorem. (This is a scope *narrowing* of the emitted object, permitted; it enlarges nothing. Support-restricted reporting, if ever wanted, is a separate downstream operation outside the frozen theorem and is not claimed here.)

*Consequence:* $\mathsf A(F,\cdot)$ is well-defined and $W_1$-closed for every measurable $F:Z\to\Delta_m$ (FC-5 class map on a valid band $B(z)F(z)\in\mathbb B$), so no nonemptiness or post-restriction-stability obligation arises — the Route-B obligation the audit described is moot under Route A.

## Base objects (self-contained; the audit's undefined-symbol list, closed)

**AL-2 [definition].** Carried verbatim from the closure package and repeated here so this file is auditable alone: fixed deployment $\mathcal D=(z_H^0,B(\cdot),\Delta_m,\mu,h)$; compact statistic domain $(Z,d_Z)$; finite context map $\kappa:Z\to C_\kappa$; value interval $V=[a_{\min},a_{\max}]$, $D_V^{\mathrm{val}}=a_{\max}-a_{\min}$; law space $(\Delta(V),W_1)$; coefficient simplex $C=\Delta_m$; band rule $B(z)=[\beta_0(z)|\beta_1|\cdots|\beta_m]$, $\beta_0(z)=b^{\mathrm{pop}}_{\kappa(z)}$, anchors fixed, $\kappa_B=\sup_z\|B(z)\|_{\mathrm{op}}<\infty$; band polytope $\mathbb B$ (closed/open CDF-band convention ⇒ $W_1$-closed classes); class map $K$; loss $L$ (convex, $L_{\mathrm{Lip}}$-Lipschitz, bounded $\bar L$ in the band argument).

**AL-3 [definition] The three declared assumptions, in-folder (replacing external A-STAT/A-CONT references).**
- **(S-IID)** the tasks $T_1,\dots,T_N$ and the current task are i.i.d. draws from the observable task law $P_T$ (this package uses **only** the IID branch — Item 5).
- **(S-CONT)** the base conditional risk has a continuous version: $z\mapsto L_0(z,\beta):=\mathbb E[L(\beta,Y_T)\mid\zeta=z]$ is continuous on $Z$ uniformly in $\beta$, with declared modulus $\varpi_\ell$ (so $L_0$ is everywhere-defined and measurable).
- **(S-GRID)** the Route-B output grid on $V$ at mesh $h$ is fixed (part of $\mathcal D$).
No other assumption (no C-IID, no DE-T3, no external interpolation theorem) is invoked anywhere in this package.

## Item 4 — the calibration theorem for $\mathsf A$, with correct probability typing

**Definitions AL-4 [definition].** $\zeta=z(S_T,Q_T,\gamma)$, law $\mu_\zeta$ on $Z$; $J_\mu(z,p)=L_0(z,B(z)p)+\tfrac\mu2\|p\|^2$ ($\ge\mu$-strongly convex on $\Delta_m$); target $g^\star_\mu(z)=\arg\min_{p\in\Delta_m}J_\mu(z,p)$ (unique, everywhere-defined); population risk $R_\mu(F)=\mathbb E_\zeta[J_\mu(\zeta,F(\zeta))]$; excess risk $\mathcal E_\mu(F)=R_\mu(F)-R_\mu(g^\star_\mu)\ge0$; operator metric $d_{\mathbb M}(F,g^\star_\mu)(z)=d_H^{W_1}\big(\mathsf A(F,z),\mathsf A(g^\star_\mu,z)\big)$ — **now between the declared outputs $\mathsf A$, which are the unrestricted classes, so metric and codomain coincide**; transfer constant $D_V=D_V^{\mathrm{val}}\kappa_B$ (the $\|B(z)\|$ factor $\kappa_B$ included); $\Phi(t)=D_V\sqrt{2t/\mu}$, continuous, $\Phi(0)=0$.

**Theorem AL-5 (calibration of the declared operator). [proved]**
For every measurable $F:Z\to\Delta_m$,
$$\big\|d_{\mathbb M}(F,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(\mathcal E_\mu(F)\big)+2h.$$
*Proof.* Strong convexity: $\tfrac\mu2\|F-g^\star_\mu\|_{L^2(\mu_\zeta)}^2\le\mathcal E_\mu(F)$ (as in the closure package, unchanged). Pointwise transfer, now to the **declared** object: $d_H^{W_1}(\mathsf A(F,z),\mathsf A(g^\star_\mu,z))=d_H^{W_1}(K(B(z)F(z)),K(B(z)g^\star_\mu(z)))\le D_V^{\mathrm{val}}\|B(z)F(z)-B(z)g^\star_\mu(z)\|_{\mathbb B}+2h\le D_V^{\mathrm{val}}\kappa_B\|F(z)-g^\star_\mu(z)\|+2h=D_V\|F(z)-g^\star_\mu(z)\|+2h$ (CDF-band stability + linear assembly; no intersection step, so no nonexpansiveness gap). Minkowski in $L^2(\mu_\zeta)$: $\le D_V\sqrt{2\mathcal E_\mu(F)/\mu}+2h$. Every object is the declared $\mathsf A$; the theorem calibrates the stated codomain. $\square$
