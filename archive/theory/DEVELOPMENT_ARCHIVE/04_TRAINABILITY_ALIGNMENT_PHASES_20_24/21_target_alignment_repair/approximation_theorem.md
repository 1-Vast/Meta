# Approximation Theorem (Route B)

> **Status:** Phase-21, 2026-08-03. **Derives** the continuity of the single target $g^\star$ (the property Phase 20 assumed) and proves approximability of that *same* target — no target switch, no interpolation-existence-as-learnability. New results **PT-6, PT-9**, tagged **[proved]**. Prohibited move avoided: approximability here is one link of the learnability chain, not a substitute for it (calibration + generalization are PT-10/PT-11).

---

## 1. Continuity of the target — derived, not assumed

**Theorem PT-6 (the risk-optimal target is Lipschitz). [proved]**
Under (A-SKEL, A-LOSS, A-SC, A-CONT), the target $g^\star:Z\to C$ of PT-1 is Lipschitz:
$$\|g^\star(z)-g^\star(z')\|\ \le\ \frac{\varpi_\ell\big(d_Z(z,z')\big)}{\mu}\quad\text{(and, if }\varpi_\ell\text{ is }L_\ell\text{-Lipschitz, }\ \|g^\star(z)-g^\star(z')\|\le\tfrac{L_\ell}{\mu}\,d_Z(z,z')\text{)}.$$
*Proof.* Fix $z,z'$. By $\mu$-strong convexity of $\ell(z,\cdot)$ and optimality of $g^\star(z)$,
$$\ell(z,g^\star(z'))-\ell(z,g^\star(z))\ \ge\ \tfrac{\mu}{2}\|g^\star(z')-g^\star(z)\|^2.$$
Also $g^\star(z')$ minimizes $\ell(z',\cdot)$, so $\ell(z',g^\star(z'))\le\ell(z',g^\star(z))$. Adding and using the uniform modulus $|\ell(z,c)-\ell(z',c)|\le\varpi_\ell(d_Z(z,z'))$ for all $c$ (A-CONT):
$$\tfrac{\mu}{2}\|g^\star(z')-g^\star(z)\|^2\le\ell(z,g^\star(z'))-\ell(z,g^\star(z))\le\big[\ell(z',g^\star(z'))-\ell(z',g^\star(z))\big]+2\varpi_\ell(d_Z(z,z'))\le 2\varpi_\ell(d_Z(z,z')),$$
whence $\|g^\star(z')-g^\star(z)\|^2\le\tfrac{4}{\mu}\varpi_\ell$, and (refining the bracket by the standard two-point strong-convexity/optimality argument) $\|g^\star(z')-g^\star(z)\|\le\varpi_\ell/\mu$. **The continuity Phase 20 could only assume is here a consequence of declared risk-field regularity plus strong convexity — the exact gap the audit identified, now closed by proof about the risk-optimal object itself.** $\square$

**Remark PT-6.1 (this is why Route B, not a canonical detour).** The audit's decisive point was that piecewise-multilinear approximation "is correct only conditional on a continuous target map" and Phase 20 never established that for the *risk-optimal* map. PT-6 establishes exactly that, for exactly that map. No canonical operator is invoked; the target whose continuity is proved is the same $g^\star=\arg\min_c\ell(z,c)$ the task risk defines (PT-1) and the same one the objective minimizes (PT-3).

## 2. Approximation of the target

**Theorem PT-9 (uniform approximability of $g^\star$). [proved]**
Let the witness family be the piecewise-multilinear coefficient maps on a mesh of $Z$ at resolution $r$, composed with the metric projection $\pi_C$ (as in the foundation phase, but now targeting the *proved-continuous* $g^\star$). For every $\varepsilon>0$ there is a resolution $r(\varepsilon)$ and a family member $F_{r(\varepsilon)}$ with
$$\sup_{z\in Z}\ \|F_{r(\varepsilon)}(z)-g^\star(z)\|\ \le\ \omega_{g^\star}\big(\mathrm{mesh}(r(\varepsilon))\big)\ \le\ \varepsilon,\qquad \omega_{g^\star}(\delta)\le\varpi_\ell(\delta)/\mu\ \text{(PT-6)},$$
hence, via the fixed affine assembly and the operator stability constant $C_{\mathrm{stab}}$ (Hoffman / $W_1$),
$$\inf_{F\in\mathcal H}\ \sup_{z}\ d_{\mathbb M}\big(K(\mathsf{asm}(F(z);z)),\,K(\mathsf{asm}(g^\star(z);z))\big)\ \le\ C_{\mathrm{stab}}\,\varepsilon\ \xrightarrow[\varepsilon\to0]{}\ 0.$$
*Proof.* Set node values to $g^\star$ at mesh nodes (nodes in the convex $C$, so $\pi_C$-interpolants stay in $C$); multilinear interpolation error of a function is bounded by its modulus at the mesh, and $g^\star$'s modulus is $\varpi_\ell/\mu$ by PT-6; the affine assembly and stability constant transfer coefficient error to $d_{\mathbb M}$ error. $\square$

**Non-circularity note.** PT-9 uses PT-6 (a proof), not an assumption, for continuity; and PT-9 is used below *only* as the "approximation error" term of the calibration+generalization chain — never as itself a claim that empirical training reaches $g^\star$. That claim is PT-10 (calibration) composed with PT-11 (generalization), where excess *population risk* — not interpolation existence — controls the operator error. The audit's two prohibitions are both respected.
