# Derived Approximation Theorem (§2–§3)

> **Status:** Phase-15, 2026-08-03. The approximation property is **derived, not assumed**: P-CAP is discharged as a corollary of a constructive theorem about the concrete family of `finite_parameter_family.md`. The four assumption layers are stated separately and never merged (mandate §3). New results **DM-5–DM-9**, tagged **[proved] / [declared]**.

---

## §3 first: the four labeled layers

- **A. Representation assumption.** (FIN-ATLAS) + the factorization $A_\phi=\Pi_{\mathrm{can}}\circ g^\star\circ E$ on the compact stratified statistic domain $Z$. *Status:* (FIN-ATLAS) declared; the factorization **proved** under it (DM-2). This is the only place the audit's counterexample bites, and it is quarantined here.
- **B. Continuity assumption.** $g^\star$ is continuous on $Z$. *Status:* **proved**, stratum-wise, in the strong form: on each stratum $g^\star$ is $\mathrm{clip}\circ(\text{affine})$ — piecewise affine with Lipschitz constant $1$ in the frequency coordinates (endpoints are frequencies shifted by per-stratum constants and clipped). No continuity across strata is needed or claimed: $Z$ is a disjoint union and the decoder reads the stratum label. The layer is stated separately because for *other* canonical constructions (Route B, §4) it must be re-verified, not because it is in doubt here.
- **C. Approximation theorem.** DM-5/DM-6 below — a proof, with explicit finite $p(\varepsilon)$, from A and B and nothing else.
- **D. Optimization assumption.** Split honestly (DM-8): *existence* of an optimal $\theta$ and *in-principle attainability* to any tolerance are **proved**; only the *efficiency* of a practical search is declared.

---

## §2: the theorem

**Theorem DM-5 (exact representation on the horizon). [proved]**
On every stratum with all fiber counts recorded exactly ($\sigma\in\{0,\dots,\bar N\}^{C_\kappa}$), the canonical endpoint map is $z\mapsto\mathrm{clip}_{[0,1]}(z\pm\eta_\sigma)$ with stratum-constant margins — piecewise affine in each coordinate with kinks only at the two clip thresholds. A grid containing the kink coordinates ($G$ chosen so that $\{\eta_\sigma,1-\eta_\sigma\}\subset\mathcal G$ per stratum, or the kinks added as nodes) makes multilinear interpolation **reproduce $g^\star$ exactly** (interpolation is exact on functions affine per grid cell in each coordinate). Hence there is $\theta^\star\in\Theta_p$ with
$$g_{\theta^\star}\big|_{\text{horizon strata}}\ =\ g^\star\big|_{\text{horizon strata}}\qquad\text{— zero approximation error, finitely parameterized.}\ \square$$

**Theorem DM-6 (derived uniform approximation — the P-CAP replacement). [proved under A]**
For every $\varepsilon>0$ choose the horizon $\bar N(\varepsilon)$ with the canonical margin $\eta(\bar N)\le\varepsilon$, and the grid as in DM-5 (adding mesh $1/G\le\varepsilon$ for the $\top$-strata). Let $\theta^\star$ take the exact values of DM-5 on horizon strata and the margin-zero map $z\mapsto z$ (grid values $=$ node coordinates; exactly representable, it is affine) on $\top$-strata. Then
$$\sup_{H}\ d_{\mathbb M}\big(A_{\theta^\star}(H),\ A_\phi(H)\big)\ \le\ \alpha\,\tfrac12\bar H\cdot\sup_H\big\|g_{\theta^\star}(E(H))-g^\star(E(H))\big\|_\infty\ \le\ \alpha\,\tfrac12\bar H\,\varepsilon,$$
because: horizon strata contribute $0$ (DM-5); $\top$-strata contribute at most the canonical margin $\eta(\bar N)\le\varepsilon$ (the canonical map differs from the margin-zero map by at most its own margin); confidence and rung coordinates contribute exactly $0$ (shared postprocessing, DM-2(iii)); and endpoint-to-operator transfer is the proved Hoffman lemma (PM-2a). The parameter count is explicit and finite:
$$p(\varepsilon)\ =\ (\bar N(\varepsilon)+2)^{|C_\kappa|}\cdot\big(G(\varepsilon)+1\big)^{\,2|\mathcal E||C_\kappa|}\cdot 2|\mathcal E||C_\kappa|,\qquad \bar N(\varepsilon)=O\!\big(\varepsilon^{-2}\ln(|\mathcal E||C_\kappa|/\delta)\big),\ G(\varepsilon)=O(\varepsilon^{-1}).$$
Therefore $\inf_{\theta\in\Theta_{p(\varepsilon)}}\sup_H d_{\mathbb M}(A_\theta(H),A_\phi(H))\le\alpha\tfrac12\bar H\varepsilon$ — **a theorem with a witness, not a capacity postulate**. $\square$

**Corollary DM-7 (P-CAP discharged). [proved]** The Phase-13 assumption P-CAP holds for the concrete family $\{A_\theta\}_{\theta\in\Theta_{p(\varepsilon)}}$ by DM-6; every Phase-13 conditional result (PM-8/PM-9) now holds unconditionally for this family, with (FIN-ATLAS) as the surviving named hypothesis. The audit's objection — "the central approximation guarantee is P-CAP restated" — is closed: the guarantee is now derived from representation (A) and continuity (B) via an explicit constructive witness. $\square$

**Theorem DM-8 (optimization layer, split). [proved / declared as marked]**
(i) *Existence [proved]:* the objective $F(\theta)=\sup_H d_{\mathbb M}(A_\theta(H),A_\phi(H))=\sup_{z\in Z}(\text{per-stratum endpoint error transferred})$ is a supremum of $\theta$-Lipschitz functions (DM-4(i) + PM-2a), hence Lipschitz on the compact $\Theta_p$; its minimum is attained (Weierstrass).
(ii) *In-principle attainability [proved]:* a $\gamma$-net of $\Theta_p$ (finite, size $(\lceil L/\gamma\rceil)^p$) contains a $\gamma$-optimal parameter; evaluating $F$ to accuracy $\gamma$ needs only a finite grid of $Z$ (the per-stratum error is Lipschitz in $z$). So for every $\gamma>0$ a finite procedure outputs $\theta$ with $F(\theta)\le\min F+\gamma$ — $\gamma^{\mathrm{opt}}$ is achievable by proof, not by hope.
(iii) *Efficiency [declared]:* that a *practical* search attains (ii)'s tolerance faster than net enumeration is the only optimization content left as an assumption — labeled **D**, echoed, and deliberately outside this program (it is the sole residue, and it is not a correctness assumption: (ii) already secures the mathematics). $\square$

**Theorem DM-9 (total error, unconditional form for the concrete family). [conditional only on the declared statistical stack + (FIN-ATLAS)]**
With $\theta_N$ any $\gamma_N$-optimizer from DM-8(ii):
$$d_{\mathbb M}\big(A_{\theta_N}(H_N),\ M^\dagger\big)\ \le\ \underbrace{\alpha\tfrac12\bar H\,\varepsilon_N+\gamma_N}_{\text{approximation+optimization (proved achievable)}}+\underbrace{\alpha\tfrac12\bar H(\eta_N+\rho)+\beta\delta_N}_{\text{statistical (PM-5)}}\ \big[+(\alpha{+}\beta{+}\gamma)\ \text{on}\ \mathrm{Miss}_N\cup\text{bad, with their probabilities}\big],$$
all schedules provably realizable, all terms separately attributed; $\to0$ a.s. under the Phase-13 conditions with $\rho=0$. The four tiers (existence / identification / statistical / approximation) each now rest on their own proofs; no tier cites another; and the approximation tier is, for the first time, a **proved** tier. $\square$
