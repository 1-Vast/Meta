# $L_p$ Definition and Probability Schedule (Items 2, 3)

> **Status:** Phase-27.1, 2026-08-03. Closes blockers 2 (undefined $L_p$) and 3 (probability schedule). New results **MR-3–MR-6**, tagged **[definition] / [proved] / [declared]**. Nothing else changed.

---

## Item 2 — the coefficient-loss Lipschitz constant $L_p$, defined

**Definition MR-3 ($L_p$). [definition]**
$$\boxed{\ L_p\ :=\ \sup_{z\in Z}\ \sup_{\substack{p_1,p_2\in\Delta_m\\ p_1\ne p_2}}\ \frac{\big|\,L_0(z,B(z)p_1)-L_0(z,B(z)p_2)\,\big|}{\|p_1-p_2\|}\ ,\quad\text{i.e. }\ \big|L_0(z,B(z)p_1)-L_0(z,B(z)p_2)\big|\le L_p\,\|p_1-p_2\|\ \ \forall z,p_1,p_2.\ }$$

**Proposition MR-4 ($L_p$ is finite, with an explicit bound). [proved]**
$L_p\le L_{\mathrm{Lip}}\,\kappa_B<\infty$, where $L_{\mathrm{Lip}}$ is the band-Lipschitz constant of the loss (AL-2) and $\kappa_B=\sup_z\|B(z)\|_{\mathrm{op}}$ (AL-2).
*Proof.* For fixed $z$: $|L_0(z,B(z)p_1)-L_0(z,B(z)p_2)|=|\mathbb E[L(B(z)p_1,Y_T)-L(B(z)p_2,Y_T)\mid\zeta=z]|\le L_{\mathrm{Lip}}\|B(z)p_1-B(z)p_2\|_{\mathbb B}\le L_{\mathrm{Lip}}\kappa_B\|p_1-p_2\|$ (loss Lipschitz in the band argument, then linear assembly). Take sup over $z,p_1,p_2$. $\square$

**Use in the approximation-to-risk conversion (the AL-12 step, now grounded). [proved]**
The excess-risk contribution of approximation is bounded via $L_p$: for the witness $F_{\omega^\star}\in\mathcal H_N$ of AL-7,
$$\inf_{\Omega_N}R_\mu-R_\mu(g^\star_\mu)\ \le\ R_\mu(F_{\omega^\star})-R_\mu(g^\star_\mu)\ \le\ \mathbb E_\zeta\big[L_0(\zeta,B(\zeta)F_{\omega^\star}(\zeta))-L_0(\zeta,B(\zeta)g^\star_\mu(\zeta))\big]+\tfrac\mu2\mathbb E_\zeta[\|F_{\omega^\star}(\zeta)\|^2-\|g^\star_\mu(\zeta)\|^2],$$
and the first expectation is $\le L_p\,\sup_z\|F_{\omega^\star}(z)-g^\star_\mu(z)\|\le L_p\,\varepsilon_{\mathrm{approx}}(N)$ by MR-3; the ridge difference is $\le\mu\,\mathrm{diam}(\Delta_m)\,\varepsilon_{\mathrm{approx}}(N)$ (Lipschitz of $\|\cdot\|^2$ on the compact $\Delta_m$). So the clean constant in AL-12 is $L_p':=L_p+\mu\,\mathrm{diam}(\Delta_m)$; writing $L_p$ for this combined coefficient-to-excess-risk constant (or carrying $L_p'$ explicitly) closes the conversion. The symbol is now a defined quantity with a finite bound, not a symbol-index name. $\square$

## Item 3 — the probability schedule, corrected

**Correction MR-5 (schedule condition for $\Gamma_N\to0$). [proved]**
Recall $\Gamma_N=C_0(\bar L+\tfrac\mu2)\sqrt{\dfrac{D_N\log(\Lambda N)+\log(1/\delta_N)}{N}}$ (AL-11). Replace the schedule requirement "$\delta_N\to0$" by
$$\boxed{\ \frac{\log(1/\delta_N)}{N}\ \longrightarrow\ 0\ }\qquad(\text{together with the existing }\tfrac{D_N\log(\Lambda N)}{N}\to0).$$
Then both terms under the root vanish and $\Gamma_N\to0$. The audit's counterexample $\delta_N=e^{-N}$ (which has $\delta_N\to0$ but $\log(1/\delta_N)/N=1$) is excluded, because it violates $\log(1/\delta_N)/N\to0$. Note $\log(1/\delta_N)/N\to0$ also implies $\delta_N\to0$ (indeed $\delta_N\ge e^{-o(N)}$), so this is the correct, strictly stronger schedule condition. A concrete admissible choice: $\delta_N=1/N$ (then $\log(1/\delta_N)/N=\log N/N\to0$), or any $\delta_N\ge e^{-o(N)}$. $\square$

**Declaration MR-6 (almost-sure clause — both conditions required). [declared / proved]**
The high-probability consistency bound (AL-12) holds with probability $\ge1-\delta_N$ under MR-5. For an **almost-sure eventual** statement $\limsup_N\|d_{\mathbb M}(F_{\hat\omega_N},g^\star_\mu)\|_{L^2(\mu_\zeta)}\le2h$, require **both**:
$$\text{(i) }\ \frac{\log(1/\delta_N)}{N}\to0\ \ (\text{so }\Gamma_N\to0,\ \text{the }\Phi\text{-argument vanishes}),\qquad\text{(ii) }\ \sum_N\delta_N<\infty\ \ (\text{Borel–Cantelli}).$$
Summability alone is insufficient (the audit's point: $\delta_N=e^{-N}$ is summable yet violates (i), so $\Gamma_N\not\to0$); the vanishing-generalization condition (i) is independently necessary. Both are declared explicitly; neither is folded into the other. A choice meeting both: $\delta_N=N^{-2}$ ($\log(1/\delta_N)/N=2\log N/N\to0$ and $\sum N^{-2}<\infty$). $\square$
