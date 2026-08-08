# Meta-Learning Contract, Fully Defined (Item 2)

> **Status:** Phase-25.1, 2026-08-03. Defines every learning symbol the consistency theorem invokes — the audit's finding 6 ($F_{\hat\omega_N},\Omega_N,\Gamma_N,L_p,\gamma^{\mathrm{opt}}_N$, plus the task sample, empirical objective, estimator, hypothesis class, and the population/excess-risk link all undefined in the freeze folder). Definitions **FC-12–FC-18**, theorem restatement **FC-19**. Tagged **[definition] / [proved] / [conditional]**.

---

## 1. The learning problem, defined

**FC-12 [definition] Task distribution and sample.** $P_T$ = the observable task law (the fixed deployment's population; A-STAT: tasks (IID), or (C-IID-$\kappa$) with fiber counts, from $P_T$; within-task noise is the frozen adversarial-bounded model). A **task** is $T_i=(S_i,Q_i)$ with observable support $S_i$ and query $Q_i$; its identified target information is $A_i$ (point-valued affinity on the point-supervised channel). The **meta-training sample** is $(T_1,\dots,T_N)$. Each task induces a statistic $\zeta_i=z(S_i,Q_i,\gamma)$ and a per-task loss $L(B(\zeta_i)\,F(\zeta_i),A_i)$.

**FC-13 [definition] Parameter space and hypothesis class at level $N$.** $\Omega_N=\Xi_N\times\Delta_m^{\,?}$... precisely: $\Omega_N\subset\mathbb R^{D_N}$ a **nonempty compact** parameter set of dimension $D_N=\dim\Omega_N$ (the coefficient-map parameters; anchors are fixed constants of $\mathcal D$, so not parameters). The **realization** $G_N:\Xi_N\times Z\to\mathbb R^{m+1}$ is jointly continuous, $L_\xi$-Lipschitz in the parameter, with declared input-modulus; the **coefficient map** is $F_\omega=\pi_{\Delta_m}\circ G_N(\omega,\cdot):Z\to\Delta_m$ ($\pi_{\Delta_m}$ the $1$-Lipschitz metric projection onto the simplex — makes $F_\omega$ $\Delta_m$-valued for every $\omega$). The **hypothesis class** is
$$\mathcal H_N=\{\,F_\omega:\omega\in\Omega_N\,\},\qquad \text{operator form }F_\omega(S,Q,\gamma;z_H^0)=\big(K(B(z)F_\omega(z))|_{\mathrm{supp}\,I(S)},\ \mathrm{conf}(z_H^0),\ \mathrm{rung}(z_H^0)\big).$$

**FC-14 [definition] Empirical and population objectives.**
$$\widehat R_{\mu,N}(\omega)=\tfrac1N\sum_{i=1}^N\Big[L\big(B(\zeta_i)F_\omega(\zeta_i),A_i\big)+\tfrac\mu2\|F_\omega(\zeta_i)\|^2\Big],\qquad R_\mu(F_\omega)\ \text{as in FC-7}.$$
$\widehat R_{\mu,N}$ is the empirical regularized risk (the average whose expectation is $R_\mu(F_\omega)$); the ridge is evaluated at $F_\omega(\zeta_i)$, consistent with the target's definition.

**FC-15 [definition] Estimator.** $\hat\omega_N\in\arg\min_{\omega\in\Omega_N}\widehat R_{\mu,N}(\omega)$ (a minimizer exists — $\widehat R_{\mu,N}$ continuous on compact $\Omega_N$; a measurable selection is fixed, so $\hat\omega_N$ is a well-defined random element). $F_{\hat\omega_N}$ is the learned coefficient map.

## 2. The three error terms, defined

**FC-16 [definition] Approximation term.** $\varepsilon_{\mathrm{approx}}(N)=\inf_{\omega\in\Omega_N}\sup_{z\in Z}\|F_\omega(z)-g^\star_\mu(z)\|$ — the best sup-coefficient error achievable in $\mathcal H_N$ against the single target $g^\star_\mu$ (finite; $\to0$ as $\dim\Omega_N\to\infty$ along the schedule, by the approximation theorem, since $g^\star_\mu$ is continuous, FC-6/Phase-24 FD-5.5).

**FC-17 [definition] Generalization term.** $\Gamma_N=C\,\bar L\sqrt{\dfrac{D_N\ln N+\ln(1/\delta_N)}{N}}$ — the uniform-deviation bound over $\mathcal H_N$: $C$ an absolute constant, $\bar L$ the loss bound (FC-6), $D_N=\dim\Omega_N$, $\delta_N\in(0,1)$ the declared confidence level. It bounds $\sup_{\omega\in\Omega_N}|\widehat R_{\mu,N}(\omega)-R_\mu(F_\omega)|$ with probability $\ge1-\delta_N$ (covering-number bound; $\Omega_N$ compact, the map $\omega\mapsto$ per-task loss Lipschitz with constant $\mathrm{Lip}(L)L_{\mathcal H_N}$, absorbed into $C$/logarithms). Under (C-IID-$\kappa$) it is fiber-relative with the missing-fiber term; under undeclared shift the DE-T3 reversal applies, tagged.

**FC-18 [definition] Optimization term.** $\gamma^{\mathrm{opt}}_N\ge0$ — the declared tolerance to which $\hat\omega_N$ approximately minimizes $\widehat R_{\mu,N}$: $\widehat R_{\mu,N}(\hat\omega_N)\le\inf_{\omega\in\Omega_N}\widehat R_{\mu,N}(\omega)+\gamma^{\mathrm{opt}}_N$. ($L_p$ of FC-6 is the coefficient-loss Lipschitz constant used to convert the sup-coefficient approximation error into an excess-risk contribution: $\inf_{\mathcal H_N}R_\mu-R_\mu(g^\star_\mu)\le L_p\,\varepsilon_{\mathrm{approx}}(N)$.)

## 3. The consistency theorem, every symbol defined

**Theorem FC-19 (fixed-deployment consistency, self-contained). [conditional on the declared schedule]**
Impose the schedule: $D_N\to\infty$ with the coefficient resolution tied to $D_N$; $\dfrac{D_N\ln N}{N}\to0$; $\delta_N\to0$ with $\dfrac{\ln(1/\delta_N)}{N}\to0$; **output mesh $h$ fixed**. Decompose the excess risk of the estimator:
$$\mathcal E_\mu(F_{\hat\omega_N})=\underbrace{\big[R_\mu(F_{\hat\omega_N})-\inf_{\Omega_N}R_\mu\big]}_{\le\,2\Gamma_N+\gamma^{\mathrm{opt}}_N}+\underbrace{\big[\inf_{\Omega_N}R_\mu-R_\mu(g^\star_\mu)\big]}_{\le\,L_p\,\varepsilon_{\mathrm{approx}}(N)}\ \le\ 2\Gamma_N+\gamma^{\mathrm{opt}}_N+L_p\,\varepsilon_{\mathrm{approx}}(N),$$
(first bracket: uniform deviation FC-17 twice + empirical optimality FC-15 + tolerance FC-18; second: FC-18/FC-16). Feed into the calibration theorem FC-11: with probability $\ge1-\delta_N$,
$$\big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ \Phi\big(2\Gamma_N+\gamma^{\mathrm{opt}}_N+L_p\,\varepsilon_{\mathrm{approx}}(N)\big)+2h,$$
and since the $\Phi$-argument $\to0$ (each term $\to0$ on the schedule) and $\Phi$ is continuous at $0$ (FC-10),
$$\limsup_{N\to\infty}\ \big\|d_{\mathbb M}(F_{\hat\omega_N},\,g^\star_\mu)\big\|_{L^2(\mu_\zeta)}\ \le\ 2h\qquad(\text{a.s.-eventually, }\delta_N\to0).$$
Every symbol — $P_T$, $T_i$, $\mathcal H_N$, $\Omega_N$, $\widehat R_{\mu,N}$, $\hat\omega_N$, $R_\mu$, $\mathcal E_\mu$, $\Gamma_N$, $\gamma^{\mathrm{opt}}_N$, $L_p$, $\varepsilon_{\mathrm{approx}}$, $\Phi$, $d_{\mathbb M}$, $h$, $g^\star_\mu$ — is defined in this package (FC-1–FC-18). $\limsup\le2h$ per the FD-2 correction (no "$\to2h$"). Single target throughout; fixed deployment; no continuum, no ranking, no varying $z_H$. $\square$
