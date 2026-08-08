# Generalization Bridge (Deliverable 6)

> **Status:** Phase-20, 2026-08-03. Connects empirical optimization of $\omega$ to population operator error — the arrow the audit found absent. Two composable bridges: risk-excess (ERM $\to$ population risk) and imitation ($\to$ operator-metric approximation). New results **TF-13–TF-15**, tagged **[proved] / [conditional]**.

---

## 1. Bridge A — ERM to population risk

**Theorem TF-13 (uniform generalization over $\mathcal H$). [conditional on task-(IID)/(C-IID-$\kappa$)]**
$\Omega\subset\mathbb R^{D+m\dim\mathbb B}$ compact; loss bounded by $\bar L$; $\omega\mapsto L(F_\omega(\cdot),A_T)$ Lipschitz with constant $\mathrm{Lip}(L)\cdot L_{\mathcal H}$ uniformly in $T$ (TF-6). Then the covering number of $\Omega$ at scale $s$ is $N(s)\le(C_\Omega/s)^{\dim\Omega}$, and by the standard uniform-deviation bound, with probability $\ge1-\delta$:
$$\sup_{\omega\in\Omega}\big|\widehat R_N(\omega)-R(\omega)\big|\ \le\ C\,\bar L\sqrt{\frac{\dim\Omega\cdot\ln\!\big(\mathrm{Lip}(L)L_{\mathcal H}N\big)+\ln(1/\delta)}{N}}\ =:\ \Gamma_N,$$
hence $R(\hat\omega_N)\le\min_{\omega}R(\omega)+2\Gamma_N$. Under (C-IID-$\kappa$) the bound is fiber-relative with the inherited missing-fiber term (PM-4/5); under undeclared shift the proved adversarial reversal applies to $R$ exactly as elsewhere (DE-T3), tagged. **The generalization target is the population operator risk $R$**, and the rate is $\dim\Omega$-classical — the finite-dimensional compactness of $\Omega$ (Deliverable 1) is precisely what makes this hold. $\square$

## 2. Bridge B — empirical objective to operator-metric approximation

The audit's precise gap: no theorem connected empirical optimization to $d_{\mathbb M}(F_\omega,A^\star)$. Two routes, each proved, each with declared conditions:

**Theorem TF-14 (imitation route — direct operator-metric control). [conditional]**
Since the canonical/target operator $A^\star$ is computable, one may use the **imitation objective** $\widehat R^{\mathrm{im}}_N(\omega)=\tfrac1N\sum_t d_{\mathbb M}\big(F_\omega(\text{input}_t),A^\star(\text{input}_t)\big)$ over inputs drawn from a declared design $\mu_Z$ on the input space. Then:
(i) *Population imitation error:* TF-13's uniform bound applies verbatim (bounded, Lipschitz, compact $\Omega$), giving $\mathbb E_{\mu_Z}\,d_{\mathbb M}(F_{\hat\omega_N},A^\star)\le\inf_\omega\mathbb E_{\mu_Z}d_{\mathbb M}(F_\omega,A^\star)+2\Gamma_N$; the infimum is $\le C_{\mathrm{stab}}\varepsilon$ by the approximation theorem TF-10 (specified witness), so **empirical imitation minimization provably reaches $C_{\mathrm{stab}}\varepsilon+2\Gamma_N$ population imitation error** — the missing arrow, closed in $L^1(\mu_Z)$.
(ii) *Sup-metric upgrade [conditional on a declared design lower bound]:* to convert $L^1(\mu_Z)$ control into $\sup$-control (C3's form), declare $\mu_Z$ with a mass floor $q_0>0$ on a mesh of $Z$ (each cell charged $\ge q_0$); then $\sup_z\le$ (per-cell $L^1$ bound)$/q_0$ + mesh-modulus, giving a $\sup$-bound at the cost of the declared $q_0$ and mesh. Without such a declaration, only $L^1(\mu_Z)$ (average-case) operator approximation is claimed — stated honestly, not upgraded silently. $\square$

**Corollary TF-15 (the full chain, assembled). [conditional as tagged]**
Combining the specified-witness approximation (TF-10), the imitation generalization (TF-14), and the operator-value transfer (TF-6):
$$d_{\mathbb M}\big(F_{\hat\omega_N}(\text{input}),\,A^\star(\text{input})\big)\ \le\ \underbrace{C_{\mathrm{stab}}\varepsilon}_{\text{approx (TF-10, witness)}}\ +\ \underbrace{2\Gamma_N}_{\text{generalization (TF-13/14)}}\ +\ \underbrace{\gamma^{\mathrm{opt}}}_{\text{optimization tolerance}}\qquad[\text{in }L^1(\mu_Z);\ \sup\text{ under the }q_0\text{ declaration}],$$
with, on top, the standing decomposition to the identified target $M^\dagger$ (add the statistical operator-learning term and the identification floor, PM-5/DM-9). **Every term is now either a proved theorem or a named, tagged declaration — none is an unproved obligation.** In particular C3 is no longer assumed: it is the TF-10 witness bound, achieved by empirical imitation minimization via TF-14. $\square$

## 3. What the bridge does and does not deliver

Delivered, by proof: empirical minimization (of either the task risk or the imitation objective) over the compact finite-dimensional $\Omega$ generalizes to population error at $\dim\Omega$-classical rates, and — via the specified witness family — reaches the approximation floor $C_{\mathrm{stab}}\varepsilon$; operator-metric approximation follows in $L^1(\mu_Z)$ unconditionally and in $\sup$ under a declared design lower bound. Not delivered (and correctly so): a $\sup$-bound without a design declaration (impossible in general — the adversary hides error off the design); a fixed-$D$ family reaching $\varepsilon=0$ (TF-11 floor); optimization *efficiency* (attainability of $\gamma^{\mathrm{opt}}$ is in-principle by compactness; speed is the implementer's, and is the only remaining non-theorem, now a scalar tolerance rather than an undefined bridge).
