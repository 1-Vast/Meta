# Learning Objective (§4)

> **Status:** Phase-18, 2026-08-03. The objective, its spaces, its metric, its generalization target, and the calibration semantics — with the certificate firewall restated where it binds. New results **MI-10–MI-12**, tagged **[proved] / [conditional] / [declared]**.

---

## 1. The objective

$$\theta^\star\ \in\ \operatorname*{arg\,min}_{\theta\in\Theta}\ \ R(\theta)\ =\ \mathbb E_{T\sim p(T)}\Big[\ \mathrm{operator\_loss}\big(A_\theta(S_T),\ Q_T\big)\ \Big],\qquad \Theta=[0,1]\times\mathbb B^m\ \text{(compact convex)},$$
with $p(T)=\Pi_{\mathrm{obs}}$ the observable task law and $Q_T$ the task's **identified** query information (value / forced-compatible interval / admissible-order set — always observable, never a latent mark).

- **Loss space [declared class; properties proved]:** convex, Lipschitz functionals $L:\mathbb B\times\mathcal Y\to[0,\infty)$ scoring the emitted description against identified query information; canonical instance the interval/band score (width $+$ $\tfrac2\alpha\cdot$violation distance, summed over declared events/grid; censored outcomes scored against their compatible regions) — convex and Lipschitz on the compact $\mathbb B$ (MR-11).
- **Metric [carried]:** outputs live in the weighted operator value metric (Hausdorff-$d_K$ on classes $+$ confidence $+$ rung; PM-1 / MR-8); the loss is $d_K$-compatible: band-Lipschitz losses are class-Lipschitz through the stability constants (Hoffman on Route A; $\varepsilon D_V+2h$ on Route B), so risk differences are controlled by operator-metric differences — the loss and the metric tell one story.

## 2. Existence and generalization

**Theorem MI-10 (existence; convex program). [proved — carried MR-12]** $R$ is convex and Lipschitz on compact convex $\Theta$; a global minimizer exists; the minimizing set is a convex face; parameter gauge on it is irrelevant.

**Theorem MI-11 (generalization target and finite-task bound). [conditional on task-(IID)/(C-IID-$\kappa$)]**
Let $\hat\theta_N$ minimize the empirical risk over $N$ tasks. Since $\Theta$ is a compact subset of $\mathbb R^p$ (fixed $p$) and the task-wise loss is bounded by $\bar L$ and $\mathrm{Lip}_\theta$-Lipschitz in $\theta$ uniformly in $T$ (affine decoder $\times$ Lipschitz loss), a covering-number argument over $\Theta$ gives, with probability $\ge1-\delta$:
$$R(\hat\theta_N)\ \le\ \min_\Theta R\ +\ C\,\bar L\,\sqrt{\frac{p\,\ln\big(\mathrm{Lip}_\theta N\big)+\ln(1/\delta)}{N}},$$
— the **generalization target is the population operator risk $R$**, and the rate is dimension-$p$ classical; under (C-IID-$\kappa$) the statement holds fiber-wise with fiber counts, inheriting the missing-fiber accounting (PM-4/5). Under undeclared shift, DE-T3's adversarial reversal applies to the *risk* claims exactly as to everything else — tagged, never silent. $\square$

## 3. Uncertainty calibration

**MI-12 (calibration semantics — elicitation, displayed; certificates, firewalled). [proved / declared]**
(i) *[proved — classical elicitability]* The interval score at level $\alpha$ elicits central quantile bands: at the population optimum, the learned band's conditional coverage of the identified query outcome equals $1-\alpha$ per context fiber (first-order condition). **Calibration is therefore a well-defined, testable population property of the learned component**: empirical fiber-wise coverage vs $1-\alpha$, with its own confidence tags — a diagnostic of the *preference* layer.
(ii) *[declared firewall — carried MR-13]* Calibration statements are rung-tagged population claims and may never be emitted as worst-case certificates; the certificate rows ($b_{\mathrm{can}}$-derived floors, envelopes, flags) are $\theta$-invariant and untouched by training. A perfectly calibrated learned band and the outer identified band answer different questions — the ledger keeps both, per the standing Phase-7 separation. Miscalibration, correspondingly, is a *performance* defect (priced by $R$ and visible in the coverage diagnostic), never a validity defect: no $\theta$ can emit an invalid object (MR-9). $\square$
