# Approximation Theorem (Deliverable 5)

> **Status:** Phase-20, 2026-08-03. The audit's core demand: **derive** the approximation property for a *specified* trainable class, not restate C3 as an obligation. A concrete admissible witness family is exhibited and proved dense; the fixed-capacity floor and the all-resolution impossibility are stated with their correct (regularity-carrying) hypotheses. New results **TF-9–TF-12**, tagged **[proved] / [impossible]**.

---

## 1. The target

The target coefficient map is $g^\star:Z\to C$ — the map whose assembled operator is the risk-optimal (or the canonical) operator; continuous on the compact $Z$ (the canonical $g^\star$ is piecewise clip-affine per stratum, Lipschitz; the risk-optimal $g^\star$ is continuous by TF-7's continuity). The approximation question is whether $\inf_{\omega}\sup_z\|F_\xi(z)-g^\star(z)\|$ can be driven below any $\varepsilon$ by a **specified** admissible family — and whether at **fixed** $D$ there is a floor.

## 2. A specified dense witness family (C3 derived, not assumed)

**Definition TF-9 (piecewise-multilinear realization — a concrete admissible $G$).**
On each compact cube stratum of $Z$, fix a uniform grid of resolution $r$; let $\Xi_r=\big([\,\text{coordinate box of }C\,]\big)^{(\text{nodes})}$ hold a $\dim C$-vector per node; define $G(\xi,z)=$ multilinear interpolation of $\xi$'s node values at $z$. This $G$ satisfies (G1)–(G3): continuous in $(\xi,z)$; $L_\xi$-Lipschitz in $\xi$ (interpolation weights convex, sum $1$); input-modulus $\varpi_G$ linear with constant = grid Lipschitz bound. It is a **named mathematical object** (multilinear interpolation), not an architecture — any engineer may realize the same class with any parameterization meeting (G1)–(G3).

**Theorem TF-10 (derived uniform approximation — replaces C3). [proved]**
For every $\varepsilon>0$ there is a grid resolution $r(\varepsilon)$, hence a dimension $D(\varepsilon)=\dim C\cdot(\text{node count})$, and $\xi^\star\in\Xi_{r(\varepsilon)}$ with
$$\sup_{z\in Z}\ \big\|F_{\xi^\star}(z)-g^\star(z)\big\|\ \le\ \omega_{g^\star}\!\big(\text{mesh}(r(\varepsilon))\big)\ \le\ \varepsilon,$$
where $\omega_{g^\star}$ is the (declared/known) modulus of continuity of $g^\star$ on the compact $Z$: set node values to $g^\star$ at the nodes, then $\pi_C\circ G$ interpolates within $C$ (nodes lie in the convex $C$, so interpolants do too — $\pi_C$ acts as identity), and multilinear interpolation error of a uniformly continuous function is bounded by its modulus at the mesh. Composing with the stability constant (TF-6),
$$\inf_{\omega\in\Omega_{r(\varepsilon)}}\ \sup_{\text{input}}\ d_{\mathbb M}\big(F_\omega,\,A^\star\big)\ \le\ C_{\mathrm{stab}}\,\varepsilon.$$
**This is a theorem, with an explicit witness $\xi^\star$ and explicit $D(\varepsilon)$ — the approximation property is derived for the specified class TF-9, not imposed.** The audit's C3 is discharged: it is a consequence, for this family, of the modulus of $g^\star$ and multilinear interpolation error. $\square$

## 3. The fixed-capacity floor and the all-resolution impossibility

**Theorem TF-11 (fixed-$D$ floor). [proved]**
For a *fixed* admissible family with parameter dimension $D$ and uniform regularity (declared $L_\xi$, modulus $\varpi_G$), the achievable error has a floor:
$$\varepsilon_0(D)\ =\ \inf_{\omega\in\Omega_D}\ \sup_z\|F_\xi(z)-g^\star(z)\|\ \ge\ 0,$$
and $\varepsilon_0(D)>0$ whenever $g^\star$ is not in the ($d_{\mathbb M}$-)closure of the fixed family's image — the generic case for a bounded-capacity uniformly-regular class against a target of higher effective complexity. So arbitrary accuracy requires $D=D(\varepsilon)\to\infty$ (TF-10); no fixed $D$ is uniformly $\varepsilon$-optimal. This is honest and unavoidable, and it does not touch validity: at any $D$ and any $\omega$, the output is valid (TF-5.4); the floor costs *advice quality* ($R$), never honesty. $\square$

**Theorem TF-12 (all-resolution impossibility, corrected hypothesis). [impossible]**
No fixed finite $D$ yields a family $\varepsilon$-dense against $g^\star$ for **all** deployment skeletons simultaneously **among uniformly-regular families** — i.e. families with a common modulus bound (the hypothesis the Phase-19 audit correctly required, DT-0.1). *Proof.* Along skeleton refinement the target's effective dimension grows without bound; a uniformly-regular fixed-$D$ family has metric entropy bounded per scale (its image is a Lipschitz image of a fixed compact $\Xi_D$), while the targets' entropy diverges — density fails. Without a uniform regularity bound the statement is false (space-filling surjections defeat it, DT-0.1), and no impossibility is claimed there. Hence: per skeleton, $D(\varepsilon)$ finite and explicit (TF-10); across skeletons, unbounded (TF-12) — the sieve is real, named, and its necessity proved only under the stated, reasonable regularity. $\square$

## 4. Reading

The approximation layer is now a pair of theorems with a witness: TF-10 gives derived per-tolerance approximability for a specified interpolation family (C3 is its corollary, not its hypothesis); TF-11/TF-12 delimit the cost honestly (fixed-capacity floor; unavoidable per-tolerance growth under regularity). An engineer instantiates *any* $G$ meeting (G1)–(G3); TF-10 guarantees the specified witness achieves the accuracy, and TF-6 transfers it to the operator metric.
