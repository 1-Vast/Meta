# Parameterization of the Adaptation Operator

> **Status:** Phase-4 derivation, 2026-08-02. Question: is the three-stage decomposition
> $$z=\Phi(\mathbb T),\qquad z_S=U(z,S,\varepsilon),\qquad (\hat y,\hat\rho)=R(z_S,x)$$
> — shared prior state / task adaptation / query-conditional readout — mathematically justified? All claims refereed; the $\varepsilon$-placement and width claims below are the corrected, stratified versions (the referee refuted the unstratified originals with explicit counterexamples).

---

## 1. Verdict up front

The decomposition is **justified as exact (non-lossy), forced in structure, but not unique and not minimal in general**:

- **Existence/exactness:** always, when $\mathcal X$ is finite-dimensional — the support set itself (locations *and* values, with $\varepsilon$) is a finite-dimensional $x$-independent intermediate, so nothing is lost by the pipeline shape (DM-4(a)).
- **$\Phi$ exists by minimal sufficiency** (CP-1/CP-2): the family enters the problem only through the window system; $\Phi$ is the representation of `latent_operator_theory.md`, with all its dimension floors, gauge, truncation, and outer-admissibility constraints.
- **Structure is forced** in the specific senses of §3 — but several familiar-looking "must"s are only *without-loss-of-generality*s, and are labeled as such.
- **Minimality is regime-dependent** (§4): there is no single formula for the minimal intermediate width; only regime-stratified bounds.

---

## 2. Existence and non-lossiness (DM-4(a))

Let $\mathcal X\subseteq\mathbb R^p$. Then $z_S=(x_1,\dots,x_k,\tilde y_1,\dots,\tilde y_k,\varepsilon)\in\mathbb R^{k(p+1)+1}$ is an $x$-independent finite-dimensional intermediate from which, together with $z=\Phi(\mathbb T)$, the exact section — hence center and radius — at every query is reconstructed. The pipeline is therefore *never* the bottleneck; the mathematical content lies in (i) how small $z_S$ can be, (ii) where $\varepsilon$ must enter, (iii) what each stage is forced to preserve.

---

## 3. What is forced, and what is only WLOG

1. **Permutation symmetry: composite-level theorem, factor-level WLOG (refereed correction).** A1 constrains the *composite* operator; a factorization with order-dependent $U$ and order-discarding $R$ realizes the same invariant composite. The provable statement: an invariant factorization always exists (symmetrize $U$); invariance of the factor itself is forced only for minimal-width factors after quotienting. "U must be a set function" overstates; "U may be taken to be one, at no cost" is the theorem.
2. **Support-location dependence is irreducible.** $V=\operatorname{span}\{t\}$, $k=1$: the value $y=1$ observed at $t=1$ versus $t=2$ gives $c=1$ versus $c=\tfrac12$ — different predictions at every query $x\ne0$ (refereed wording). Any intermediate discarding locations is wrong; "task state" is a function of the *labeled* sample, not the value multiset.
3. **Joint equivariance (A2).** The affine action on values threads through all three stages: $U$ and $R$ must be jointly equivariant (including reflections) for the composite to be; this couples the stages and is a genuine constraint on any factorization, with the equivariance-vs-width tension recorded in DM-2(c).
4. **The certificate channel flows through all stages.** $R$ outputs the pair $(\hat y,\hat\rho)$ with $\hat\rho\in[0,+\infty]$ (compactified codomain — the radius genuinely attains $+\infty$ at validity boundaries); $z_S$ must carry whatever data-dependent certificate content the regime demands (§5).
5. **No continuity may be imposed on $R\circ U$ in the data** (MP-4), and **no monotone-update structure** (OP-8): both would contradict the optimal operator.

---

## 4. Width of the intermediate: stratified bounds, no equality (DM-4(b,c); refereed correction)

The original "minimal width $=$ query-relative rank" is **false as an equality**; what survives:

- **Lower bound (proven).** Any continuous $x$-independent intermediate from which a readout recovers all identifiable-query values of a query set $Q$ has $\dim z_S\ge r_Q=\dim\operatorname{span}\{\phi(x):x\in Q\ \text{identifiable}\}$ — the data-to-values map has an $r_Q$-dimensional linear image through which $z_S$ must factor injectively (invariance of domain; the slice argument was refereed).
- **No matching upper bound in general.** The intermediate must serve *all* queries at once, so it is governed by the summary sandwich (F19), not by per-query ranks. **Lipschitz witness:** each query's envelopes depend on $\sim2$ active support points (per-query rank $O(1)$), but different queries activate different points, and the $x$-independent intermediate must determine both envelopes at every $x$ — forcing width $\Theta(k)$ (all non-dominated support points). Per-query economy does not compose into global economy.
- **Exact width in the surjective-trace linear regime** ($\operatorname{rank}G=k$; see §5): $z_S=G^+\tilde y$ of dimension $r=\operatorname{rank}G$ suffices for all identifiable queries, with readout $\phi(x)^\top z_S$ — here lower and upper bounds meet.

---

## 5. Where $\varepsilon$ must enter: the stratified deferral theorem (DM-5; refereed and corrected)

Write "deferral" for the property that $U$ may be $\varepsilon$-free with $\varepsilon$ entering only the readout. The truth is a four-way stratification — the referee refuted the unstratified claim with numerical counterexamples:

| Regime | Center | Radius | Deferral? |
|---|---|---|---|
| **(i) Linear, surjective trace** ($\operatorname{rank}G=k$) | $\phi(x)^\top G^+\tilde y$ — $\varepsilon$-free | $\varepsilon\Lambda_*(x)$ — data-free | **Yes**: $\varepsilon$ enters only $R$; certificate needs no data |
| **(ii) Linear, overdetermined** ($\operatorname{rank}G<k$) | $\varepsilon$-**dependent** | data- and $\varepsilon$-dependent | **No.** Counterexample ($d=1$, $f_\theta(t)=\theta t$, $D=\{1,2\}$, $\tilde y=(0.1,0)$): section $[0,0.05]$ at $\varepsilon=0.1$ but $[-0.1,0.1]$ at $\varepsilon=0.2$ — the *center* moves with $\varepsilon$; moreover $G^+\tilde y$ is not even sufficient (two data vectors with equal $G^+\tilde y$, different centers — numerically verified) |
| **(iii) 1-Lipschitz class** | envelope average — $\varepsilon$-free (the $\pm\varepsilon$ shifts cancel) | ($\varepsilon=0$ width$)+2\varepsilon$ — $\varepsilon$-affine | **Yes**, up to the $\varepsilon$-dependent *partiality domain* (the realizability threshold moves with $\varepsilon$) |
| **(iv) General** | — | — | **Canonical deferral fails**: the $\varepsilon=0$ envelopes cannot serve as the deferred state (single-member family: all off-trace data have $S_0=\emptyset$ yet different $\varepsilon$-thresholds). *Refereed scoping:* this rules out the canonical $\varepsilon$-free state, not every $\varepsilon$-free factorization (a distance-to-trace statistic restores deferral in that example) |

**Interface consequence:** since overdetermined designs ($\operatorname{rank}G<k$) are *typical* at budget $k=5$ with task dimension $d\le2$, the safe general form has $\varepsilon$ entering $U$; $\varepsilon$-late factorizations are correct only in declared stratum (i)/(iii) regimes.

**Certificate content of $z_S$ (refereed scoping):** data-free exactly in stratum (i) — the section width over realizable data is the constant $2\varepsilon\Lambda_*(x)$ (translation of the $\varepsilon$-slab; valid *only* at $\operatorname{rank}G=k$: for constants on two points the width is $2\varepsilon-|\tilde y_1-\tilde y_2|$, genuinely data-dependent). In strata (ii)–(iv) the intermediate must carry data-dependent certificate information.

---

## 6. The additive form

$\;\hat y=B(x)+\Delta(z_S,x)$ is available exactly where prior sections are bounded and data realizable (MP-2), with $B=\operatorname{cen}T_x$, $|\Delta|\le\tfrac12\operatorname{diam}T_x$, and guarantee gain $\tfrac12(\operatorname{diam}T_x-\omega(2\varepsilon))$; where prior sections are unbounded, any anchor is a gauge choice (MP-3) and the additive form is a declared convention, not a derived structure. The unconstrained linear regime — stratum (i) — has *unbounded* priors: the additive-around-baseline reading and the exact-width reading of §5(i) are **different regimes**, not to be conflated.

---

## 7. Degenerate cases

$k=0$: $U(z,\emptyset,\varepsilon)$ must reduce to the baseline object (MP-2 where bounded; the gauge caveat MP-3 where not). Empty section (unrealizable support): $U$ must emit the misspecification flag, not a state (A10). Query in support: floor $\varepsilon$ (F4) — $R$ may not return certainty at support points under noise. Validity regions of measure zero: with ambient $d>k$ in stratum (i), $\{x:\phi(x)\in\operatorname{row}(G)\}$ is typically a null set of queries — the honest $R$ outputs $\hat\rho=+\infty$ almost everywhere, and finite-radius training targets exist only on a null set (recorded as a learnability constraint in `learnability_conditions.md`).

---

## 8. Conclusion

The candidate form $z=\Phi(\mathbb T)$, $z_S=U(z,S,\varepsilon)$, $(\hat y,\hat\rho)=R(z_S,x)$ is justified: it can represent the canonical operator exactly, and the theory pins down what each stage must preserve (§3), how wide the middle must be (§4, stratified), and where $\varepsilon$ lives (§5, stratified). What the theory *refutes* is any claim that this factorization is unique, that its middle can be a fixed small dimension across regimes, that $\varepsilon$ can generally be deferred to the readout, or that the intermediate can forget support locations, the certificate channel, or the partiality flags.
