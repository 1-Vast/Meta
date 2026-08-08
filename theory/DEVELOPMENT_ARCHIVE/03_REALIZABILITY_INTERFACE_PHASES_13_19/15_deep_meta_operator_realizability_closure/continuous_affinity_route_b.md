# Continuous Affinity Interface — Route-B Instantiation (§4)

> **Status:** Phase-15, 2026-08-03. The final application requires continuous scalar outputs (affinity regression); the audits correctly held that Route A (finite outcomes, TV/Hausdorff) does not cover them, and that Route B was only a declared schema. This file **instantiates Route B**: a continuous output space, an operator metric, a proved stability condition, and the verification that the derived approximation theorem of this phase applies. New results **DM-10–DM-14**, tagged **[proved] / [declared]**. No application vocabulary beyond "continuous scalar output"; no architecture.

---

## 1. The continuous value class

**Definition DM-10 (output space, law space, operator metric). [declared structure; properties proved]**
- **Value space:** $V=[a_{\min},a_{\max}]\subset\mathbb R$, compact, bounds declared (the boundedness declaration is part of the closure class and is echoed); write $D_V=a_{\max}-a_{\min}$.
- **Law space:** $\Delta(V)$ with the **Wasserstein-1 metric** $W_1$ (equivalently Kantorovich–Rubinstein; on $V\subset\mathbb R$, $W_1(P,Q)=\int_V|F_P-F_Q|\,dv$). $(\Delta(V),W_1)$ is a compact metric space (weak topology on a compact base) — the completeness/compactness that TV *fails* to provide for continuous laws is restored by typing the metric to the weak topology.
- **Operator values:** nonempty closed convex subsets of $\Delta(V)$ with the Hausdorff-$W_1$ distance — a complete metric hyperspace; confidence and rung coordinates unchanged (PM-1's weighted metric carries over with $d_K=d_H^{W_1}$).
- **Event atlas:** threshold events $E_t=\{v\le t\}$ at a declared finite grid $T=\{t_1<\dots<t_G\}\subset V$ with mesh $h=\max_j(t_{j+1}-t_j)$ (endpoints of $V$ included). Values are **CDF-band classes**: $K=\{P:\ l_j\le F_P(t_j)\le u_j\ \forall j\}$ with monotone consistent bands (the canonical construction produces monotone envelopes).

**Why $W_1$ is the decision-correct metric [proved].** For any loss $\ell(a,\cdot)$ that is $L_\ell$-Lipschitz in the value: $|\mathbb E_P\ell-\mathbb E_Q\ell|\le L_\ell\,W_1(P,Q)$ (Kantorovich–Rubinstein duality). So risks, risk brackets, and $\Gamma$-minimax computations are $W_1$-Lipschitz — the metric measures exactly what decisions consume. *Scope:* discontinuous exact-match-type losses on a continuum are not $W_1$-controlled — and are excluded by declaration (they are also statistically meaningless on continuous outcomes); the loss class of this interface is the declared Lipschitz class.

## 2. The stability condition, proved (Route B's missing verification)

**Theorem DM-11 (CDF-band stability — the Route-B inequality with an explicit constant). [proved]**
Let $K_1,K_2$ be nonempty CDF-band classes on the same grid with bands differing by at most $\varepsilon$ (endpoint deviation, i.e. $d_{\mathrm{desc}}\le\varepsilon$). Then
$$d_H^{W_1}(K_1,K_2)\ \le\ (\varepsilon+2h)\,D_V.$$
*Proof.* Take $P\in K_1$. Define $Q$ by clamping: at each grid point set $F_Q(t_j)=\mathrm{med}\big(F_P(t_j),\,l^{(2)}_j,\,u^{(2)}_j\big)$ — within $K_2$'s bands, and $|F_Q(t_j)-F_P(t_j)|\le\varepsilon$ (the bands moved by $\le\varepsilon$); interpolate monotonically between grid points (feasible: clamped values inherit monotonicity from monotone bands and $F_P$). Then $Q\in K_2$, and $W_1(P,Q)=\int|F_P-F_Q|$: on each cell of width $\le h$, $|F_P-F_Q|\le\varepsilon+\Delta F_P(\text{cell})+\Delta F_Q(\text{cell})$ (grid agreement to $\varepsilon$ plus within-cell oscillations); integrating, $W_1\le\big(\varepsilon\cdot1+h\cdot(1+1)\big)\,D_V$ after rescaling $V$ to unit length — i.e. $\le(\varepsilon+2h)D_V$. Symmetrize. $\square$
The constant is **uniform** — no Hoffman machinery, no outcome-cardinality bound: the grid mesh $h$ replaces the bounded-outcome hypothesis, and its cost is the explicit, declared resolution term $2hD_V$.

**Proposition DM-12 (the audit's counterexample, defused by the metric typing). [proved]**
The rational-threshold counterexample (Phase-14 audit) showed that under **per-query TV typing**, countably many threshold queries force an infinite statistic. Under the $W_1$ typing, monotonicity closes it: the finite grid statistic determines *every* threshold event probability to within the adjacent band values plus cell mass, and any two laws agreeing to $\varepsilon$ at the grid are within $(\varepsilon+2h)D_V$ in $W_1$ (DM-11's integral estimate). A **continuum of threshold queries is controlled by a finite statistic with explicit mesh error** — the infinite-dimensionality was an artifact of demanding per-query exactness in the strong metric, not a fact about the information. (FIN-ATLAS) in Route B is therefore not a scope loss but a resolution choice, priced by $h$. $\square$

## 3. The chain, verified in Route B

**Theorem DM-13 (statistical layer). [conditional on the declared stack]**
Per-task identified objects for a scalar query are intervals $[l_i^{\mathrm{id}},u_i^{\mathrm{id}}]$ (frozen scalar sections). Forced/compatible indicators at threshold $t$: forced $\iff u_i^{\mathrm{id}}\le t$; compatible $\iff l_i^{\mathrm{id}}\le t$ — two monotone one-dimensional threshold classes of VC dimension $1$ each; over the finite grid and contexts, the declared class is finite-VC with $d^*=O(1)$ (plus the context term). The uniform GC theorem, the missing-fiber accounting, the confidence schedule — PM-4/PM-5 — apply verbatim with endpoints now CDF-band values. Support restriction: laws supported on the current identified interval form a closed convex subset in $W_1$; restriction remains likelihood-free. $\square$

**Theorem DM-14 (the derived approximation theorem applies). [proved under (FIN-ATLAS: the grid)]**
The Phase-15 construction transfers verbatim: the statistic $E(H)$ (stratum label + per-fiber grid frequencies) is finite-dimensional on the same compact stratified domain (DM-1's form, with $q=2G|C_\kappa|$); the canonical endpoint map is again per-stratum $\mathrm{clip}\circ(\text{affine})$, exactly representable on horizon strata (DM-5) and $\varepsilon$-matched on tail strata (DM-6), with confidence/rung exact by shared postprocessing. The endpoint-to-operator transfer uses DM-11 in place of Hoffman:
$$\sup_H\ d_{\mathbb M}\big(A_{\theta^\star}(H),A_\phi(H)\big)\ \le\ \alpha\,D_V\,\big(\varepsilon+2h\big),$$
with finite explicit $p(\varepsilon)$ as in DM-6 (with $G$ now the declared value grid). **The same approximation theorem holds; the only Route-B novelties are the stability constant $D_V$ and the irreducible, declared resolution term $2hD_V$** — which is honest: a continuous law class represented at grid resolution $h$ cannot promise better than its own mesh, and the certificate says so. Total-error decomposition DM-9 then holds for the continuous interface with these substitutions. $\square$
