# Structure-Preserving Parameterization (§3)

> **Status:** Phase-17, 2026-08-03. The valid-description polytope $\mathbb B$, the closed continuous value class, and the proofs that **every** parameter decodes to a valid operator and that interpolation preserves every constraint. New results **MR-7–MR-10**, tagged **[proved] / [declared]**. Tool choice per mandate: **convex parameterization** (with simplex/partition weights); exponential families, kernel means, normal cones, and implicit layers are *not* needed — convexity alone suffices, and minimality forbids the rest.

---

## 1. The valid-description polytope $\mathbb B$

**Definition MR-7 (lifted linear description). [proved: compact convex; explicit]**
Under (SKEL), a band vector stacks per (context, index) the event endpoint pairs $(l(E),u(E))_{E\in\mathcal E}$ (Route A) or the grid CDF bands $(l_j,u_j)_{j\le G}$ (Route B). Define the lifted polytope in $(b,P)$-variables — $P$ a witness law, finitely represented (Route A: a point of $\Delta(\Omega)$ per index; Route B: its grid values $F_j\in[0,1]$):
- **ranges:** $0\le l\le u\le1$ per event/grid point;
- **witness feasibility:** $l(E)\le P(E)\le u(E)$ (Route A); $l_j\le F_j\le u_j$ with $F_j\le F_{j+1}$, $F_G=1$ (Route B) — all affine;
- **coherence ties (Route A):** for every declared coarsening $h$ and coarse event $E'$: $l'(E')\le l(h^{-1}E')$ and $u(h^{-1}E')\le u'(E')$ — affine; these *imply* $h_*K(b_{\mathrm{fine}})\subseteq K(b_{\mathrm{coarse}})$ for the decoded classes (for $P\in K(b_{\mathrm{fine}})$: $h_*P(E')=P(h^{-1}E')\in[l(h^{-1}E'),u(h^{-1}E')]\subseteq[l'(E'),u'(E')]$) — the audit's coherence counterexample is excluded **by constraint, for every parameter**, not by the provenance of canonical estimates;
- **band monotonicity (Route B):** $l_j\le l_{j+1}$, $u_j\le u_{j+1}$ — affine.
$\mathbb B$ = the projection onto the $b$-coordinates: a **projection polytope** — compact, convex, with an explicit finite linear description in the lifted variables (a *closed feasible representation*, the mandate's requirement); nonemptiness of every $K(b)$, $b\in\mathbb B$, holds by the carried witness, and convex combinations of $(b,P)$-points certify convexity of feasibility (MR-1). Emptiness of $\mathbb B$ itself is a declaration-consistency error checked once at deployment (LP feasibility) and flagged. $\square$

## 2. The Route-B value class, closed (the codomain repair)

**Definition/Theorem MR-8 (closed band classes in $W_1$). [proved]**
Replace the audited non-closed convention by the **closed-constraint convention**: for grid points $t_j$,
$$K(b)\ =\ \big\{\,P\in\Delta(V)\ :\ P\big([a_{\min},t_j]\big)\ \ge\ l_j\ \ \text{and}\ \ P\big([a_{\min},t_j)\big)\ \le\ u_j\ \ \forall j\,\big\}.$$
Each constraint is $W_1$- (weak-) closed by the Portmanteau theorem: for closed sets, $\limsup_n P_n(F)\le P(F)$ makes $\{P(F)\ge l\}$ closed; for open sets, $\liminf_n P_n(G)\ge P(G)$ makes $\{P(G)\le u\}$ closed. Hence $K(b)$ is a closed convex subset of the compact $(\Delta(V),W_1)$ — compact convex, a genuine element of the declared hyperspace. *Audit witness re-run:* upper band $0$ at $t$ now reads $P([a_{\min},t))\le0$; every $\delta_{t+1/n}$ satisfies it, and the $W_1$-limit $\delta_t$ **also** satisfies it ($\delta_t([a_{\min},t))=0$) — the sequence converges inside the set; closedness restored. The semantic price is the left-limit convention $F(t^-)\le u_j$, which changes no integral quantity ($F$ is altered on a null set) and no decision value for the declared Lipschitz loss class.
**Stability, re-typed [proved]:** the monotone-clamping bound is dimensionally corrected to
$$d_H^{W_1}\big(K(b_1),K(b_2)\big)\ \le\ \varepsilon\,D_V\ +\ 2h,$$
$\varepsilon=$ band sup-deviation (dimensionless), $h=$ absolute mesh in value units: $W_1=\int|\Delta F|\,dv\le(\text{grid deviation }\varepsilon)\cdot D_V+\sum_{\text{cells}}(\text{cell length})(\Delta F_P+\Delta F_Q)\le\varepsilon D_V+2h$ — the audit's typing correction adopted. $\square$

## 3. Validity and constraint-preservation of the full decoder

**Theorem MR-9 (every parameter is valid; interpolation is structure-preserving). [proved]**
With $\Theta=[0,1]\times\mathbb B^m$ and the MR-3 decoder:
(i) for every $(\theta,H)$, the decoded band is a convex combination of points of $\mathbb B$ ($b_{\mathrm{can}}(H)$ included — the canonical vector carries its own empirical witness and satisfies the ties), hence lies in $\mathbb B$; its class $K(\cdot)$ is nonempty, compact, convex, coherent, and (Route B) $W_1$-closed: **$A_\theta(H)\in\mathbb M$ always** — the audited feasibility ($l=0.9>u=0.1$), coherence, and monotonicity counterexamples are all unconstructible, because the coordinates they used independently are now jointly constrained by the type of $\Theta$;
(ii) the interpolation weights $\varphi_j(z)$ act *inside* the convex set — the operation that broke the cube preserves every constraint on $\mathbb B$ by definition of convexity; no repair projection, no implicit layer, no feasibility post-hoc step is needed (and none is used: minimality);
(iii) $(\theta,H)\mapsto A_\theta(H)$ is a Carathéodory map into $\mathbb M$ — affine (hence continuous) in $\theta$; measurable in $H$ ($b_{\mathrm{can}}$ and $z$ measurable); and now the codomain membership holds globally, so joint measurability *into $\mathbb M$* is established (the audit's objection was membership, not the map);
(iv) parameter continuity in the operator metric holds globally: bands move affinely in $\theta$, and the stability theorems (Hoffman on Route A between same-pattern nonempty polytopes — always the case here; MR-8 on Route B) transfer band motion to $d_{\mathbb M}$ with uniform constants. $\square$

**Remark MR-10 (tool minimality).** The mandate's tool list was consulted; only convexity + simplex/partition coordinates are load-bearing. Exponential families would impose smooth interiors the bands don't need; kernel means and normal cones re-encode what the lifted polytope already states linearly; implicit layers would reintroduce a feasibility computation the type system makes unnecessary. Each is rejected as not mathematically necessary. $\square$
