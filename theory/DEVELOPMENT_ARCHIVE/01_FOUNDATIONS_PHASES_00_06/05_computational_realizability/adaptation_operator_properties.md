# Properties of the Adaptation Operator U — and What Threads Through It

> **Status:** Phase-5 consolidation, 2026-08-02. Characterizes the mathematical properties any approximation of $U$ (and, where the property is joint, of $U$ and $R$ together) must satisfy. **Nothing here is an architectural component**; every item is a theorem-backed extensional constraint with its certificate. The plan audit required two additions the Phase-4 U-file deferred: query dependence and joint equivariance thread through $U$ — they cannot be delegated wholly to $R$.

---

## AP-1 Permutation symmetry
Theorem at the composite level (OP-7(i)); WLOG at the factor level (symmetrization) — an order-dependent $U$ with order-discarding $R$ realizes the same invariant composite; invariance of the factor is forced only for minimal-width factors after quotienting. *(Refereed scoping, Phase 4.)*

## AP-2 Support-location dependence
Irreducible: the task state is a function of the *labeled* sample. Witness: $V=\operatorname{span}\{t\}$, $y=1$ at $t=1$ vs $t=2$ — different predictions at every query $x\ne0$. Value-only pooling is provably wrong. *(IB-2.)*

## AP-3 Query dependence (joint with $R$ — plan-audit addition)
Three certified couplings that constrain what $U$ must *make available* to $R$: the validity region $\{x:\phi(x)\in\operatorname{row}(G)\}$ is support-dependent (so $U$'s output must determine it); the sensitivity profile is the $x$-dependent $w(x)=(G^+)^\top\phi(x)$ on the validity region; the certificate is query-dependent, and **data-dependent outside the surjective-trace stratum** — so $z_S$ must carry data-dependent certificate content in general (DM-5, strata (ii)–(iv)). A $U$ that emits only a point summary starves $R$ of the certificate channel. *(IB-3, MP-6, DM-5.)*

## AP-4 Equivariance (joint with $R$ — plan-audit addition)
The affine action on values ($\alpha\ne0$ including reflections, $\beta$; $\varepsilon\mapsto|\alpha|\varepsilon$) must commute with the *composite*; since it acts on $S$, on $z_S$'s content, and on the outputs, $U$ and $R$ must be **jointly** equivariant — the constraint threads through the factorization and cannot be localized to the readout. Derived tension: equivariant representations may require strictly more width (DM-2(c)). *(A2/IB-5.)*

## AP-5 $\varepsilon$-awareness (stratified; with the two refereed qualifications)
$\varepsilon$ enters $U$ in general. Deferral to the readout is valid exactly in the declared strata: (i) surjective-trace linear ($\operatorname{rank}G=k$): $z_S=G^+\tilde y$ $\varepsilon$-free, radius $\varepsilon\Lambda_*(x)$ data-free; (iii) Lipschitz: envelope average $\varepsilon$-free, radius $\varepsilon$-affine — *qualified by the $\varepsilon$-dependent partiality domain* (the realizability threshold moves with $\varepsilon$). In stratum (ii) (overdetermined linear — typical at $k=5$, $d\le2$) the *center* is $\varepsilon$-dependent and $G^+\tilde y$ is insufficient. In general, what fails is *canonical* deferral (the $\varepsilon=0$ envelopes as state); non-canonical $\varepsilon$-free factorizations may exist (distance-to-trace statistics). The strata are definable cells of the joint decomposition (CR-1), so stratum membership is a flagged, piecewise object. *(DM-5 with both refereed qualifications.)*

## AP-6 Uncertainty propagation
$U$ propagates enclosures, not points: upstream error $h$ transfers *inside the modulus* ($\tfrac12\omega(2\varepsilon+2h)+h$), erodes the selector on the margin collar $\{\mu<2h\}$ (CR-4 — $U$ must support the margin computation or the collar flag), and is priced outward by $\eta$-inflation (CR-6) — never absorbed silently. On empty sections $U$ emits the misspecification flag, not a state. *(IB-7/IB-8, CR-4, CR-6.)*

## AP-7 Capacity limits (with the stratified width theory — plan-audit addition)
Ceilings: at most $k$ continuous dimensions of task identity (F20/CP-3; P2); the size-$(k{+}1)$ window truncation (DM-3; NP-1). Width of $z_S$: bounded below by the query-relative rank $r_Q$ (proven, invariance of domain), **not equal to it in general** — the refereed refutation stands: the Lipschitz witness needs width $\Theta(k)$ (all non-dominated support points) though per-query rank is $O(1)$; exact width $r=\operatorname{rank}G$ only in stratum (i). Width claims must be stratified, never global. *(DM-4.)*

## AP-8 Regularity obligations
No imposed continuity across certified transitions (MP-4/CR-3: soft selection does not escape); no monotone-update structure on the *center* (OP-8) — while radii are *required* monotone ($\varepsilon$- and support-monotonicity, A4/A5); selection coherence: reported center $1$-Lipschitz in reported envelopes (IB-10); reproduction on exact identifiable data (IB-12); sensitivity-sum on constants-including classes (IB-11).

---

## Degenerate cases (mandated; consolidated)

| Case | Required behavior of $U$ | Certificate |
|---|---|---|
| $k=0$ | reduce to the baseline object (MP-2 where prior sections bounded; gauge-anchor caveat MP-3 where not) | MP-2/3 |
| query in support | no certainty under noise: floor $\varepsilon$ | F4 |
| empty section | misspecification flag, not a state | IB-8 |
| DM-5 stratum boundary | flagged stratum transition (definable cell boundary) | CR-1 |
| forced-zero archive branch | the single off-coverage exception; vanishes at $\varepsilon>0$ | C13/F18 |
| rigid classes ($N\le d$, $d=0$) | $z$ trivial; $U$ pure selection against a point class | DM-2 degenerate |

**Closing statement.** These properties characterize $U$ extensionally and exhaust what the mathematics requires of it. Anything satisfying AP-1…AP-8 *is* an admissible adaptation map, whatever its internal form; anything violating one of them is certified wrong before any evaluation is run.
