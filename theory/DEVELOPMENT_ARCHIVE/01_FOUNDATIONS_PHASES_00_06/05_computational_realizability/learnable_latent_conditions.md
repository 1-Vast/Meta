# Learnable Latent Conditions

> **Status:** Phase-5 derivation, 2026-08-02. Question: does a finite-dimensional latent state exist, with what is necessarily preserved and necessarily lost? Builds on DM-2/DM-3/DM-8 (Phase 4); new results CR-5, CR-9 refereed with corrections incorporated (category labels, full-class hypotheses, completeness discipline on characteristic-class obstructions, and the $(m,L)$-joint phrasing of entropy floors).

---

## 1. Existence: necessary and sufficient conditions — with the honest gap

**Necessary.**
- Latent dimension $m\ge$ class dimension, in **both** the continuous and definable categories (DM-2(a),(b); the definable category is the right home — the floor survives discontinuity there, and only merely-measurable encoders evade it, at the price of C8-collapse and no stability).
- Topological excess where the class is obstructed: $m_{\min}=4>2=\dim$ for the full $\operatorname{Gr}(1,3)$ class (continuous category; full-class hypothesis — excluding degenerate configurations can void the obstruction).
- Entropy floors, **jointly in $(m,L)$** (CR-9, refereed correction): the DM-8 bound $N(\mathcal C,2h)\lesssim(Lr/h)^m$ constrains the *pair* (latent dimension, decoder Lipschitz constant). "No finite $m$" follows only under a declared Lipschitz/norm budget; a finite-parameter family with unbounded decoder constants evades the floor. All impossibility statements of this phase are therefore $(m,L)$-tradeoffs, not bare dimension bounds.

**Sufficient.**
- Compact smooth classes: Whitney embedding at $m=2\cdot\dim$ (with an extension lemma making the decoder total, consistent-valued, and outer-admissible off the embedded image — the three range constraints of the representation definition).
- Definable classes of finite dimension: definable choice/cell structure supplies representations with definable (possibly discontinuous) encoders at $m=\dim$.

**The gap (stated honestly; refereed).** Necessity is entropy-based; sufficiency is proven only for smooth-compact or definable classes — both strictly stronger than "finite entropy dimension." A class of finite metric-entropy dimension that is neither carries **no sufficiency theorem** in this corpus. The existence question is settled on the two proven regimes and open on the gap class.

---

## 2. Minimal information preserved

Exactly the **closed, size-$(k{+}1)$-truncated, covered-window quotient** of the declared class (DM-3 + CP-2): nothing off coverage, nothing beyond the truncation, nothing below closure at $\varepsilon>0$, and nothing about the indexing. Any latent preserving less is lossy on identifiable content (measured in the operational topology of DM-1); any latent claiming more is representing its prior.

## 3. Unavoidable information loss

- **Forward:** representation error $h$ enters the prediction *inside the modulus*: $\le\tfrac12\omega(2\varepsilon+2h)+h$ — nonlinear, possibly infinite (C3); with the margin lemma CR-4, $h$ additionally erodes selector reliability on the definable collar $\{\mu<2h\}$.
- **Converse:** the $(m,L)$-entropy floor (CR-9) bounds achievable $h$ from below for infinite-dimensional classes: loss is not an artifact of a bad encoder but a counting fact.
- **Certificate direction:** loss must be spent *outward* — $\eta$-inflation of envelopes (CR-6) — or certificates become false; there is no loss-free finite representation of an infinite class, only honestly-priced loss.

---

## 4. Identifiable information versus representation choices (CR-5)

The separation demanded by the mandate is exact, and it has a new topological teeth:

**Identifiable:** the decoder image $\beta(z)$ — the induced (truncated, closed, covered) window system — and *only* that. Comparisons of trained systems are meaningful only in the induced-operator metric on the double gauge quotient (DM-9 / NP-6).

**Non-identifiable representation choices:** the latent coordinates themselves (decoder-gauge), the choice of basis/chart, the anchor of any additive decomposition where prior sections are unbounded (MP-3).

**Theorem CR-5 (no global continuous gauge fixing; refereed with disciplines).** Setting: the exactly-$d$ linear class $\cong\operatorname{Gr}(d,N)$, $1\le d<N$ (all $N\ge2$ — the referee removed my $N\ge3$ restriction).
(a) A *basis-style latent* — a continuous global assignment $V\mapsto(b_1(V),\dots,b_d(V))$ of an ordered basis, i.e. a continuous section of the (ordered-basis or, via Gram–Schmidt, orthonormal) Stiefel bundle — **does not exist**: it would trivialize the tautological bundle $\gamma_d$, forcing all Stiefel–Whitney classes to vanish, but $w_1(\gamma_d)$ generates $H^1(\operatorname{Gr}(d,N);\mathbb Z/2)=\mathbb Z/2\ne0$.
(b) Elementary $d=1$ proof: a continuous unit spanning vector on $\mathbb{RP}^{N-1}$ would be a section of the double cover $S^{N-1}\to\mathbb{RP}^{N-1}$; sections of covering maps have clopen image (closed as an equalizer in a Hausdorff total space, open via evenly-covered sheets), contradicting connectedness of $S^{N-1}$ — valid already at $N=2$.
(c) **Category contrast (the refereed phrasing):** latents that factor through a Stiefel lift cannot be globally continuous; gauge-*invariant* embeddings (e.g. $V\mapsto P_V$) can. This is a contrast of representation styles, not a theorem that only invariants are continuous.
(d) **Disciplines:** continuous-category only (definable discontinuous gauges exist — cut along a curve); full-class hypothesis (contractible charts admit continuous local bases — bundles over contractible paracompact bases are trivial); and $w_1$ is a *complete* obstruction only for $d=1$ — for $d\ge2$, vanishing Stiefel–Whitney classes are necessary, not sufficient, for triviality, so CR-5 is an obstruction theorem, not a characterization.

**Consequence for latent design (constraint, not architecture):** globally continuous canonical latents must be gauge-invariant objects; coordinate-style latents are forced to be local (chart-wise, with discrete chart transitions — which is again the branch-plus-selector structure of CR-1, now at the meta level).

---

## 5. Summary table

| Question | Answer | Status |
|---|---|---|
| Does a finite latent exist? | iff the class is finite-dimensional in the proven regimes (smooth-compact / definable); gap class open | refereed |
| Minimal preserved information | closed, $(k{+}1)$-truncated, covered-window quotient | proven (DM-3, CP-2) |
| Minimal dimension | class dimension (cont & def), with topological excess (full-class hypothesis); $(m,L)$-entropy floors for infinite classes | proven (DM-2, CR-9) |
| Unavoidable loss | inside-the-modulus transfer + collar erosion + entropy converse; priced outward or false | proven (DM-8, CR-4, CR-6) |
| Identifiable vs choice | decoder image vs gauge; no global continuous basis-style gauge (CR-5) | proven with disciplines |
