# XP2-G — Mathematical Interface Audit Of The Candidate Section Statistic

Audited against `theory/FINAL_FROZEN_THEORY/` chapters 01–07 and the deployed
`model/config.py` / `model/meta_operator.py`. **This audit is conditional.** It
establishes what the candidate statistic *would* have to satisfy and where it
does or does not satisfy it. It is not an admission, and XP2's empirical Gate
result (see `XP2_FINAL_REPORT.md`) governs whether the statistic may be used at
all.

## 1. The candidate

```
z_section = ( section_center,        continuous, scalar
              outer_radius,          continuous, scalar
              support_rank,          finite, {0 … d}
              query_coverage,        continuous, [0,1]
              inverse_conditioning,  continuous, [0,1]
              validity_flag,         finite, {0,1}
              declared_context )     finite
```

with

```
section_center(S,Q,gamma) = clip( < uhat(Q), vhat(S) >, -c, +c )
vhat(S)                   = ridge solve of the support residuals on the
                            gamma-declared loading design
uhat(Q)                   = chi_u(chemistry of Q), chi_u frozen in gamma
```

## 2. Requirement-by-requirement

| # | frozen-theory requirement | verdict | note |
|---|---|---|---|
| 1 | observable from `S, Q, gamma` | **PASS** | every coordinate is a function of the support multiset, the query structure and gamma-declared maps |
| 2 | no query label | **PASS** | `uhat(Q)` uses the query's *structure*; the query affinity never enters |
| 3 | support permutation invariant | **PASS** | ridge solution, SVD spectrum, rank and projector are all functions of the multiset `{(u_s, y_s)}` |
| 4 | bounded / compact `Z` | **PASS, conditional on declaration** | requires declared `c` (centre clip), `R` (radius cap) and the natural `[0,1]` ranges; `Z` is then a finite union of compact cubes indexed by the finite coordinates |
| 5 | at most `k` continuous support-identified d.o.f. | **PASS** | the support contributes at most `min(k-1, d)` identified directions (§3), and `z` exposes only the scalar `section_center` plus the scalar radius |
| 6 | gauge-invariant outputs | **PARTIAL — requires a declared gauge** | see §4; `section_center` and `support_rank` are invariant, `query_coverage` and `inverse_conditioning` are **not** |
| 7 | explicit epsilon dependence | **PASS, conditional** | the numerical rank tolerance `eps` must be frozen in `gamma`; `support_rank` is a step function of `eps` |
| 8 | section radius / conservative outer enclosure | **PASS, conditional** | requires the two-term radius of §5; a ridge point estimate alone is **not** an enclosure |
| 9 | partiality and abstention | **PASS, and natively representable** | see §6 |
| 10 | fixed deployment state | **PASS, conditional** | `chi_u`, `chi_alpha`, the gauge, `eps`, `tau`, `B` must all be frozen into `gamma`/`D` **before** meta-training `F`, and estimated on a corpus disjoint from the meta-training task law |
| 11 | no modification of CSMO / Band / `K` / mesh | **PASS** | the operator `A(F,z) = K(B(z)F(z))` is untouched; only `d_z` and the view map are re-declared, which `model/config.py` already documents as an engineering choice carrying no theorem |

## 3. The identifiability ledger (measured, not assumed)

The support intercept is unpenalised, so it absorbs the mean of the support
loadings. Only the **centred** support design carries identifiable interaction
directions:

```
identified section dimension  =  rank( U_S - mean(U_S) )  <=  min(k - 1, d)
```

XP2-C measured exactly this:

| `k` | identified dim (measured) | mean query coverage | consequence |
|---|---|---|---|
| 1 | **0.00 / 3** | 0.000 | ridge returns `v = 0`; the arm is *identically* the additive arm |
| 2 | 1.00 / 3 | 0.343 | one direction |
| 3 | 2.00 / 3 | 0.672 | two directions |
| 4 | 3.00 / 3 | ~1.0 | fully identified at `d = 3` |
| 5 | 3.00 / 3 | ~1.0 | fully identified at `d = 3` |

**This is the mathematically decisive constraint of the whole stage.** A frozen
requirement of `k <= 5` caps the identifiable section at `d <= 4`. At `k = 1` no
interaction direction exists at all, and any nonzero prediction would be the
ridge prior speaking, not the support. The registration's instruction —
*"ridge regularization must not be described as identification of unobserved
directions"* — is therefore not a stylistic preference; it is the difference
between `d = 0` and `d = 3` at `k = 1`.

## 4. Gauge

The factorisation `Gamma = U V^T` is invariant under `U -> U A`,
`V -> V A^{-T}` for any invertible `A`. Consequences:

- `section_center = <uhat, vhat>` is **invariant**. Safe.
- `support_rank` is **invariant** (rank survives invertible maps). Safe.
- `query_coverage` and `inverse_conditioning` are defined through an orthogonal
  projector and a singular-value ratio, which are **not** `GL(d)`-invariant.
  They are meaningful only once a gauge is fixed.

**Required remedy:** freeze the gauge in `gamma` by whitening the loading map on
the declared corpus, `Cov(U_train) = I`, and record the whitening matrix as part
of the deployment state. Under that declaration both coordinates become
well-defined functions of `(S, Q, gamma)`.

**Consequence for interpretation:** individual coordinates of `u` and `v` carry
no meaning even after whitening, because whitening fixes the gauge only up to an
orthogonal rotation. Only (i) fixed named features and (ii) gauge-invariant
objects — the scalar `section_center`, the rank, the subspace `span(U_S)` and its
projector — may be interpreted. **No coordinate of the latent section may be
labelled hydrogen-bonding, hydrophobic, DFG or any other biological name.** XP1
already showed the phenomenon empirically: coordinate-wise loading `R2` in a
fold-local gauge was ≈0 or negative while the gauge-invariant reconstruction `R2`
was clearly positive.

## 5. The radius must have two terms

A ridge point estimate silently sets the unidentified component of `v` to zero.
A conservative outer enclosure requires both an estimation term and an
identification term:

```
outer_radius(Q,S) = z_{1-a} * sigma_hat * sqrt( uhat^T (D^T D + lam I)^{-1} uhat )      (estimation)
                  + || (I - P_S) uhat || * B                                            (unidentified)
```

where `P_S` projects onto the centred support span and `B` is a declared bound on
`||v||` over the deployment population. The second term is what makes the
statement honest at `k <= 3`, where coverage is 0.34–0.67: **a third to two
thirds of the query direction is not identified by the support at all**, and the
enclosure must widen accordingly rather than pretend to a point value.

## 6. Abstention is representable inside the frozen operator

`B(z) = [beta_0(z) | beta_1 | … | beta_m]` with `beta_0(z) = b^pop_{kappa(z)}`.
When `validity_flag = 0`, the coefficient map may place all mass on `p_0`, and
the emitted class is exactly the population band for the declared context. So
abstention needs **no** new output object, no new operator and no change to
`K`. It is the `p = e_0` vertex of the existing simplex.

**Placement rule.** `(S-CONT)` requires `L_0` to be uniformly continuous in the
continuous part of `z`. `validity_flag` and `support_rank` are discontinuous by
construction, so they must enter through the **finite context map `kappa`**, not
through the continuous coordinates consumed by the multilinear sieve. The
deployment already has this machinery: `model/meta_operator.py::context_index`
maps declared coordinates into a finite context. This is the correct home for
them, and it keeps `(S-CONT)` intact.

## 7. Dimensional budget

`Gamma_N = O( sqrt( D_N log(Lambda N) / N ) )` with `D_N = (m+1) nu_N` and
`nu_N ~ r_N^{-dim Z}`. The candidate adds **two** continuous coordinates
(`section_center`, `outer_radius`) plus two more if coverage and conditioning are
carried continuously, and two finite coordinates that belong in `kappa`. With
`view_res = 6`, one 2-D CSMO view costs `(6+1)^2 = 49` nodes and
`49 * 8 = 392` parameters. Adding one view for `(section_center, outer_radius)`
is affordable. Carrying the full seven-tuple continuously is **not** the intended
reading and is not recommended.

## 8. What this audit does not establish

It does not establish that the statistic is worth admitting. Interface legality
is necessary and not sufficient; the empirical Gate in `PREREG_XP2.md` §10
governs. It also does not evaluate the emitted probability law: nothing in XP2
scores `K(B(z)F(z))`, its coverage, its width, or its calibration. Interaction
`R2` is not law calibration, and the two must not be conflated.
