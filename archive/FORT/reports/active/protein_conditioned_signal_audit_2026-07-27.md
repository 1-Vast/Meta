# Obtaining a stable, identifiable protein-conditioned signal: audit, literature review, and stop decision

**Date:** 2026-07-27. **Type:** literature audit + dependency audit + substrate audit + power audit +
re-analysis of already-executed runs. **No new labels were read; no model was trained; no download or API
call was made.** `sealed_test_consumed=false`; `panel_development_labels_read=false`;
`confirmation_labels_read=true` (pre-existing).

**Headline:** the recommended next step is **`MODEL_SIDE_ROUTE_NOT_IDENTIFIABLE_STOP`**. Both stated
reopening requirements fail their audit, and the power audit shows the substrate is short of the required
independent-component count by **8×–18×**. This also forces a **correction to my own PARC M0 reading** and
to HQ-GBMA Stage D: both gates were underpowered, so neither is a refutation.

---

## Part A — Evidence classes used in this report

| tag | meaning |
| --- | --- |
| **[LIT]** | literature-supported fact, cited primary source |
| **[THY]** | theoretical inference from a cited result |
| **[PRJ]** | project-specific measured evidence in `history.md` / `reports/active/` |
| **[NEW]** | measured in this session (audits and re-analysis of stored per-component values) |
| **[HYP]** | untested hypothesis |
| **[ENG]** | engineering feasibility |

Nothing below is predictive evidence. Every project number cited is train-only mechanism evidence.

---

## Part B — Audit of the two stated reopening requirements

`task.md` §8.10 states that reopening the pocket-coordinate question requires **(a)** a verified 3Di or
equivalent structure-token representation **and (b)** a suitable multi-family substrate. Both were audited
rather than assumed.

### B.1 Requirement (a) — structure tokens: **NOT SATISFIED** [NEW][ENG]

| probe | result |
| --- | --- |
| `foldseek` python module / binary on PATH | **absent / NOT FOUND** |
| `mmseqs` binary | NOT FOUND |
| `biotite`, `Bio` (biopython), `fair_esm`, `esm` | **all absent** |
| local `.pdb` files under `dataset/` | **0** |
| local `.cif` / `.mmcif` files under `dataset/` | **0** |
| `dataset/structure/alphafold` | **missing** (removed in the 2026-07-25 cleanup) |
| HuggingFace hub cache | ESM-2 (8M/150M/650M), ChemBERTa, prot_bert, MiniLM, assorted GGUF — **no SaProt, no ProSST, no ProstT5** |

Producing a 3Di or quantized-structure coordinate would require **two separate downloads**: (i) ~10³
AlphaFold monomers, and (ii) a Foldseek binary or a structure-token PLM. Neither exists locally, and both
are user-authorised actions. Requirement (a) is **open**, not met.

Note also **[LIT]**: both leading structure-token PLMs require 3D structure *at inference*, so neither
removes the dependency. SaProt combines Foldseek 3Di tokens with residues into a structure-aware
vocabulary over ~40M AlphaFold structures
([bioRxiv 2023.10.01.560349](https://www.biorxiv.org/content/10.1101/2023.10.01.560349v5.full)); ProSST
quantizes residue-level local micro-environments into discrete tokens with disentangled attention over
18.8M structures ([NeurIPS 2024](https://proceedings.neurips.cc/paper_files/paper/2024/file/3ed57b293db0aab7cc30c44f45262348-Paper-Conference.pdf)).
Both are benchmarked on **mutation-effect and generic protein tasks**, not on target-conditioned ligand
reordering, and the independent
[VenusMutHub](https://www.sciencedirect.com/science/article/pii/S2211383525001650) evaluation finds the
ranking between them is dataset-dependent and that ProSST is preprocessing-sensitive (one assay
0.54 → 0.08). Under the project's own rule — a generic protein-benchmark improvement is not DTA evidence —
neither model currently carries admissible evidence for `g(t,d)`.

### B.2 Requirement (b) — multi-family substrate: **NOT SATISFIED** [NEW][PRJ]

Local substrates with a built registry:

| substrate | families | targets | components | median ligands/target | interaction identified? |
| --- | --- | ---: | ---: | ---: | --- |
| ChEMBL-37 dual-cold registry | **multi-family** | 1,045 | 965 | 86 (723 targets ≥40 lig) | **NO** — `P0_CYCLE_A` projected label SD **0.356 < 0.5** floor, 51.45% of residual energy on 1% of ligands |
| Metz `panel_metz` | kinase only | 112 | 101 | 43 | **YES** — `PA2` 0.663, `PA5` p=0.000488 |
| Davis `panel_davis` | kinase only | 116 | 102 | 12 | sealed, `MDE80 0.1596` |
| Reinecke `panel_reinecke` | kinase only | 109 | 104 | 5 | `MDE80 +0.0668` |
| KirHub | kinase only | 377 | 324 | 44 | single-source, one-concentration |
| SPD (multi-family) | GPCR/kinase/channel/enzyme/NR/transporter | 101 genes | ~80–100 | **14** | `STANDARDIZED_PANEL_NOT_DUAL_COLD_CAPABLE` |
| ToxCast, Papyrus | multi-family | — | — | — | `TCOPA_G0_...STOP`, `PAPYRUS_F0_...STOP` |

**The binding fact is a clean dichotomy: every substrate on which the interaction is identified is
kinase-only, and the only multi-family substrate with per-target depth fails the identifiability floor.**
Requirement (b) is not met, and this is the same measurement-design limit the program has already reduced
to — not a new one.

---

## Part C — Failure diagnosis (§6.1), including a correction to my own result

### C.1 PARC M0 was **underpowered**. Its stop is not a refutation. [NEW]

The power audit uses the per-component values stored by the executed run (`reports/active/parc_m0.json`),
with the component as the inference unit.

```
contrast                  n    mean      raw SD   MDE80(raw)  |  MAD-sigma  MDE80(robust)
parc - shared_global     77   -0.0956    0.8436     0.2693    |   0.2663      0.0850
esm  - shared_global     77   +0.5652    6.7723     2.1622    |   0.2621      0.0837
parc - random_positions  77   +0.2712    1.6841     0.5377    |   0.2349      0.0750
```

Components required to resolve the **frozen 0.02** advantage:

```
effect 0.02  ->  n = 13,966 (raw SD)      1,391 (robust scale)
effect 0.03  ->  n =  6,207               618
effect 0.05  ->  n =  2,234               223
available on the only identified substrate (Metz):  77
```

**PARC M0's gate was short by a factor of 18× (robust) to 181× (raw) in component count.** A genuine
+0.02 containment effect could not have been detected. The verdict code
`PARC_M0_COORDINATE_NOT_LOAD_BEARING_STOP` stands as a *frozen procedural outcome* — the gates were
preregistered and they failed — but **its scientific reading must be downgraded to
`UNDERPOWERED`**. I recorded it as a refutation of the coordinate hypothesis; that reading was wrong and
is corrected here.

The same arithmetic applies to **HQ-GBMA Stage D**, which used the same estimand, the same substrate and
~78 components. Its `FAIL_STOP` is likewise an underpowered gate, not a refutation.

**Limitation of my own G0 control, disclosed:** the synthetic estimator-sensitivity control used
`V_t = sigma^2 I` with a single matched `sigma`, so its denominator was well-conditioned. It proved the
estimator is sensitive *under well-behaved noise*; it did **not** prove sensitivity under the real
heavy-tailed signal distribution. G0 is therefore weaker evidence than I claimed.

### C.2 What *is* resolvable: sign and signed-rank statistics [NEW]

Mean-based contrasts are hopeless here, but ordinal statistics on the same 77 components are far better
powered. Recomputed from the stored per-component values (no new labels):

| contrast | median | fraction positive | sign test p | Wilcoxon p |
| --- | ---: | ---: | ---: | ---: |
| **`esm − shared_global`** | **−0.1165** | 0.31 | **0.0013** | **0.0003** |
| `parc − shared_global` | −0.0243 | 0.39 | 0.068 | 0.079 |
| `parc − esm` | +0.0253 | 0.55 | 0.494 | 0.418 |
| `parc − random_positions` | +0.0317 | 0.57 | 0.254 | **0.0171** |
| `parc − wrong_target` | +0.0269 | 0.58 | 0.171 | 0.059 |
| `parc − random_features` | +0.0246 | 0.58 | 0.171 | 0.232 |
| `parc − pocket_composition` | +0.0079 | 0.53 | 0.649 | 0.857 |
| `esm − esm_wrong_target` | +0.0129 | 0.52 | 0.820 | 0.694 |

Three findings survive at this resolution, and they are the substantive scientific content of the round:

1. **Pooled ESM-2 is significantly *worse* than the protein-free shared-global basis** (median −0.117,
   sign p = 0.0013, Wilcoxon p = 0.0003). Stage D's *direction* is confirmed with a valid test even though
   its magnitude is not.
2. **Pooled ESM-2 is statistically indistinguishable from its own exposure-matched derangement**
   (p = 0.82). It carries no recoverable target-specific subspace information on this substrate.
3. **Aligned position adds nothing over pocket amino-acid composition** (p = 0.86). Whatever weak pocket
   signal exists is *compositional*, not *positional* — which specifically undercuts the residue-locality
   hypothesis that motivated PARC.

The one positive whiff — `parc − random_positions`, Wilcoxon p = 0.0171 — does **not** survive Holm
correction across the eight contrasts (0.0171 × 8 = 0.137) and disagrees with its own sign test
(p = 0.254). It is not evidence.

### C.3 Two estimator defects, one of which is new [NEW][LIT]

**Defect 1 — mean of unbounded ratios.** `containment = inside/signal`, `signal = beta_hat^T beta_hat −
tr(V_t)`, min positive signal 0.00309 vs median 0.14730 (**47.7×**). One component reached 56.09 and moved
the ESM arm's mean from 0.151 to 0.878 alone. **[LIT]** The comparative study of ratio uncertainty
concludes in favour of the **ratio of means (ratio of sums)** over the mean of ratios on lower statistical
uncertainty ([arXiv:1409.4896](https://arxiv.org/pdf/1409.4896)); classical survey sampling agrees, with
an `O(1/n)` bias that is diagnosable through the denominator's coefficient of variation. Fieller's method
is the exact interval and, correctly, returns an **unbounded** interval when the denominator is not
significantly non-zero ([Franz, arXiv:0710.2024](https://arxiv.org/pdf/0710.2024)) — an honest signal
rather than a numerical failure. Guidance is to prefer Fieller once the denominator CV exceeds ~15%.

**Defect 2 (new, and worse) — selection on a noisy denominator.** `held_containment` retains a target only
if `signal[t] > 0`. **22.5% of targets have a non-positive noise-corrected signal** (positive-signal
fraction 0.7748; the signal 5th percentile is −0.640 and its minimum is −8.669). Conditioning on
`signal > 0` is selection on a *noisy* variable that is mechanically correlated with the numerator, so the
retained containment fractions are **biased upward by construction**, and the bias differs by arm because
each arm's numerator differs. This is a bias, not just variance, and it affects every containment number
in `hqgbma_stage_d.json` and `parc_m0.json`.

**Consequence:** no containment mean anywhere in this program is a usable effect size, and the fix is not
only aggregation but also removing the denominator-based selection (e.g. by aggregating
`sum_t inside_t / sum_t signal_t` over *all* targets, which is well-defined because negative
noise-corrected signals cancel in expectation rather than being discarded).

### C.4 The mechanical cause of the protein-map failure: degrees of freedom [NEW][THY]

| object | free parameters | training targets | params per target |
| --- | ---: | ---: | ---: |
| `ProteinGrassmann(32 → 64, r=6, hidden=64)` (Stage D & M0) | **27,072** | ~89 | **≈ 304** |
| shared-global subspace `Gr(64, 6)` | 348 (manifold dim) | ~89 | ≈ 3.9 |
| MLP correction `Delta = W2 σ(W1 u)`, `32→16→6` (my Part 8 §8.2) | 630 | ~89 | ≈ 7.1 |
| **rank-1 gated correction** `w̄ + γ (a^T u_t) v` | **39** | ~89 | **≈ 0.44** |

The protein-conditioned arm was fitting **~304 free parameters per training target**. The shared-global
arm is not a "weaker" model — it is the only one in the comparison whose estimator is regularized (a
spectral eigendecomposition of a pooled scatter). **[THY]** This alone predicts the observed ordering
without any biology: an over-parameterized amortized map trained on ~89 units will lose to a pooled
spectral estimate. **My own Part 8 §8.2 design (630 params) is also over-parameterized for this substrate**
and would have inherited the same problem. That is a design error in the blueprint I proposed, corrected
in Part E.

### C.5 Full decomposition of the PARC failure (§6.1)

| cause | present? | evidence |
| --- | --- | --- |
| **Estimator instability** | **yes, dominant** | 47.7× denominator range; one component = 56.09; mean sign flips under one-target removal |
| **Denominator selection bias** | **yes, new** | 22.5% of targets discarded on a noisy `signal > 0` filter |
| **Small effective component count** | **yes, dominant** | 77 components vs 618–1,391 required |
| **Excessive output degrees of freedom** | **yes** | 27,072 params / 89 targets ≈ 304 per target |
| **Static target-only coordinate** | untested | ligand-conditioned coordinate never evaluated |
| **Representation misspecification** | **partly refuted** | positional pocket adds nothing over composition (p = 0.86) |
| **Noisy side information** | **yes** [THY] | dirty-IMC: inexact row features restore linear-in-`n` d.o.f. |
| **Kinome-only generalization limit** | **yes, unresolved** | every identified substrate is kinase-only |
| **Measurement-design limitation** | **yes, binding** | no multi-family identified substrate exists |

---

## Part D — Literature review and candidate comparison (§6.2)

### D.1 Noisy / misspecified side information [LIT]

* Exact side information gives `O(log n)` sample complexity for two-sided cold start (Xu et al. 2013;
  Jain & Dhillon 2013). With **noisy or partial** side information the degrees of freedom remain **linear
  in n** and the advantage is lost — "dirty IMC"
  ([Chiang, Hsieh & Dhillon, NeurIPS 2015](https://www.cs.utexas.edu/~inderjit/public_papers/dirtyIMC_nips15.pdf)),
  which models the target as *side-information interaction plus a low-rank residual* — i.e. exactly the
  `G_shared + G_bio + G_latent` decomposition proposed in the brief.
* [Yang & Ma, arXiv:2605.17189](https://arxiv.org/abs/2605.17189) close the noisy-IMC gap and show the
  error degrades **with the level of subspace misspecification**, proposing a **penalized interpolation
  between IMC and ordinary MC** — the mathematical justification for a shrinkable `gamma`.
* [Fine-grained generalisation analysis of IMC (NeurIPS 2021)](https://proceedings.neurips.cc/paper/2021/file/d6428eecbe0f7dff83fc607c5044b2b9-Paper.pdf)
  gives bounds depending only on the side-information dimension `d`, and notes they become **vacuous** as
  the features lose information — the formal statement of "biological features must be shrinkable to zero".
* Jalali et al. (2010) "dirty model" (shared + task-specific components) is the multi-task ancestor;
  ["Multi-Task Learning Via Sharing Inexact Low-Rank Subspace"](https://ieeexplore.ieee.org/iel7/9413349/9413350/09414782.pdf)
  is the direct analogue of our situation — an inexact shared subspace with per-task deviation.

**[THY] Conclusion:** the `G_shared + γ·G_bio + G_latent` decomposition is the *theoretically correct*
form. It does not, however, reduce the number of independent units required to estimate `γ`.

### D.2 Partial pooling and shared/task-specific subspaces [LIT]

Reduced-rank regression is the statistical core (Anderson 1951; Izenman 1975), with nuclear-norm
adaptivity (Yuan et al. 2007; Chen et al. 2013) and joint rank/variable selection (Bunea et al. 2011/2012).
`SMART`/common-mechanism regression assumes a **shared left low-rank component with per-task right
components**, "dramatically reducing the number of samples needed", solved by a non-iterative spectral
algorithm ([arXiv:2604.20161](https://arxiv.org/pdf/2604.20161)) — a direct, better-conditioned alternative
to an amortized neural map.

**The decisive sample-complexity fact [LIT]:**
[meta-learning of shared linear representations](https://arxiv.org/html/2501.18975) shows that in the
low-data-per-task regime, recovering a shared subspace requires **the number of tasks to scale
exponentially in the subspace dimension**. With `r = 6` and 77 components, we are far outside the regime
where any of these estimators is consistent.

### D.3 Manifold-valued regression [LIT]

Grassmannian geodesic regression (Hong et al., ECCV 2014) fits a geodesic parameterized by a base point
and a velocity — formally the `R_t = Exp_{R_0}(γ Δ_t)` form in the brief. Uncertainty-aware variants exist
([Gaussian Process Subspace Regression, arXiv:2107.04668](https://arxiv.org/pdf/2107.04668)), which is
explicitly data-efficient. But the field itself records that regression with manifold-valued *responses*
is far less developed than with manifold-valued predictors, precisely because there is no linear structure
for averaging responses.

**[THY] Recommendation: drop the Grassmann formulation.** The manifold-valued response is the *source* of
the ratio pathology — containment is a normalized energy fraction, i.e. a ratio with an unbounded
denominator, and it exists only because the estimand is a subspace rather than a vector. A Euclidean
residual on the coefficient itself (`w_t` vs `w̄`) admits ordinary paired inference with no ratio, no
denominator selection and no manifold. Retaining Grassmann because the previous experiment used it would
be exactly the error the brief warns against.

### D.4 Structure-informed protein representation [LIT]

Covered in B.1. Additionally, the closest sister field has converged on the *opposite* of what the pocket
hypothesis assumes about pooling: enzyme-specificity work argues substrate scope is set by "a select few
residues within the active site" so whole-enzyme encoders "stifle their signal"
([arXiv:2607.05084](https://arxiv.org/html/2607.05084v1)) — yet our own measurement finds **position adds
nothing over composition** (p = 0.86) on the identified panel. **[HYP]** These are reconcilable only if the
relevant residues are ligand-dependent (D.5) or if the kinase panel's pocket variation is too small
(mean pairwise pocket identity 0.373 across 111 distinct pockets — i.e. there *is* variation to exploit).

### D.5 Ligand-conditioned protein coordinates [LIT]

The literature does **not** support treating unsupervised cross-attention as a way to recover the relevant
residues. ArkDTA
([Bioinformatics 39:i448](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i448/7210465))
shows earlier attention-based DTA models could not separate active from inactive residues, and adds
explicit **non-covalent-interaction supervision** to fix it — i.e. it works *because* of pose-derived
labels, which this program does not have and whose route is closed. Cold-start DTA papers reporting large
gains (CS-DTA; XAttn-DTA, "MSE reductions up to 79.0%") evaluate on Davis/KIBA/BindingDB splits that
[target mirroring](https://www.biorxiv.org/content/10.64898/2026.06.29.735309v1) and
[Graber et al.](https://www.nature.com/articles/s42256-025-01124-5) show are leaky; their authors also
disclaim that attention hotspots are "prioritized candidates … rather than definitive assertions of a
precise physical binding pose". **[THY]** Not admissible evidence, and it adds degrees of freedom exactly
where we are starved.

### D.6 Robust statistics and small-component inference [LIT]

Prefer **ratio of sums**; diagnose the denominator CV; use **Fieller** when the CV is large (accepting
unbounded intervals as honest); and note that **cluster/multilevel bootstraps for ratios have documented
coverage failures** — a systematic evaluation found all four bootstrap procedures gave inaccurate standard
errors and below-nominal coverage for cluster-specific ratios
([PMC9564085](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9564085/)). This last point matters: the
program's component bootstrap has been applied to ratio estimands throughout, and its coverage has never
been validated for that use.

### D.7 Ranked candidate comparison (§6.2)

Effective d.o.f. are per §C.4; "units needed" uses the robust MDE scale from §C.1 at a 0.02 effect.

| rank | candidate | formulation | eff. d.o.f. | units needed | behaviour under wrong biology | cold-target inference | nested null? | verdict |
| ---: | --- | --- | ---: | ---: | --- | --- | --- | --- |
| **1** | **Stop (no model-side route)** | — | 0 | 0 | — | — | — | **ADOPT** |
| 2 | Rank-1 gated Euclidean correction | `w_t = w̄ + γ(a^T u_t)v` | 39 | ~1,391 | `γ→0` exactly | fully available | **yes, exact at γ=0** | best *if* a substrate existed |
| 3 | Dirty partially-pooled residual | `w̄ + γΔ(u_t) + λ_t` (latent train-only) | 630 + n·r | ≥1,391 | `γ→0`; latent absorbs misfit | latent **unavailable** → reduces to 2 | yes | theoretically correct, **not estimable here** |
| 4 | Simple regularized bilinear residual | `g = φ(d)^T B u_t` | 32×6 = 192 | ~1,391 | no shrink-to-zero unless gated | available | only if gated | acceptable fallback |
| 5 | Tangent-space / Grassmann correction | `R_t = Exp_{R_0}(γΔ_t)` | ≥348 | ≫1,391 | `γ→0` | available | yes | **REJECT** — manifold response *causes* the ratio pathology |
| 6 | Ligand-conditioned coordinate | `z_{t|d} = Σ_j α_j(H_t,d)H_{t,j}` | 10³–10⁴ | ≫1,391 | no guarantee | available | no | **REJECT** — needs NCI supervision (ArkDTA); closed route |
| 7 | Structure-distilled coordinate (SaProt/ProSST) | replace `u_t` | same as 2–4 | ~1,391 | — | needs 3D at inference | — | **BLOCKED** (B.1) + generic-benchmark evidence only |
| 8 | Full amortized subspace map (Stage D / M0 as run) | `R_t = f(u_t)` on Stiefel | 27,072 | ≫≫ | overfits | available | no | **CLOSED** — 304 params/target |

---

## Part E — Final decision (§6.3) and the design that *would* be correct

### E.1 Decision: **`MODEL_SIDE_ROUTE_NOT_IDENTIFIABLE_STOP`**

The §6.4 stop rule fires: *"If the current substrate cannot distinguish the proposed model from
shared-global at the frozen effect floor, stop before implementation."* It cannot, by 8×–18× in units.
Requirements (a) and (b) both fail. **No further model-side implementation is authorised.** The
scientifically justified next step remains the data track — LEXOR, or a prospective multi-family factorial
panel.

### E.2 The design to run **if and only if** a qualifying substrate appears [HYP]

Recorded so it is preregisterable later, not to be built now. Two innovations, both acting during
training; shared-global is an exact nested null; the biological term is exactly shrinkable to zero.

* **Protein representation:** `u_t ∈ R^{d_u}`, `d_u ≤ 32`, from a structure-token encoder *if* B.1 is
  satisfied; otherwise pocket composition (position is not load-bearing, §C.2).
* **Ligand representation / interaction operator:** unchanged — `φ(d)` from the exact ORRC observed-edge
  projection `M_X^W`, `b(d)` untouched.
* **Shared component:** `w̄`, spectral estimate on pooled cross-fitted coefficients.
* **Biological correction (innovation 1):** **rank-1 and scalar-gated**,
  `w_t = w̄ + γ·(a^T u_t)·v`, `a ∈ R^{d_u}`, `v ∈ R^r`, `‖v‖ = 1`. **39 free parameters**, not 27,072.
  `γ = 0` recovers shared-global *exactly*.
* **Train-only latent nuisance:** `λ_t` with a strong prior, present in the training likelihood only and
  **structurally absent at inference**, so it cannot create a cold-start shortcut (dirty-IMC form).
* **Regularization/prior:** `γ` under a horseshoe-type global–local prior; `λ_t` under an
  error-corrected empirical-Bayes prior with the sandwich covariance `V_t` subtracted.
* **Estimator (innovation 2):** all component-level ratio estimands aggregated as **ratio of sums over all
  targets** (no `signal > 0` selection), reported with denominator CV, Fieller interval, median, trimmed
  mean, positive-component fraction, and leave-one-component-out sensitivity. Primary inference by
  **signed-rank / sign test** on components, not by a mean.
* **Objective:** weighted least squares on the projected residual + within-target pairwise ranking
  auxiliary; cross-fitted on held homology components.
* **Inference / `k=0`:** `y = b(d) + φ(d)^T(w̄ + γ̂(a^T u_t)v)`; no support labels; the Part 1.3 contract
  holds vacuously. **Few-shot:** deferred; `Δ_info = +0.0154 [−0.0155,+0.0464]` says the channel is below
  the usable floor.

### E.3 Stability analysis for that design (§6.4) [NEW]

Even at 39 parameters, the binding constraint is unchanged: `γ` is a **scalar contrast between two nested
models**, and its resolution is governed by the component count and the per-component dispersion, not by
the parameter count. At the measured robust dispersion, **~1,391 independent components** are needed for a
0.02 effect and **~618** for 0.03. **Minimum admissible substrate: ≥ 600 independent homology components,
multi-family, with the interaction identified (`PA2 ≥ 0.5 pK`) and ≥ 40 scaffold-diverse query ligands per
target after all firewalls, including the binding-profile-correlation axis.** No such substrate exists
publicly; this is a stricter and better-quantified version of the standing `NO_OPEN_POWERED_INDEPENDENT_PANEL`
requirement, and it is the concrete specification the prospective panel must satisfy.

---

## Part F — Destruction controls to carry forward (§6.5)

Retained for any future gate; `A2`, `A8`, and the two new ones marked ★ are the ones that were missing:

```
true protein vs target shuffle                    | true protein vs matched wrong target
structure tokens vs sequence-only matched capacity| native pocket vs same-protein random positions  (A2)
ligand-conditioned vs target-only pooling         | side-information arm vs latent-only dirty control
native coordinate vs random features              | binding-profile-correlation firewall
provenance-family-disjoint sensitivity            | synthetic gamma-recovery under REAL noise  ★ (A8+)
                                                  | denominator-selection sensitivity: with vs without
                                                  |   the signal>0 filter                       ★
```

★ `A8+` corrects my G0: the synthetic control must reuse the **empirical `V_t`**, not `sigma^2 I`.

---

## Part G — Routes that remain closed (§9.7)

Unchanged and not reopened: physical pose / pocket / docking rescues (4 framings); atom–residue and
six-channel interaction fields; posterior and support-kernel rescues; free expert routing; unrestricted
meta-learning; unregistered Transformer growth; family-selectivity anchors (SAFSA); MMP transformation
anchors; ToxCast target-contrast pretext; Papyrus aggregated anchors; within-kinome group/family
resolution (TR); pan-family single-cold under document isolation (PFSC); aggregated-database
cross-document replication (RECRO L0); protein-conditioned precision and covariance priors; the amortized
Stiefel subspace map. **Newly closed here:** reparameterizing `B` on pooled ESM-2 (already §2.5); the
Grassmann/containment estimand as a *gating* statistic; and any model-side gate on a substrate with
< 600 independent components.

---

## Part H — What would change this decision

1. A substrate with ≥ 600 multi-family independent components at `PA2 ≥ 0.5` — via LEXOR L4/L5 or a
   prospective panel. **This is the only route that changes the answer.**
2. Authorised acquisition of (i) AlphaFold monomers for a multi-family target set and (ii) a Foldseek
   binary or ProstT5/SaProt/ProSST weights — which would satisfy B.1 but **not** B.2, and therefore still
   would not authorise a gate.
3. A different estimand with materially lower per-component dispersion. The ratio-of-sums fix improves
   variance but does not create 600 components.
