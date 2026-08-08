# XP1 — Crossed-Panel Identification Of The Protein-By-Ligand Affinity Interaction

Date: 2026-08-08
Preregistration: `research/crossed_panel_identification/PREREG_XP1.md`
Environment: `drug` conda env, Python 3.11.15, numpy 1.26.4, torch 2.6.0+cu124,
CUDA on RTX 4060 Laptop (used for the ESM-2 encoder only; every inference in this
report is a linear-algebra computation on CPU).

---

## 0. One-paragraph answer

The bottleneck is **not** missing data and **not** the frozen mathematics. On a
properly crossed panel the protein-by-ligand affinity interaction `gamma_ij` is
large (≈60% of the affinity variance), low-rank (rank 1–3 carries most of it),
reproducible across disjoint compound halves (`r = 0.885`) and across two
completely independent measurement platforms (`r = 0.565`), and it **transfers to
proteins from entirely unseen kinase groups** — but only when the protein's
interaction coordinate is identified from a handful of labelled support
observations. Every protein *representation* tested — including the production
ESM-2 t30 encoder, the structure-aligned 85-residue KLIFS binding pocket, pocket
physicochemistry, KLIFS conformational-state availability, family and group
labels, and homolog-kernel averaging — recovers **essentially none** of that
coordinate once near-homologs are removed from training, while the same features
recover **more than half** of it when near-homologs are present. The failure is
therefore located precisely: MetaSieve's zero-shot protein pathway carries
homolog-interpolation information (partner compatibility), not target-specific
affinity direction; the affinity direction is a **support-identifiable
low-dimensional section**, which is exactly the object the frozen theory's
`z(S,Q,gamma)` is designed to hold.

**Final decision: `OBJECTIVE_OR_PARAMETERIZATION_FAILURE`**, with a secondary
`REPRESENTATION_FAILURE` verdict confined to the zero-shot protein-feature
pathway. `H2` (data identifiability failure) is rejected.

---

## 1. Current diagnosis — what is and is not identified

Retained unchanged from the repository ledger, now with XP1 evidence attached:

```text
GEOMETRY_IDENTIFIED                                  (P1B, unchanged)
PARTNER_COMPATIBILITY_PARTIALLY_IDENTIFIED           (unchanged; XP1 localises it)
COMPONENT_ALGEBRA_IMPLEMENTED                        (unchanged)
```

New, established by XP1:

```text
CROSSED_INTERACTION_EXISTENCE_IDENTIFIED             (XP1-A, two panels, two platforms)
INTERACTION_IS_LOW_RANK_r1_TO_r3                     (XP1-A2, XP1-B S2)
SUPPORT_IDENTIFIED_INTERACTION_SECTION_TRANSFERS     (XP1-B A4, group closure)
ZERO_SHOT_PROTEIN_FEATURE_INTERACTION_MAP_NOT_IDENTIFIED
PROTEIN_SEQUENCE_CARRIES_HOMOLOG_INTERPOLATION_ONLY  (XP1-D1)
TRUNCATION_ARTEFACT_FALSIFIED                        (XP1-E)
```

Still **not** established, and not claimed here:

```text
AFFINITY_ENERGETICS_NOT_IDENTIFIED
UNSEEN_LIGAND_LOADING_NOT_TESTED
BIOLOGICAL_Z_NOT_ADMITTED_TO_PRODUCTION
END_TO_END_DTA_CLAIM_NOT_SUPPORTED
```

The critical distinction the project asked for is now measured rather than
argued. Writing `V` for the protein-side interaction coordinate:

| question | answer | evidence |
|---|---|---|
| does `gamma` exist? | yes, ≈60% of variance | XP1-A1 |
| is `gamma` reproducible? | yes | XP1-A3 (`r=0.885`), XP1-A4 (`r=0.565` cross-platform) |
| is `gamma` low-dimensional? | yes, rank 1–3 | XP1-A2, XP1-B rank sweep |
| can `V` be read off protein features? | **no** at group closure | XP1-B A3, XP1-D1 |
| can `V` be read off protein features when a homolog is in training? | **yes** | XP1-D1 (`R2 = 0.52` leave-one-protein-out) |
| can `V` be identified from `k=16` support labels? | **yes**, and it is target-specific | XP1-B A4 vs A6 |

---

## 2. Theory constraints — what the frozen theory permits and requires of `z`

Read from `theory/FINAL_FROZEN_THEORY/` (chapters 01–07) and
`model/config.py` / `model/meta_operator.py`.

**2.1 What `z` must be.** `z = z(S,Q,gamma) in Z` with `Z` a compact metric
domain represented by a finite union of compact cubes. `z` is computed from the
observable support set `S`, the query `Q` and the declared specification
`gamma` — at meta-training *and* at inference. The query label never enters.
The statistic must be deterministic, bounded, measurable, and permutation
invariant in `S`.

**2.2 Sufficiency is not required; it sets the risk floor.** The theory defines
`L_0(z, beta) = E[L(beta, Y) | zeta = z]` and targets
`g_mu^star(z) = argmin_p J_mu(z,p)`. A `z` that discards affinity-relevant
information is still *legal*: it simply moves the achievable risk floor. Only
`(S-CONT)` — uniform continuity of `L_0` in `z` — is assumed. **Hence
"informative feature" and "admissible statistic" are different predicates, and
neither implies the other.** Admission is an empirical question, which is why
the repository's Gate discipline is correct.

**2.3 Dimension is the binding resource.** The hypothesis class is a multilinear
sieve on a mesh of `Z`: `nu_N` nodes, `D_N = (m+1) nu_N` parameters, and
`Gamma_N = C_0 (Lbar + mu/2) sqrt((D_N log(Lambda N) + log(1/delta_N))/N)`.
Since `nu_N ~ r_N^{-dim Z}`, **generalization cost is exponential in `dim(Z)`**.
The deployed configuration (`model/config.py`: `d_z=28`, `m=7`, `M=32`,
`view_res=6`) only survives this because CSMO never builds the literal 28-D mesh;
it mixes fixed low-dimensional views (`DEFAULT_VIEWS` = six 2-D projections).
So an admissible biological coordinate must be **one or two bounded scalars that
can be given their own CSMO view**, not a 640-D or 1280-D embedding. XP1's
finding that the interaction section is rank 1–3 is therefore not a convenience:
it is the difference between an admissible and an inadmissible statistic.

**2.4 What `B(z)` can and cannot do.** `B(z) = [beta_0(z) | beta_1 | ... |
beta_m]` with `beta_0(z) = b^pop_{kappa(z)}` and `kappa` finite-valued; all
continuous `z`-dependence flows through `F(z) in Delta_m`. The emitted law class
is always inside the convex hull of `m+1 = 8` fixed bands. A biological
coordinate can only **move mixture weights**. It cannot introduce a new shape.
This is compatible with an affinity-*location* statistic (which is what
`gamma_ij` is, in log units) and incompatible with using `z` to smuggle in an
unrestricted regression.

**2.5 A protocol constraint that follows, and that XP1 makes concrete.** The
statistic recommended in §9 uses a fixed ligand-loading matrix `U`. `U` is
estimated from data, so it must be **frozen into the declared specification
`gamma` / the deployment `D` before meta-training `F`**, and estimated on a
corpus disjoint from the meta-training task law. If `U` were re-estimated on the
same tasks used to fit `F`, the fixed-deployment premise of chapters 01 and 04
would no longer hold and none of the calibration chain would apply. This is a
hard requirement, not a hygiene preference.

---

## 3. Dataset audit

### 3.1 Why the existing corpora could not answer the question

`E-AFF-X0-FEAS` already proved the point structurally: an X0 rectangle needs two
proteins inside one document-keyed panel, and D1 homology-document closure unions
every pair of targets sharing a document, so **crossing and document-disjointness
are produced by opposite kinds of study**. ChEMBL-style corpora are unions of
single-target SAR papers; they are wide in ligands and shallow in crossing. No
amount of label access repairs that.

The construction that *does* produce crossing is a **panel**: one laboratory, one
assay, one complete ligand × protein rectangle. XP1 therefore acquired panels.

### 3.2 Acquired releases (all checksums frozen in the preregistration)

| ID | content | shape | crossing | verdict |
|---|---|---|---|---|
| `P-METZ` | Metz 2011 kinase `pKi` | 704 compounds × 172 kinases, **complete** | 100% nominal | **primary**; single lab, single assay family ⇒ assay identity is constant by construction |
| `P-KLAEGER` | Klaeger 2017 kinobeads `pKd_app` | 222 drugs × 343 kinases, complete | 100% nominal but **93.3% at the floor** | replication only (binarised) |
| `P-PDSP` | NIMH PDSP `Ki` | 98,678 rows, 274 targets, 16,698 ligands | 32,673 human cells, 3,008 replicated | **secondary**, independent protein class + assay technology |
| `A-KLIFS` | 85-residue aligned kinase pockets, family, group, UniProt | 521 human kinases | — | protein annotation |

Excluded by governance and never downloaded: DAVIS (prohibited), Anastassiadis
and PKIS2 (consumed development panels). No ChEMBL37 affinity value was read;
the `X1`/`X2` prohibition is untouched.

### 3.3 Censoring — the one serious threat, and how it was handled

- `P-METZ` is left-censored at `pKi = 4.0`: **59.16%** of all cells sit exactly at
  the floor. Analysis block `BLK-METZ-60` (greedy peel to ≥60% uncensored) =
  **704 compounds × 82 kinases, 34,764 uncensored cells**.
- `P-KLAEGER` is 93.27% at its floor — unusable as a continuous panel.
- `P-PDSP` carries an explicit flag: 34,300 rows `>` and 55 `<`, all dropped.
  Human uncensored crossed core: **2,344 ligands × 45 targets, 10,701 cells**,
  16 MMseqs2 40%-identity clusters.

Conditioning on `y > floor` is conditioning on the outcome and **can by itself
manufacture apparent non-additivity**. XP1-E (§6.5) is the destructive control
for exactly this, and it clears the result.

### 3.4 Measurement-noise ceiling

PDSP contains 2,490 cells with ≥2 independent literature reports. Splitting the
reports in half: `sd(h1-h2) = 1.0103` log units ⇒ **per-report `sigma ≈ 0.714`
log units**, `r(h1,h2) = 0.727`. Literature `Ki` is far noisier than a single
panel; this is why PDSP effect sizes in §6.4 are small and why Metz is primary.

---

## 4. Estimand and arm structure

```
y_ij = mu + alpha_i (ligand) + beta_j (protein) + gamma_ij (interaction) + eps_ij
```

The target is `gamma_ij`. Every arm predicts the **same** held-out cells of the
**same** held-out proteins, shares the **same** `mu` and `alpha` (fitted on
training proteins only, stage 1 of a two-stage fit so that no arm gets a
rank-dependent main effect), and shares the **same** ligand-loading basis `U`
(fitted on training proteins only, stage 2). Arms differ **only** in how the
held-out protein's interaction coordinate `v_j` is obtained.

| arm | `v_j` from | control role |
|---|---|---|
| `A0` | — (`mu` only) | population |
| `A1` | — (`mu + alpha_i`) | **ligand-only `f(l)`** (control A) |
| `A2` | — (`+ support intercept`) | additive + support location — the baseline for interaction; also the **protein-ID shortcut control** (control E) |
| `A3::<rep>` | ridge map from protein representation `<rep>` | **zero-shot correct protein `f(l,p)`** (control B) |
| `A3B::knn` | homolog-kernel average of training interaction columns | strongest natural zero-shot baseline |
| `A4` | ridge solve on `k` support residuals | **few-shot interaction section** |
| `A5::<rep>` | features of a **deranged** held-out protein | **wrong-protein control** (control C) |
| `A6` | support residuals of a **deranged** held-out protein | **permuted-support control** (control D) |
| `A7` | random Gaussian protein embedding | feature-content null |
| `AO1` | ridge solve on **all** of the protein's own cells | oracle ceiling at rank `r` |

Splits: 5-fold over protein closure components. **Primary = KLIFS group** (8
components — the strictest possible kinase split, holding out whole branches of
the kinome). Also reported: KLIFS family (38 components) and single-linkage
pocket-identity ≥0.60 clusters (44 components). Support and test ligands are
disjoint within each held-out protein. 5 seeds. All intervals are **cluster
bootstrap over held-out closure components** (2,000 resamples) — the correct
independence unit.

---

## 5. Literature-backed candidate mechanisms actually tested

Only mechanisms with a stated reason why geometry/sequence alone cannot recover
them were tested, per the instruction not to stack modules.

| candidate | mechanism | why sequence/geometry may not recover it | availability for unseen protein | leakage risk |
|---|---|---|---|---|
| **aligned pocket residue identity** (KLIFS 85-mer) | the residues that physically contact the ligand set steric/H-bond/electrostatic complementarity | it *is* the direct structural determinant — the strongest a-priori candidate | yes, for any kinase | none |
| **pocket physicochemistry** (hydropathy, charge, volume, HB donors/acceptors at each of the 85 aligned positions) | encodes desolvation, electrostatics and steric fit rather than residue identity, so it can generalize across non-identical residues | identity one-hot cannot express "similar chemistry, different residue" | yes | none |
| **conformational-state availability** (KLIFS DFG-in/out, αC-in/out fractions, subpocket occupancy) | type-II / allosteric binding requires an accessible DFG-out or αC-out state; this is a property of the conformational ensemble, **not** of the aligned sequence | a sequence encoder sees one string, not an ensemble | **only 66/82 kinases have any solved structure** | **high** — labels are derived from inhibitor co-crystals, so "has DFG-out structures" partly reflects past medicinal chemistry |
| **ESM-2 t30 full sequence** | the production MetaSieve protein encoder | evolutionary/statistical context beyond the pocket | yes | none |
| **ESM-2 t30 on the 85-mer pocket** | pocket-restricted version of the same | — | yes | none |
| **family / group label** | phylogenetic prior | — | yes | none |
| **homolog-kernel averaging** | non-parametric nearest-homolog transfer | basis-free; avoids any low-rank parameterization artefact | yes | none |

Deliberately **not** pursued: docking-derived energetics, explicit water models,
MD-derived flexibility, and any new neural module. XP1-D1 shows why — see §6.6.

---

## 6. Results

### 6.1 XP1-A1 — variance decomposition (`BLK-METZ-60`, 34,764 cells)

| component | sd (log units) | variance share |
|---|---|---|
| ligand main effect `alpha` | 0.508 | **0.283** |
| protein main effect `beta` | 0.330 | **0.131** |
| interaction + noise | 0.699 | **0.596** |
| total `y` | 0.905 | 1.000 |

`BLK-METZ-70` (704 × 41, 70% uncensored): 0.327 / 0.093 / 0.590 — stable.

The quantity MetaSieve needs is the largest single component of the panel.

### 6.2 XP1-A2 — how much of the interaction is reproducible, and at what rank

Random-cell 5-fold CV on `BLK-METZ-60`, additive baseline vs additive + rank `r`:

| rank | held-out RMSE | `R2_gamma` | per-fold MSE gain [95% CI] |
|---|---|---|---|
| additive | 0.7174 | — | — |
| 1 | 0.6392 | **+0.2060** | +0.1060 [+0.0992, +0.1128] |
| 2 | 0.6236 | +0.2444 | +0.1258 [+0.1207, +0.1308] |
| 3 | 0.6088 | +0.2799 | +0.1440 [+0.1359, +0.1522] |
| 5 | 0.5900 | +0.3236 | +0.1665 [+0.1573, +0.1758] |
| 8 | 0.5755 | +0.3564 | +0.1834 [+0.1772, +0.1896] |
| 12 | 0.5663 | +0.3769 | +0.1940 [+0.1843, +0.2036] |
| 20 | 0.5648 | **+0.3801** | +0.1956 [+0.1853, +0.2060] |
| 30 | 0.5712 | +0.3660 | +0.1884 [+0.1728, +0.2040] |

**≈38% of the interaction residual is reproducible structure, saturating around
rank 12–20 and over-fitting by rank 30, with rank 1 alone delivering 54% of it.**
Implied reproducible interaction sd = `sqrt(0.7174² − 0.5648²) = 0.442` log units
against a residual sd of 0.699 — i.e. the reproducible interaction is comparable
in magnitude to the noise, and larger than the protein main effect (0.330).

### 6.3 XP1-A3 / A4 — is the interaction geometry a stable property of the protein?

The kinase × kinase correlation matrix of interaction-residual columns is
basis-free, so it can be compared between datasets with **no shared compounds**.

| comparison | pairs | Pearson | Spearman | label-permutation null \|r\| (mean / p95) | p |
|---|---|---|---|---|---|
| Metz compound-half 1 vs half 2 (disjoint compounds) | 3,321 | **0.8849** | 0.8688 | 0.0133 / 0.0331 | <5e-4 |
| same, on raw columns (main effects retained) | 3,321 | 0.8570 | — | — | <5e-4 |
| **Metz `pKi` vs Klaeger kinobeads `pKd_app`** (different lab, different technology, different endpoint, **disjoint compounds**) | 1,128 (48 shared kinases) | **0.5650** | 0.4938 | 0.0237 / 0.0589 | <5e-4 |

The interaction geometry replicates across an independent measurement platform.
This is the strongest available evidence that `gamma` is biology and not panel
idiosyncrasy.

### 6.4 XP1-B — transfer to unseen proteins

**Primary: `BLK-METZ-60`, KLIFS group closure, rank 8, k=16, 5 seeds.**
`R2_gamma` is measured against `A2` (additive + support location).

| arm | RMSE | `R2_gamma` vs `A2` [95% CI] |
|---|---|---|
| `A0` population | 0.9224 | −0.3223 [−0.4149, −0.2694] |
| `A1` ligand-only | 0.8571 | −0.1417 [−0.1756, −0.0965] |
| **`A2` additive + support location** | **0.8021** | *(baseline)* |
| **`A4` few-shot section** | **0.7352** | **+0.1600 [+0.1086, +0.1945]** |
| `A6` permuted support | 0.7810 | +0.0521 [−0.0105, +0.1000] |
| `A34` combined | 0.7516 | +0.1221 [+0.0949, +0.1355] |
| **`AO1` oracle rank-8** | **0.6315** | **+0.3802 [+0.3185, +0.4191]** |
| `A3::pocket_identity_kernel` | 0.7940 | +0.0201 [−0.0159, +0.0581] |
| `A3::pocket_onehot` | 0.7949 | +0.0181 [−0.0209, +0.0547] |
| `A3::pocket_physchem` | 0.7913 | +0.0268 [+0.0115, +0.0491] |
| `A3::group_onehot` | 0.7988 | +0.0083 [−0.0066, +0.0292] |
| `A3::family_onehot` | 0.8043 | −0.0055 [−0.0222, +0.0146] |
| `A3::esm2_t30_fullseq` | 0.8017 | +0.0011 [−0.0305, +0.0220] |
| `A3::esm2_t30_pocket85` | 0.8151 | −0.0325 [−0.0709, +0.0185] |
| `A3::klifs_conformation` | 0.7978 | +0.0109 [−0.0017, +0.0301] |
| `A3B::knn_pocket` | 0.8412 | −0.0997 [−0.1725, −0.0515] |
| `A5::pocket_physchem` (wrong protein) | 0.7953 | +0.0170 [+0.0063, +0.0321] |
| `A5::esm2_t30_fullseq` (wrong protein) | 0.8006 | +0.0038 [−0.0210, +0.0244] |
| `A7` random features | 0.9193 | −0.3133 [−0.4580, −0.2282] |

Registered contrasts (paired MSE gain in log units², cluster bootstrap):

| contrast | value [95% CI] | reading |
|---|---|---|
| `Delta_protein = A1 - A4` | **+0.1941 [+0.1131, +0.2743]** | correct protein beats ligand-only |
| `Delta_protein = A1 - A3` | +0.1041 [+0.0559, +0.1410] | but this is almost entirely `A2`'s support intercept |
| `Delta_interaction = A2 - A4` | **+0.1029 [+0.0570, +0.1488]** | support-identified section beats additive+location |
| `Delta_interaction = A2 - A3` | +0.0129 [−0.0103, +0.0340] | **zero-shot adds nothing** |
| `Delta_specific = A6 - A4` | **+0.0695 [+0.0524, +0.0874]** | the section is target-specific |
| `Delta_specific = A5 - A3` | +0.0046 [+0.0010, +0.0088] | 15× smaller; statistically non-zero, experimentally negligible |
| `Delta_oracle = A2 - AO1` | +0.2446 [+0.1635, +0.3203] | headroom that exists and is not being reached |

**Decomposition of `A4`'s gain.** `A6` (correct protein's intercept, *wrong*
protein's section) recovers `R2_gamma = +0.052` — a shared, non-specific residual
ligand correction whose CI includes zero. The genuinely protein-specific part is
`+0.160 − 0.052 = +0.108`.

**Closure level — the central diagnostic.** `R2_gamma` vs `A2`, rank 8, k=16:

| closure | components | `A4` few-shot | `A3` zero-shot pocket kernel | `A3B` homolog kNN | `AO1` oracle |
|---|---|---|---|---|---|
| **group (strict)** | 8 | **+0.160 [+0.109, +0.195]** | **+0.020 [−0.016, +0.058]** | −0.100 [−0.173, −0.051] | +0.380 [+0.318, +0.419] |
| family | 38 | +0.179 [+0.132, +0.228] | **+0.105 [+0.048, +0.158]** | +0.019 [−0.065, +0.110] | +0.402 [+0.350, +0.456] |
| pocket ≥60% identity (permissive) | 44 | +0.172 [+0.119, +0.214] | **+0.084 [+0.051, +0.111]** | +0.037 [−0.050, +0.112] | +0.397 [+0.351, +0.433] |

This is the single most informative table in the study. **`A4` is essentially
invariant to closure level (0.160 → 0.179), while `A3` collapses from 0.105 to
0.020 the moment near-homologs are removed from training.** A genuinely
target-specific mechanism should not care how the split is drawn; a
homolog-interpolation mechanism must. `A4` behaves like the former, `A3` like the
latter.

**Rank (= required `dim` of the biological coordinate), group closure, k=16:**

| rank | `A4` `R2_gamma` | `AO1` `R2_gamma` |
|---|---|---|
| 1 | +0.148 [+0.084, +0.186] | +0.265 [+0.191, +0.311] |
| 2 | +0.157 [+0.089, +0.200] | +0.292 [+0.220, +0.340] |
| 3 | +0.156 [+0.098, +0.195] | +0.314 [+0.250, +0.356] |
| 5 | +0.166 [+0.106, +0.209] | +0.353 [+0.285, +0.400] |
| 8 | +0.160 [+0.109, +0.195] | +0.380 [+0.318, +0.419] |
| 12 | +0.154 [+0.104, +0.188] | +0.408 [+0.353, +0.443] |

**Rank 1 already delivers 93% of the few-shot-attainable signal**, and `A4` is
flat from rank 1 to 12 while the oracle keeps climbing. Higher-rank interaction
structure exists but 16 support observations cannot identify it. The number of
biological coordinates *identifiable at realistic support size* is therefore
**1 to 3** — precisely the regime the frozen sieve can afford (§2.3).

**Support size `k` (identification curve), group closure, rank 8:**

| k | `A4` `R2_gamma` | `AO1` `R2_gamma` | `A4`/`AO1` |
|---|---|---|---|
| 4 | +0.052 [+0.026, +0.066] | +0.479 [+0.434, +0.509] | 11% |
| 8 | +0.092 [+0.068, +0.110] | +0.418 [+0.361, +0.454] | 22% |
| 16 | +0.160 [+0.109, +0.195] | +0.380 [+0.318, +0.419] | 42% |
| 32 | +0.221 [+0.155, +0.265] | +0.361 [+0.297, +0.403] | 61% |
| 64 | +0.258 [+0.186, +0.305] | +0.349 [+0.282, +0.392] | 74% |

A monotone identification curve: the section is *progressively identified* as
support grows, approaching the oracle. `A3` is flat at +0.018…+0.021 across the
entire range — it never identifies anything, at any support size, because it does
not use the support. This is the shape expected of a genuine few-shot
identification problem and not of a capacity or calibration effect.

**Censoring sensitivity (`BLK-METZ-70`, 30% truncated instead of 40%):**
`A4 = +0.123 [+0.076, +0.147]`, `A3 = +0.012 [−0.005, +0.033]`,
`AO1 = +0.354 [+0.281, +0.391]`. The pattern is unchanged on a block with less
truncation and only 41 kinases.

**Secondary panel: `BLK-PDSP-H`** (GPCRs/transporters, radioligand `Ki`,
MMseqs2 40% homology closure, rank 6, k=16):

| arm | RMSE | `R2_gamma` vs `A2` [95% CI] |
|---|---|---|
| `A1` ligand-only | 1.1211 | −0.0595 [−0.1231, −0.0134] |
| `A2` | 1.0892 | *(baseline)* |
| `A4` few-shot section | 1.0762 | **+0.0237 [+0.0074, +0.0520]** |
| `A6` permuted support | 1.0858 | +0.0063 [−0.0063, +0.0222] |
| `AO1` oracle | 1.0086 | **+0.1425 [+0.1074, +0.2015]** |
| `A3::esm2_t30_fullseq` | 1.0912 | −0.0036 [−0.0246, +0.0067] |
| `A3::kmer3_kernel` | 1.0889 | +0.0006 [−0.0031, +0.0045] |
| `A3::family_prefix` | 1.0888 | +0.0007 [−0.0029, +0.0043] |

`Delta_interaction (A2-A4) = +0.0281 [+0.0079, +0.0737]`;
`Delta_specific (A6-A4) = +0.0207 [+0.0102, +0.0424]`;
`Delta_interaction (A2-A3esm) = −0.0043 [−0.0260, +0.0082]`.

**The same qualitative pattern replicates in a different protein class with a
different assay technology.** The effect size is much smaller, consistent with
the 0.714-log-unit per-report noise measured in §3.4, and `A4`'s `R2_gamma` of
0.024 is **below the preregistered 0.05 negligibility threshold** on this panel.

### 6.5 XP1-E — destructive control for left-censoring

Synthetic panels with the **real** `mu, alpha, beta`, the **real** residual sd
(0.699), and the **same 40% truncation fraction**, then run through the identical
XP1-B pipeline. The planted interaction sd of 0.456 was fixed from the rank curve
before §6.2 was recomputed under the two-stage fit, which puts the reproducible
interaction sd at 0.442 — a 3% difference that does not affect any reading below.

| scenario | `A4` `R2_gamma` | `Delta_specific (A6-A4)` | `AO1` `R2_gamma` | `A3::pocket` |
|---|---|---|---|---|
| **additive only (no interaction)** | −0.033 | **+0.00004 [−0.0005, +0.0011]** | +0.074 | −0.0001 |
| planted rank-1, sd 0.456 | +0.087 | +0.032 | +0.264 | +0.0007 |
| planted rank-3, sd 0.456 | +0.052 | +0.024 | +0.261 | +0.0003 |
| **real `BLK-METZ-60`** | **+0.160** | **+0.069** | **+0.380** | +0.020 |

Readings:

1. **Falsification passed.** Truncation alone produces `Delta_specific = +0.00004`
   with a CI straddling zero — it cannot manufacture the observed `+0.0695`.
2. **Positive control passed.** The pipeline detects planted interaction and its
   specificity, so the near-zero zero-shot arms are meaningful measurements, not
   an inert pipeline.
3. **An honest deduction.** Truncation *does* inflate the oracle by about
   `+0.074`. The real oracle `+0.380` should therefore be read as roughly
   `+0.31` of genuine interaction headroom.
4. The real effect exceeds even the matched planted interaction, because real
   `gamma` is spectrally concentrated (rank 1–3) while the planted one is
   isotropic.

### 6.6 XP1-D — can any protein statistic reach the coordinate?

**D1. Direct predictability of the interaction coordinate `V`** (rank 8, whitened),
leave-one-closure-component-out `R2`, ridge with inner leave-one-component-out
`lambda` selection. This removes the low-rank parameterization and the affinity
readout from the question and asks only: *is `V` a function of the protein
representation?*

| representation | group | family | pocket ≥60% | leave-one-protein-out |
|---|---|---|---|---|
| pocket identity kernel (KLIFS 85-mer) | **−0.087** | +0.140 | +0.133 | **+0.516** |
| pocket one-hot | −0.100 | +0.148 | +0.135 | +0.509 |
| pocket physicochemistry | −0.105 | +0.121 | +0.093 | +0.516 |
| kinase group one-hot | −0.131 | +0.073 | +0.082 | +0.229 |
| kinase family one-hot | −0.147 | −0.052 | +0.012 | +0.463 |
| **ESM-2 t30 full sequence** (production encoder) | −0.134 | −0.059 | −0.050 | +0.293 |
| ESM-2 t30 on the 85-mer pocket | −0.136 | −0.064 | −0.040 | +0.316 |
| KLIFS conformational state | −0.147 | +0.008 | −0.016 | +0.081 |
| **random Gaussian (null)** | **−0.152** | −0.058 | −0.063 | −0.025 |

Three readings, each load-bearing:

1. **At group closure every representation is statistically indistinguishable
   from the random-feature null** (−0.087 … −0.152 vs −0.152). Not one of them
   carries cross-branch information about the interaction coordinate.
2. **At leave-one-protein-out the aligned pocket reaches `R2 = 0.52`.** The
   information is emphatically there — it is homolog-interpolation information.
   The gap between +0.516 and −0.087 *is* the difference between partner
   compatibility and target-specific affinity direction, measured on one axis.
3. **The structure-aligned 85-residue pocket beats the production ESM-2 encoder
   by a wide margin** (+0.516 vs +0.293) even though neither generalizes across
   groups. If MetaSieve keeps a zero-shot protein pathway at all, the aligned
   pocket is the better representation — but §6.4 shows that pathway is worth
   `R2_gamma = 0.020` at strict closure either way.

**D2. Learning curve** — is the zero-shot failure a sample-size failure? Ridge
map from protein features to `V`, single-group holdout, training proteins
subsampled:

| representation | n=20 | n=35 | n=50 | n=65 |
|---|---|---|---|---|
| pocket identity kernel | +0.001 | +0.004 | +0.007 | +0.004 |
| ESM-2 t30 full sequence | −0.208 | −0.395 | −0.420 | −0.278 |
| KLIFS conformational state | −0.001 | −0.045 | −0.022 | −0.033 |

Flat at zero across a 3.25× increase in training proteins, while the *same*
features on the *same* 82 proteins reach `R2 = 0.52` when homologs are allowed.
This bounds — it does not logically exclude — a large-`n` effect: it says the
cross-group map is not being learned slowly, it is not being learned at all in
the 20–65 protein regime, so a kinome-scale panel would have to change the
*qualitative* behaviour, not merely extend a rising curve.

**D3. What is `V`?** In-sample explained variance of the interaction coordinate:

| predictor | in-sample `R2` | `R2` expected by chance at this df |
|---|---|---|
| kinase group (8 levels) | 0.351 | 0.086 |
| kinase family (36 levels) | 0.806 | 0.432 |
| best single KLIFS pocket position (pos 42) | 0.366 | 0.232 |
| mean over all 85 pocket positions | 0.166 | 0.232 |

`V` carries real phylogenetic structure (group 0.351 vs 0.086 chance) that
nevertheless does not extrapolate across groups — the signature of within-branch
structure. Only six of 85 pocket positions exceed the chance level: 42, 54, 50,
20, 22, 33. Positions 42/50/51/54 flank the **gatekeeper (KLIFS position 45)** and
**hinge (46–48)**; 18/20/22 lie in the **β3–αC** region. (Numbering verified
against CDK2/ABL1/EGFR/SRC, where DFG sits at 81–83 and the CDK2 gatekeeper Phe80
maps to position 45.) These are the textbook determinants of ATP-site
selectivity, which is reassuring — but the association is **in-sample, univariate
and barely above chance**, and it does not survive into any out-of-group
prediction. It is a sanity check, not a mechanism.

---

## 7. Failure analysis — which hypothesis each result falsifies

| falsified claim | falsifying evidence |
|---|---|
| `H2`: the data cannot identify protein-specific interaction | XP1-A1/A2/A3/A4 and `AO1 = +0.38`: `gamma` is 59.6% of variance, 38% of it reproducible, replicates across platforms, and transfers to unseen kinase *groups* when the coordinate is known |
| the apparent interaction is a censoring artefact | XP1-E additive-only control: `Delta_specific = +0.00004` |
| the apparent interaction is a protein-ID / potency shortcut | `A2` already contains the support-estimated protein location; `A4` beats it by `+0.103` |
| the apparent interaction is a generic (non-specific) ligand correction | `A6` permuted support recovers only `+0.052` of `A4`'s `+0.160`, and `Delta_specific` CI excludes zero |
| the apparent interaction is an assay artefact | `BLK-METZ-60` is one laboratory and one assay family; assay identity is constant by construction. Replicated on PDSP where assay differs entirely |
| **protein sequence/pocket features carry target-specific affinity direction** | XP1-B `A3` ≈ 0 for every representation at group closure, and `A5` (wrong protein) matches `A3` |
| ESM-2 (the production encoder) carries it | `A3::esm2_t30_fullseq` `R2_gamma = +0.001 [−0.031, +0.022]` (Metz), `−0.004 [−0.025, +0.007]` (PDSP) |
| the zero-shot failure is merely a small-sample artefact | XP1-D2 learning curve flat at +0.001…+0.007 over n=20→65, while the same features reach `R2 = 0.52` leave-one-protein-out on the same 82 proteins |
| the zero-shot failure is an artefact of the low-rank parameterization | XP1-D1 predicts `V` directly, with no affinity readout, and `A3B::knn` bypasses the basis entirely; both fail at group closure |
| conformational-state availability rescues it | `A3::klifs_conformation` `R2_gamma = +0.011 [−0.002, +0.030]`; additionally **unavailable for 16/82 kinases** and carrying a high circularity risk, so it fails the "available for unseen proteins" criterion regardless |
| homolog averaging rescues it | `A3B::knn_pocket` is **worse** than `A2` at group closure (`−0.0997`) |

Two mechanisms are therefore **closed** by this work: the zero-shot pocket-feature
surface and the homolog-kernel surface. Conformational-state availability is
**not closed on the merits** but is disqualified on availability and circularity;
retesting it would require apo-structure or prediction-derived ensembles, and
that is not justified until the ligand-side question in §10 is answered.

---

### 7.1 Limitations that bound these claims

1. **Unseen ligands are not tested.** `P-METZ` publishes compound indices, not
   structures, so `u(L)` is a lookup table over ligands shared between training
   and held-out proteins. XP1 establishes transfer to unseen *proteins* only.
   This is the subject of `XP2` (§10) and is the largest open risk.
2. **`U` is panel-internal.** Within each fold `U` is fitted on training proteins
   only — there is no leakage into held-out proteins — but it is fitted on Metz
   kinases. §2.5 requires `U` frozen from a disjoint corpus before meta-training;
   whether `U` survives that transfer is untested.
3. **Eight independence units.** Group-closure intervals are cluster bootstraps
   over 8 kinase groups. They are the correct unit, but 8 is small; the family
   (38) and pocket (44) closures agree, which is the main reassurance.
4. **Kinase-dominated.** The primary panel is one protein family. PDSP replicates
   the direction in a different class but at `R2_gamma = 0.024`, below the
   preregistered negligibility floor.
5. **Censoring conditions the estimand.** All cells at the `pKi = 4.0` floor are
   excluded, so the estimand is the interaction *among measurable affinities*.
   XP1-E shows truncation does not manufacture the effect; it does not make the
   excluded region observable.
6. **PDSP assay caveat.** In PDSP the radioligand is nearly a function of the
   target, so a radioligand-specific competition effect would be absorbed into
   `beta_j`; but a radioligand-by-compound interaction would not be, and could
   contribute to PDSP's `gamma`. `BLK-METZ-60` is immune (one assay family).
7. **`A2`'s intercept.** `A2` uses the intercept-only support mean while `A4`
   uses the joint ridge intercept. `A6` (joint intercept, *wrong* section) bounds
   the entire intercept-plus-shared-component contribution at `+0.052`; the
   `A4 − A6 = +0.108` specificity is unaffected by this choice.
8. **Gate thresholds are self-set.** `R2_gamma >= 0.05`, CI lower `> 0.02` were
   fixed in the preregistration before any arm was scored, but they are a
   judgement, not an external standard.
9. **Klaeger is used only binarised** (93.3% floor), so §6.3's cross-platform
   `r = 0.565` compares a continuous geometry against a hit-pattern geometry and
   is, if anything, an underestimate.

---

## 8. Biological interpretation

The signal `A4` recovers is **true ligand–protein interaction**, not any of the
alternatives, and the arm structure separates them explicitly:

- not **ligand main effect**: `alpha_i` is fitted on training proteins and given
  identically to every arm; `A1` is 0.086 log units *worse* than `A2`;
- not **protein main effect**: `A2` already carries the support-estimated protein
  location and `A4` beats it;
- not **raw affinity**: the metric is RMSE on the residual after both main
  effects;
- not **compatibility**: compatibility is what `A3`/`A3B` measure — the ability to
  tell that a protein and ligand belong together based on protein features — and
  those arms carry the signal **only when a homolog is in training**. That is the
  quantitative content of the project's existing
  `PARTNER_COMPATIBILITY_PARTIALLY_IDENTIFIED` state;
- not **geometry**: pocket residue identity *is* the geometric determinant, and it
  fails at group closure. This directly confirms the project's framing that a
  validated geometry bridge is not an affinity bridge.

Biologically, the recovered coordinate is close to what medicinal chemists call a
kinase's **selectivity phenotype**: a low-dimensional description of which
chemotypes it accepts. XP1-D3 shows it is only weakly aligned with phylogeny, and
the cross-platform replication (§6.3) shows it is a property of the protein, not
of a panel. The reason sequence cannot reach it across kinase groups is that the
same phenotype is realised by different residue combinations in different
branches — an *epistatic*, not additive, function of the pocket — and 65 training
proteins cannot identify such a function.

---

## 9. Mathematical interpretation and the admission verdict

For the two candidate statistics, evaluated against the registered Gate (§9 of
the preregistration), on the strict group closure:

| Gate condition | `z_feature = phi(protein)` (zero-shot) | `z_section = <u_Q, vhat(S)>` (support-identified) |
|---|---|---|
| 1. `R2_gamma >= 0.05`, CI lower `> 0.02` | **FAIL** (best +0.027, lower +0.012) | **PASS** (+0.160, CI [+0.109, +0.195]) |
| 2. `Delta_interaction` CI `> 0` | **FAIL** (+0.013 [−0.010, +0.034]) | **PASS** (+0.103 [+0.057, +0.149]) |
| 3. `Delta_specific` CI `> 0` | technically pass, **experimentally negligible** (+0.0046) | **PASS** (+0.0695 [+0.0524, +0.0874]) |
| 4. `Delta_protein` CI `> 0` | pass, but attributable to `A2`'s intercept | **PASS** (+0.194 [+0.113, +0.274]) |
| 5. survives the feature-content null | pass | n/a; the analogue `A6` is passed |
| 6. dimension compatible with §2.3 | n/a | **PASS** (rank 1–3) |

Classifying each candidate on the ladder the project asked for:

- **`z_feature` (ESM-2, pocket sequence, pocket physicochemistry, conformation,
  family):** *merely correlated*, and only within homology neighbourhoods. Not
  informative about `gamma` across kinase groups. **Not admissible.**
- **`z_section = <u_Q, vhat(S_j)>`:** *informative* and *identifiable* — it is
  recovered from `k=16` labelled support observations, it is target-specific
  under derangement, it transfers across the strictest kinase split, and it
  survives the truncation falsification. It is **not sufficient**: the oracle
  retains `+0.38` (≈`+0.31` after the truncation correction) versus `A4`'s
  `+0.16`, so roughly half of the reproducible interaction remains unreached at
  `k=16`. It satisfies every structural requirement of §2.1–2.4: computable from
  `(S, Q, gamma)`, bounded after clipping, permutation-invariant in `S`,
  query-label-free, and **one-dimensional as it enters `z`**.

**Verdict: `z_section` passes the registered admission Gate on the primary panel
and replicates directionally on an independent panel, but is NOT promoted to
production `z`.** Per `task.md`, production admission additionally requires a
separately registered independent-source replication and a sealed novel-target
transfer Gate. XP1's PDSP replication is directionally correct but its effect
size (`R2_gamma = 0.024`) is below the preregistered negligibility floor, so the
honest status is:

```text
SUPPORT_IDENTIFIED_INTERACTION_SECTION_ADMITTED_ON_PRIMARY_PANEL
INDEPENDENT_REPLICATION_DIRECTIONAL_ONLY_EFFECT_BELOW_FLOOR
BIOLOGICAL_Z_NOT_YET_ADMITTED_TO_PRODUCTION
```

### 9.1 Concrete interface recommendation (no change to the frozen operator)

`model/component_statistic.py` currently implements

```
prediction(P, L; S) = biological_surface(P, L) + location(S)
```

with `location(S)` a bounded **scalar** — algebraically the **rank-0** case of
what XP1 identifies. XP1 says two things about this decomposition, and both are
measurable:

1. `biological_surface(P, L)` — the zero-shot protein-dependent term — contributes
   `R2_gamma ≈ 0` at strict closure. This is *why* the registered F6I total
   verdict was `NOT_ADMISSIBLE`; the present work supplies the mechanism rather
   than only the verdict.
2. The scalar `location(S)` should be generalised to a **rank-`r` section**,
   `r in {1,2,3}`:

```
prediction(P, L; S) = mu + alpha(L) + b(S) + <u(L), v(S)>
```

where `u(L)` is a frozen ligand-loading vector belonging to the deployment
declaration `gamma` (see §2.5), and `(b, v)` are the ridge solution on the
support residuals. The single admissible new `z` coordinate is the bounded scalar
`<u(L), v(S)>`, which needs **one** CSMO view — no change to `Z`'s cube
structure, `B(z)`, `kappa`, `Delta_m`, the ridge `mu`, the mesh `h`, or the
operator `A(F,z) = K(B(z)F(z))`.

---

## 10. Final decision and the single next experiment

```text
FINAL DECISION: OBJECTIVE_OR_PARAMETERIZATION_FAILURE
  secondary:    REPRESENTATION_FAILURE  (confined to the zero-shot protein pathway)
  rejected:     DATA_IDENTIFIABILITY_FAILURE
```

The necessary signal is present in accessible data, is reproducible across
platforms, and is recoverable from the biological information available at
inference — but only through a **support-identified low-dimensional section**,
not through the feature-conditioned protein surface the current interface
assumes. Redesign the interface between the biological statistic and the frozen
operator; do not enlarge the network, and do not modify CSMO, Band, the geometry
bridge, or the mathematics.

**The single most decision-relevant untested question** is the mirror image of
what XP1 answered. XP1 held out *proteins* and reused *ligands*, so `u(L)` was a
lookup table. Production needs `u(L)` for unseen ligands, from chemistry alone.
If the ligand loading is no more predictable from structure than the protein
coordinate was from sequence, the section is identifiable but not deployable, and
the correct conclusion flips back toward representation failure — on the ligand
side. That experiment (`XP2`) is cheap: PDSP and Klaeger both carry structures,
and the arm structure is the transpose of XP1-B.

**Not authorized by this report:** DAVIS access, recipient labels, ChEMBL
`X1`/`X2`, `model/` promotion, changes to the frozen operator, P2–P4, or any
end-to-end DTA claim.

---

## 11. Reproduction

```bash
D:/anaconda/envs/drug/python.exe research/crossed_panel_identification/xp1a_existence.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_identification/xp1b_sweep.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_identification/xp1c_pdsp.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_identification/xp1d_statistic.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_identification/xp1e_truncation_control.py
```

Acquisition (run once): `dl_pdsp.py`-equivalent PDSP export, `dl_kinase.py`,
`dl_klifs.py`, `build_protein_features.py`, `build_conformation_features.py`,
`pdsp_build.py`. All release SHA-256 values are asserted at load time by
`panels.verify_releases()`; every reported run verified them.

Artifacts in `report/crossed_panel_identification/`:

| file | content |
|---|---|
| `xp1a_existence.json` | variance decomposition, rank curve, geometry replication, PDSP noise ceiling |
| `xp1b_sweeps.json` | every arm of every closure / rank / support-size / density configuration |
| `xp1c_pdsp.json` | PDSP replication |
| `xp1d_statistic.json` | direct `V` predictability, learning curve, `V` structure |
| `xp1e_truncation_control.json` | additive-only and three planted-interaction synthetic panels |
| `xp1*_console.txt` | verbatim console output of each run |

`panels.verify_releases()` asserts the three frozen SHA-256 values at load time
and was green in every run reported here. Raw panels are not redistributed
through Git, per `DATA_AVAILABILITY.md`; `dataset/raw/crossed_panels/*/acquisition_manifest.json`
records URL, byte count and SHA-256 for each.
