# XP2 — Deployability Of The Crossed-Panel Interaction Section

Date: 2026-08-08
Preregistration: `research/crossed_panel_deployability/PREREG_XP2.md` (written
before any XP2 arm was scored)
Environment: `drug` conda env, Python 3.11.15, numpy 1.26.4, scipy 1.17.1,
pandas 2.3.3, rdkit 2023.09.6, torch 2.6.0+cu124 (CUDA on RTX 4060 Laptop, used
only for the frozen ESM-2 and ChemBERTa encoders), transformers 4.46.3.
Every inference reported here is deterministic linear algebra.

Label-read counters, asserted in every artifact: DAVIS `0`, recipient `0`,
ChEMBL37 affinity `0`, PKIS2 `0`, Anastassiadis `0`.

---

## 0. Terminal verdict

```text
XP2 TERMINAL VERDICT
  CROSSED_INTERACTION_REPRODUCED
  K_LE_5_SECTION_NOT_IDENTIFIED
  PANEL_LOCAL_LOW_RANK_META_LEARNING
  BIOLOGICAL_LANDING_NOT_IDENTIFIED
  EXTERNAL_REPLICATION_FAILED
```

Not claimed, and explicitly refused: `DEPLOYABLE_SECTION_STATISTIC_IDENTIFIED`,
`DOUBLE_HELD_OUT_SECTION_IDENTIFIED`. XP1's evidence reproduced cleanly, and the
ligand side turned out *better* than feared — but the section fails the frozen
non-negligibility floor at `k <= 5`, loses target specificity entirely once both
closures are enforced together, gains nothing from any protein representation,
and does not replicate on an independent platform.

**This is not a deep biological meta-learning model and must not be described as
one.** Of the five layers in §1, two are established, three are not.

The one bounded verdict that does **not** apply is
`LIGAND_SIDE_DEPLOYMENT_REPRESENTATION_FAILED`. The ligand loading is
transferable to unseen scaffolds. That was the defect XP2 was principally
designed to test, and it is resolved in the affirmative.

---

## 1. What each layer established

The five layers the registration asked to be kept apart:

| layer | status | decisive evidence |
|---|---|---|
| **interaction existence** | **ESTABLISHED** | XP1 reproduced from immutable artifacts (XP2-A, 18/18 checks); interaction is 59.6% of panel variance, geometry replicates across disjoint compound halves (`r = 0.885`) and across an independent platform (`r = 0.565`) |
| **support identifiability** | **ESTABLISHED ONLY UNDER LIGAND REUSE** | with ligands reused, the section is target-specific from `k = 2` upward (every specificity control CI-positive); the identifiable dimension is exactly `min(k-1, d)`, so `k <= 5` caps it at 4 and `k = 1` gives **zero**; magnitude stays below the frozen floor at every `k <= 5` |
| **ligand loading recoverability** | **OBSERVED** | gauge-invariant loading transfer to unseen scaffolds, ECFP `R2 = +0.199 [+0.133, +0.261]` at `d = 3` against random-feature `+0.025` and mean-loading `+0.024` |
| **joint deployability** (both closures at once) | **NOT ESTABLISHED** | under simultaneous protein-group and ligand-scaffold closure the derangement, permutation, zero-adaptation and random-correction controls **all** have intervals containing zero, and random ligand features reproduce the entire gain (§6.1–6.2) |
| **protein-side biological landing** | **NOT IDENTIFIED** | every zero-shot protein representation sits within `0.002` of a random protein embedding (§7) |
| **external replication** | **FAILED** | direction transfer to Klaeger gives `Delta_interaction = +0.00346 [−0.00234, +0.00863]` and `Delta_deploy = −0.0766` (§6.3) |
| **probability-law admission** | **NOT EVALUATED** | nothing in XP2 scores `K(B(z)F(z))`; interaction `R2` is not law calibration and is not offered as such |

---

## 2. XP2-A — evidence and artifact audit

`XP1_REPRODUCTION_AUDIT.json`, verdict **`XP1_EVIDENCE_REPRODUCED`**, 0 failures.

Release integrity: all five pinned SHA-256 values match. Recomputed from the
immutable releases: `BLK-METZ-60` shape `704 x 82` / `34,764` cells; variance
shares `0.2831 / 0.1309 / 0.5962` against reported `0.283 / 0.131 / 0.596`;
split-half interaction geometry `r = 0.8849` against reported `0.8849`; an
independently recomputed rank-3 fold-0 gain of `0.14059` falls inside the
archived CI `[0.13593, 0.15215]`.

Construction facts were read out of the XP1 source itself, not from prose:
ligands reused across train and test proteins **confirmed**; no ligand-scaffold
closure in XP1 **confirmed**; `U` fitted on training proteins only **confirmed**;
support/test disjoint within protein **confirmed**; `k_support` default 16
**confirmed**.

### 2.1 A correction to XP1 that the audit produced

XP1 described `BLK-METZ-60` as left-censored at a single `pKi = 4.0` floor. The
journal supplement in fact encodes **three** cell states: 103,118 measured
values, 154,175 left-censored strings at **50 distinct thresholds** spanning
`4.0`–`6.2`, and 405,482 untested blanks. The mirror CSV collapses the latter two
to `4.0`.

The audit shows XP1's analysis set was nevertheless **exactly right**: its mask
admitted 49,457 cells, of which 49,457 are genuinely measured — zero censored,
zero untested — and those values match the supplement with max `|diff| = 0.0`.
Only the *censoring model* used in XP1's destructive control was an
approximation (single nominal floor instead of 50 per-cell thresholds). No XP1
conclusion changes; the description does.

Two initial audit failures were defects in the audit harness itself — a
positional column slice that mislocated the kinase block after `set_index`, and a
substring scan that flagged the audit file's own exclusion assertions. Both were
replaced with **stricter** checks (explicit metadata-column naming; a
data-access-call detector rather than a substring detector) and rerun.

---

## 3. The XP2 panel

XP2 promotes `metz.xls` — the journal supplement — from provenance-only to
primary source, because it carries `Canonical_Smiles`, per-cell censoring
notation and 3,858 compounds that the derived matrix discards.

`BLK-METZ-XP2`: **928 compounds x 147 kinases, 32,849 measured cells, density
0.241, 518 Bemis-Murcko scaffolds merged into 258 scaffold components at ECFP4
Tanimoto >= 0.5, 8 KLIFS groups.** Index SHA-256
`7bcb2c05daa4aa5a…`, panel SHA-256 `beb97a0e125ccabf…`.

This is the panel XP1 could not have: it supports **simultaneous** protein-group
and ligand-scaffold closure, which is what the whole stage turns on.

---

## 4. XP2-B — ligand landing (the defect XP2 was built to test)

`LIGAND_LANDING_AUDIT.json`. Scaffold-component folds; the interaction basis and
main effects are fitted on training ligands only; a held-out ligand's target
parameters are its own least-squares `(alpha, u)` against that frozen basis; the
predictor sees chemistry only.

Gauge-invariant interaction-reconstruction `R2` on unseen scaffolds:

| arm | d=1 | d=2 | d=3 | d=5 |
|---|---|---|---|---|
| **L-ECFP** (Morgan r=2, 1024 bit) | +0.205 [−0.029, +0.292] | +0.144 [+0.036, +0.260] | **+0.199 [+0.133, +0.261]** | +0.199 [+0.168, +0.247] |
| L-CHEMBERTA (frozen encoder) | +0.147 | +0.077 | +0.107 [+0.075, +0.171] | +0.079 |
| L-DESC (10 descriptors) | −0.028 | +0.019 | +0.053 [+0.002, +0.102] | +0.039 |
| L-RANDOM (control) | +0.008 | +0.003 | +0.025 | +0.017 |
| L-MEAN (control) | +0.008 | +0.003 | +0.024 | +0.017 |

At `d = 3`: sign agreement ECFP `0.644` vs random `0.586`; within-ligand Spearman
`0.351` vs `0.193`. The ligand main effect `alpha` is predictable from chemistry
at `R2 = +0.206` (ECFP) against `−0.003` for random features.

**Conclusion: `LIGAND_LOADING_RECOVERABILITY_OBSERVED`.** `u(L)` is not a pure lookup table. This is a recoverability observation, not a deployability claim: XP2-D shows the recovered loading does **not** make the section deployable. A frozen fingerprint plus a linear
ridge recovers a fifth of the interaction energy for compounds whose scaffolds
were never seen. The smallest arm (10 descriptors) is much weaker, and the
learned chemical encoder is weaker than plain ECFP — capacity was not the binding
constraint, chemical resolution was.

### 4.1 The gauge result, stated plainly

Coordinate-wise loading `R2` in a fold-local gauge is ≈0 or negative in the same
runs where the gauge-invariant reconstruction `R2` is clearly positive (e.g.
`+0.100` coordinate-wise vs `+0.199` gauge-invariant at `d = 3`). Latent
coordinates are not stable objects; only gauge-invariant functionals are. This is
the empirical form of defect #6 and it is why §8 forbids naming latent factors.

---

## 5. XP2-C — true few-shot audit at `k <= 5`

`K5_SECTION_AUDIT.json`. Protein-group closure, ligands reused (the
XP1-comparable regime), `d = 3`, five seeds, cluster-bootstrap intervals taken as
the **wider** of protein-component and scaffold-component clusterings.

| `k` | identified dim | median cond | query coverage | SEC `R2_gamma` vs ADD | vs foreign support | vs permuted support | `Delta_deploy` vs ligand-only |
|---|---|---|---|---|---|---|---|
| 1 | **0.00 / 3** | ∞ | 0.000 | **+0.0000** | +0.0000 | +0.0000 | **−0.431 [−0.577, −0.287]** |
| 2 | 1.00 / 3 | 1.0 | 0.343 | +0.0088 [+0.0016, +0.0121] | +0.0030 [+0.0006, +0.0046] | — | **−0.115 [−0.247, −0.030]** |
| 3 | 2.00 / 3 | 2.8 | 0.672 | +0.0187 [+0.0087, +0.0242] | +0.0040 [+0.0006, +0.0069] | — | −0.026 [−0.156, +0.086] |
| 4 | 3.00 / 3 | 9.2 | 1.000 | +0.0230 [+0.0144, +0.0275] | +0.0068 [+0.0023, +0.0108] | — | +0.040 [−0.089, +0.142] |
| 5 | 3.00 / 3 | 4.8 | 1.000 | **+0.0248 [+0.0114, +0.0321]** | +0.0084 [+0.0016, +0.0128] | +0.0202 [+0.0038, +0.0297] | +0.078 [−0.054, +0.188] |

Rank sweep at `k = 5`: `d=1` `+0.0193`, `d=2` `+0.0242`, `d=3` `+0.0248` —
saturating immediately, reproducing XP1's rank finding on a different panel with
a different construction. Oracle at `k = 5` is `+0.264`.

Three readings:

1. **The identifiability ledger is exact.** The unpenalised support intercept
   absorbs the mean of the support loadings, so the identifiable section
   dimension is `rank(U_S − mean U_S) = min(k−1, d)`. Measured: `0, 1, 2, 3, 3`
   for `k = 1…5`. At `k = 1` the ridge returns `v = 0` and the section arm is
   **identically** the additive arm — `R2_gamma` is `+0.0000`, not
   approximately zero. Any nonzero prediction there would have been the ridge
   prior speaking, not the support.
2. **The section is real and target-specific.** Every specificity control is
   CI-positive from `k = 2` upward: foreign support, permuted support, zero
   adaptation and a norm-matched random correction all lose to the correct
   support.
3. **And it is too small to deploy.** `R2_gamma` maxes at `+0.0248`, half the
   preregistered `0.05` floor, and `Delta_deploy` — value over ligand-only
   chemistry — never clears zero at any `k <= 5`. At `k <= 2` it is decisively
   *negative*: a support intercept estimated from one or two observations injects
   more variance than the section removes.

The threshold was not moved. `K_LE_5_SECTION_NOT_IDENTIFIED`.

---

## 6. XP2-D — double held-out

`DOUBLE_HELD_OUT_RESULT.json`. Protein groups **and** ligand scaffold components
held out simultaneously. Test cells are `G_f x S_f`; training is
`(not G_f) x (not S_f)`; support for a test protein comes from that protein's
cells on **training scaffolds only**; the query ligand's `alpha` and `u` come
from chemistry. Nothing from `S_f` enters training, feature fitting,
hyperparameter selection or support. The test set was opened once and no model,
threshold, rank, feature family or support-selection rule was changed after.

Support-size ladder at `d = 3`, `L-ECFP`:

| `k` | identified dim | coverage | SEC `R2_gamma` vs ADD | vs **foreign support** | `Delta_deploy` vs ligand-only |
|---|---|---|---|---|---|
| 1 | 0.00 / 3 | 0.000 | +0.0000 | +0.0000 | **−0.534 [−0.670, −0.437]** |
| 2 | 1.00 / 3 | 0.338 | +0.0175 [+0.0104, +0.0219] | +0.0009 [−0.0031, +0.0023] | **−0.181 [−0.284, −0.129]** |
| 3 | 2.00 / 3 | 0.667 | +0.0268 [+0.0105, +0.0390] | −0.00001 [−0.0019, +0.0019] | **−0.095 [−0.174, −0.046]** |
| 4 | 2.99 / 3 | 0.998 | +0.0219 [+0.0085, +0.0315] | +0.0023 [−0.0026, +0.0050] | **−0.060 [−0.130, −0.009]** |
| **5 (primary)** | **3.00 / 3** | **1.000** | **+0.0199 [+0.0076, +0.0283]** | **+0.0019 [−0.0048, +0.0055]** | **−0.033 [−0.091, +0.023]** |

Rank sweep at `k = 5`:

| `d` | SEC `R2_gamma` | vs foreign support | `Delta_deploy` |
|---|---|---|---|
| 1 | +0.0060 [−0.0010, +0.0101] | −0.0009 [−0.0057, +0.0016] | −0.046 [−0.099, +0.017] |
| 2 | +0.0165 [+0.0074, +0.0219] | +0.0008 [−0.0058, +0.0043] | −0.036 [−0.093, +0.024] |
| 3 | +0.0199 [+0.0076, +0.0283] | +0.0019 [−0.0048, +0.0055] | −0.033 [−0.091, +0.023] |
| 5 | +0.0267 [+0.0090, +0.0391] | +0.0034 [−0.0015, +0.0055] | −0.027 [−0.089, +0.027] |

**Three findings, in order of importance.**

1. **Target specificity disappears.** Under protein-only closure the
   foreign-support control was beaten with a CI clear of zero from `k = 2`
   upward. Under the double closure it is not beaten at **any** `k` or **any**
   `d <= 5`: every derangement interval straddles zero. Using a foreign
   protein's section is statistically indistinguishable from using the correct
   one. Whatever `SEC` gains over `ADD` is a generic correction that any
   protein's `v` supplies equally well.
2. **This is not an identifiability artefact.** At `k = 5, d = 3` the support
   design has rank `3.00/3` and query coverage `1.000` — the section is fully
   determined numerically. The failure is substantive.
3. **`ADD` is the wrong yardstick here, and saying so matters.** Under the double
   closure a support intercept estimated from `k <= 5` cells on an unseen protein
   injects more variance than it removes, so plain ligand-only chemistry beats
   the additive baseline outright (`LIG` `R2_gamma` vs `ADD` = `+0.381` at
   `k = 1`, `+0.187` at `k = 2`, `+0.123` at `k = 3`). `R2_gamma` vs `ADD` can
   therefore look healthy purely because its denominator is weak. The honest
   deployability number is the preregistered `Delta_deploy` against ligand-only,
   and it is negative at every `k <= 4` with a CI clear of zero, and negative
   with a CI spanning zero at `k = 5`.

### 6.1 Every arm at the registered primary configuration

`closure = double, d = 3, k = 5, L-ECFP`, 705 tasks, five seeds:

| arm | RMSE | `R2_gamma` vs ADD |
|---|---|---|
| `P0` population mean only | 0.9461 | +0.0271 [−0.1609, +0.0945] |
| **`LIG` ligand-only chemistry** | **0.9320** | **+0.0560 [−0.0111, +0.1016]** |
| `ADD` additive + support intercept | 0.9592 | *(baseline)* |
| `ZERO` joint intercept, `v = 0` | 0.9506 | +0.0180 [+0.0071, +0.0261] |
| **`SEC` the section** | **0.9496** | **+0.0199 [+0.0076, +0.0283]** |
| `SEC-UHATSUP` strict, support `u` also predicted | 0.9548 | +0.0093 [+0.0057, +0.0151] |
| `PERM` permuted support labels | 0.9503 | +0.0185 [+0.0039, +0.0289] |
| `FOREIGN` foreign protein's section | 0.9506 | +0.0179 [+0.0076, +0.0223] |
| `RANDCORR` norm-matched random correction | 0.9514 | +0.0163 [+0.0053, +0.0233] |
| `ORACLE-TRSC` all training-scaffold cells | 0.8906 | +0.1380 [+0.1176, +0.1962] |

Read the middle block: `SEC`, `ZERO`, `PERM`, `FOREIGN` and `RANDCORR` lie within
**0.002 RMSE of one another**. Permuting the support labels, substituting a
foreign protein's section, or replacing the interaction with a norm-matched
random vector all cost essentially nothing. The registered contrasts confirm it:

```
Delta_specific_zero    = +0.00181 [-0.00106, +0.00528]
Delta_specific_foreign = +0.00185 [-0.00477, +0.00552]
Delta_specific_perm    = +0.00135 [-0.00127, +0.00713]
Delta_randcorr         = +0.00331 [-0.00028, +0.00652]
Delta_deploy (vs LIG)  = -0.03314 [-0.09108, +0.02326]
Delta_oracle           = +0.12695 [+0.10567, +0.16818]
```

Two further observations that a reader should not miss. **The population mean
alone (0.9461) beats every support-adapted arm**, and ligand-only chemistry
(0.9320) beats all of them by a clear margin — at `k = 5` the support-based
adaptation is not merely uninformative, it is a net cost. And the oracle at
`+0.1380` shows the double-held-out cells *do* contain recoverable interaction
structure; it is simply unreachable from five support measurements.

### 6.2 The random-feature control settles it

Ligand feature family, double closure, `k = 5`, `d = 3`:

| ligand features | SEC `R2_gamma` vs ADD | vs foreign support |
|---|---|---|
| `L-ECFP` | +0.0199 [+0.0076, +0.0283] | +0.0018 [−0.0048, +0.0055] |
| `L-CHEMBERTA` | +0.0194 [+0.0032, +0.0282] | −0.0001 [−0.0027, +0.0016] |
| `L-DESC` | +0.0158 [+0.0003, +0.0254] | −0.0000 [−0.0033, +0.0013] |
| **`L-RANDOM` (control)** | **+0.0154 [+0.0024, +0.0241]** | −0.0002 [−0.0012, +0.0003] |

**Random ligand features reproduce essentially the entire `SEC` gain over `ADD`.**
The residual `+0.02` is not chemistry and not biology: it is a variance-reduction
artefact of fitting `(b, v)` jointly on the support instead of fitting the
intercept alone, and any random direction serves equally well. Combined with the
derangement result, the double-held-out section carries **no** identifiable
protein-specific or ligand-chemical information.

This also retrospectively explains why §4's genuine ligand-landing signal
(`R2 = +0.199` for ECFP against `+0.025` for random) does not survive into §6:
the loading prediction is informative about the interaction *energy*, but the
part of it that distinguishes one protein from another is smaller than the noise
it carries into `<uhat, v>`.

**The mechanism of the collapse.** XP1 obtained its specificity with the query
ligand's loading taken from a fitted table. XP2 must predict it from structure,
and although that prediction is genuinely informative (§4, `R2 = +0.199`), the
residual noise in `uhat(L)` is multiplicative with `v_j` in `<uhat, v>`. It
degrades the *target-specific* component — the part that distinguishes the
correct protein from a foreign one — faster than it degrades the *shared*
component. What survives contact with unseen scaffolds is the part of the
correction that is not protein-specific.

**Gate conditions 3, 4, 5 and 6 fail.**
`DOUBLE_HELD_OUT_SECTION_IDENTIFIED` is refused.

### 6.3 XP2-F — external replication

`EXTERNAL_REPLICATION_RESULT.json`. Klaeger 2017 kinobeads: independent
laboratory, independent technology (chemical proteomics vs radiometric),
independent compound set, never used for any XP1 or XP2 model-selection
decision. Panel `BLK-KLAEGER-XP2`: **142 drugs x 106 kinases, 2,838 measured
cells, 133 scaffolds**; 3 drugs had unresolvable structures and 126 kinobeads
targets were dropped for having no KLIFS kinase record (they are largely
non-kinase off-targets). Support drugs are scaffold-disjoint from query drugs,
and a reserved third of scaffold components is used only for the calibration
conclusion.

**Conclusion 1 — direction transfer (every source parameter frozen):**

| arm | RMSE | `R2_gamma` vs ADD |
|---|---|---|
| `P0` | 0.8768 | +0.1120 [−0.0510, +0.2510] |
| `LIG` | 0.8864 | +0.0925 [+0.0452, +0.1455] |
| `ADD` | 0.9305 | *(baseline)* |
| `SEC` | 0.9286 | +0.0040 [−0.0026, +0.0096] |
| `FOREIGN` | 0.9302 | +0.0005 [−0.0069, +0.0072] |
| `PERM` | 0.9344 | −0.0085 [−0.0144, −0.0019] |

```
Delta_interaction (ADD - SEC)       = +0.00346 [-0.00234, +0.00863]   CI spans 0
Delta_specific_foreign              = +0.00300 [-0.00305, +0.00904]   CI spans 0
Delta_deploy (LIG - SEC)            = -0.07663 [-0.12184, -0.03695]   negative
```

**Conclusion 2 — basis transfer** (kept strictly separate): the preregistered
global affine calibration fitted on the reserved external components is
`a = 4.205, b = 0.332`. A slope of `0.33` means Metz-scale predictions must be
shrunk threefold to match the kinobeads endpoint — the two endpoints are not on a
common scale, which is the expected consequence of `pKi` versus apparent `pKd`
and is reported rather than absorbed.

**Gate condition 7 fails.** `EXTERNAL_REPLICATION_FAILED`. Note also that on the
external panel the population mean alone (0.8768) beats every adapted arm, and
the Metz-trained ligand main effect actively hurts relative to it — a further
sign that nothing protein-specific crossed the platform boundary.

---

## 7. XP2-E — biological landing

`XP2E_BIOLOGICAL_LANDING.json`. Identical double closure, identical cells,
identical ligand chemistry landing. The arms differ **only** in where the
protein-side coordinate `v_j` comes from. `d = 3`, `k = 5`.

| arm | RMSE | `R2_gamma` vs ADD |
|---|---|---|
| `ADD` | 0.9493 | *(baseline)* |
| **`SUP` support only** | 0.9417 | **+0.0159 [+0.0060, +0.0205]** |
| `PROT-ESM` zero-shot ESM-2 t30 | 0.9486 | +0.0014 [−0.0028, +0.0070] |
| `PROT-POCKET` zero-shot aligned KLIFS 85-mer | 0.9482 | +0.0022 [−0.0035, +0.0091] |
| `PROT-GROUP` zero-shot kinase group | 0.9482 | +0.0023 [−0.0019, +0.0080] |
| **`PROT-RANDOM` control** | 0.9492 | **+0.0002 [−0.0018, +0.0030]** |
| `PROT-ESM+SUP` | 0.9413 | +0.0167 [+0.0057, +0.0216] |
| `PROT-POCKET+SUP` | 0.9411 | +0.0171 [+0.0056, +0.0219] |

Every zero-shot protein representation has a confidence interval containing zero
and sits within `0.002` of the **random** protein embedding. Adding protein
features on top of the support moves `+0.0159` to `+0.0167` / `+0.0171` — a
change of `0.001`, far inside the interval width.

**Biology contributes nothing here; the dataset axes do.** This reproduces XP1-D
on a different panel, with a different construction and under a strictly harder
closure. The correct label for whatever does work is
`PANEL_LOCAL_LOW_RANK_META_LEARNING`, and §6.2 shows even that is generous,
because random ligand features reproduce the support arm's gain. This is
**not** deep biological integration and must not be described as such.

---

## 8. XP2-G — mathematical interface

Full audit in `research/crossed_panel_deployability/THEORY_INTERFACE_AUDIT.md`.
Summary: the candidate seven-tuple is interface-legal **conditional on three
declarations** that the current specification does not yet make.

1. **A declared gauge.** `section_center = <uhat, vhat>` and `support_rank` are
   `GL(d)`-invariant and safe. `query_coverage` and `inverse_conditioning` are
   **not** — they are defined through a projector and a singular-value ratio and
   are meaningful only after the loading map is whitened on a declared corpus and
   the whitening matrix is frozen into `gamma`.
2. **A two-term radius.** A ridge point estimate silently zeroes the unidentified
   component of `v`. An honest enclosure needs an estimation term *plus*
   `||(I − P_S) uhat|| * B` with `B` a declared bound on `||v||`. At `k <= 3`,
   where coverage is 0.34–0.67, a third to two thirds of the query direction is
   not identified at all and the enclosure must widen accordingly.
3. **Correct placement of the discrete coordinates.** `validity_flag` and
   `support_rank` are discontinuous, so `(S-CONT)` requires them to enter through
   the finite context map `kappa`, not through the continuous coordinates the
   multilinear sieve consumes. `model/meta_operator.py::context_index` already
   provides that machinery.

Abstention needs no new object: when `validity_flag = 0` the coefficient map can
place all mass on `p_0`, and the emitted class is exactly the population band for
the declared context — the `e_0` vertex of the existing simplex. CSMO, Band, `K`
and the mesh are untouched.

**No latent coordinate may be given a biological name.** Only fixed named
features and gauge-invariant objects — the scalar centre, the rank, the support
subspace and its projector — are interpretable. §4.1 shows why empirically.

---

## 9. Gate evaluation

Thresholds frozen in `PREREG_XP2.md` §10 before any arm was scored, and not
modified after any test set was opened. Machine-checked in
`DOUBLE_HELD_OUT_RESULT.json`.

| # | condition | required | observed | verdict |
|---|---|---|---|---|
| 1 | XP1 evidence reproduced | audit passes | 18/18 checks, `XP1_EVIDENCE_REPRODUCED` | **PASS** |
| 2 | ligand loading transferable to unseen scaffolds | `R2 > 0`, CI lower `> 0`, beats random and mean | ECFP `+0.199 [+0.133, +0.261]` vs `+0.025` / `+0.024` | **PASS** |
| 3 | `k <= 5` `R2_gamma` above floor | `>= 0.05`, CI lower `> 0.02` | `+0.0199 [+0.0076, +0.0283]` | **FAIL** |
| 4 | double held-out, same thresholds | as (3), test opened once | `+0.0199 [+0.0076, +0.0283]` | **FAIL** |
| 5a | beats zero adaptation | CI lower `> 0` | `+0.00181 [−0.00106, +0.00528]` | **FAIL** |
| 5b | beats foreign support | CI lower `> 0` | `+0.00185 [−0.00477, +0.00552]` | **FAIL** |
| 5c | beats permuted support | CI lower `> 0` | `+0.00135 [−0.00127, +0.00713]` | **FAIL** |
| 6 | `Delta_deploy` over ligand-only | CI lower `> 0` | `−0.03314 [−0.09108, +0.02326]` | **FAIL** |
| 7 | external replication | CI lower `> 0` | `+0.00346 [−0.00234, +0.00863]` | **FAIL** |
| 8 | frozen-theory interface | audit | conditional on three declarations (§8) | **CONDITIONAL** |

**Overall: FAIL.** Seven of eight scored conditions are not met (condition 8 is conditional).
`DEPLOYABLE_SECTION_STATISTIC_IDENTIFIED` is refused.

Nothing was promoted. `model/`, production `scripts/`, `contracts/` and
`theory/` were not modified at any point in XP1 or XP2 — verified by
`git status` on those paths at the end of the stage. The repository regression
suite passed `73 passed` both before and after.

---

## 10. Disposition of implementations

`history.md` records SHA-256 for all 50 XP1/XP2 code and artifact files, plus
dependency versions, seeds, upstream licences and the label-read audit.

**Removed:** `report/crossed_panel_deployability/xp2cd_smoke_protein_d3_k5_L-ECFP.json`
(a duplicate of the `C_k5` sweep entry) and
`report/crossed_panel_identification/xp1b_main_group_r8_k16.json` (superseded by
the corrected two-stage rerun inside `xp1b_sweeps.json`).

**Retained and why.** The XP2 implementations are the *evidence* for a terminal
negative, and XP2-A reproduces XP1 by reading `xp1b_transfer.py` directly, so
deleting either tree would make the negative unverifiable. Both `research/`
trees are currently **untracked in git**, which means deletion is irreversible —
there is no `8b7789e`-style recovery point for them, unlike the earlier F6I
cleanup. The deletion candidates are therefore listed for the maintainer rather
than executed:

| candidate | conclusion preserved in | note |
|---|---|---|
| `research/crossed_panel_identification/xp1c_pdsp.py` | XP1 report §6.4, `xp1c_pdsp.json` | PDSP replication below floor |
| `research/crossed_panel_identification/xp1d_statistic.py` | XP1 report §6.6, `xp1d_statistic.json` | superseded by XP2-E |
| `research/crossed_panel_identification/build_conformation_features.py` | XP1 report §5, §6.4 | conformational-state arm disqualified on availability and circularity |
| `research/crossed_panel_deployability/xp2e_landing.py` | XP2 report §7, `XP2E_BIOLOGICAL_LANDING.json` | terminal negative |
| `research/crossed_panel_deployability/xp2f_external.py` | XP2 report §6.3, `EXTERNAL_REPLICATION_RESULT.json` | terminal negative |

Nothing was promoted into `model/` or production `scripts/`. `git status` on
`model/`, `scripts/`, `contracts/` and `theory/` is clean for the whole
programme, and the regression suite is `73 passed` before and after.

---

## 11. Reproduction

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_deployability/xp2a_reproduction_audit.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_deployability/xp2_panel.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_deployability/xp2b_ligand_landing.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_deployability/xp2cd_sweep.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_deployability/xp2e_landing.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_deployability/xp2f_external.py
```

```bash
OMP_NUM_THREADS=1 D:/anaconda/envs/drug/python.exe research/crossed_panel_deployability/xp2_finalize.py
```

### Registered deviations

1. **Inner CV binning.** `PREREG_XP2` §6 specified leave-one-closure-component-out
   for hyperparameter selection. With 206 training scaffold components this is
   ~1,400 dual ridge solves per arm per fold and is computationally infeasible
   here, so components are binned into 10 grouped folds. Whole components never
   straddle an inner split, so the procedure remains group-safe; it touches
   hyperparameter selection inside the training block only, and no test cell,
   estimand or threshold. Recorded in `xp2_core.ridge_cv`'s docstring.
2. **Identifiability diagnostic corrected mid-stage.** The first implementation
   measured the rank of the *uncentred* support loadings, overstating the
   identifiable dimension by one because the unpenalised intercept absorbs the
   support mean. It was corrected to the centred design and every reported run
   was recomputed. Predictions were unaffected; only the diagnostic changed.
