# MetaSieve — Terminal Identifiability Report

Date: 2026-08-08
Stages: XP1 (crossed-panel identification), XP2 (deployability), XP3 (data and
observability census), XP4 (many-panel interaction), XP5 (fixed typed basis).
Environment: `drug`, Python 3.11.15, numpy 1.26.4, scipy 1.17.1, pandas 2.3.3,
rdkit 2023.09.6, torch 2.6.0+cu124, transformers 4.46.3.
Label-read counters: DAVIS `0`, recipient `0`, ChEMBL37 affinity `0`,
PKIS2 `0`, Anastassiadis `0`.
GPU use: frozen ESM-2 t30 and ChemBERTa **inference only**, never training. No
GPU training was performed because no stage ever reached the Stage-4 trigger
condition — no basis passed a structural or affinity Gate that a lightweight
readout then underfit.

---

## 0. Terminal outcome

```text
PUBLIC_DATA_INSUFFICIENT_FOR_IDENTIFICATION      (designs with adequate independent units)
DEPLOYMENT_INPUTS_INSUFFICIENT                   (the panel that does contain signal)
```

This is terminal outcome **B**: a rigorous stopping boundary. It is not "one
model failed". It is the conjunction of a measured data limitation and a measured
representation limitation, established on four independent designs with
preregistered Gates that were never moved.

**What is NOT claimed:** that protein-specific interaction does not exist (it
demonstrably does), or that no representation could ever capture it. The bounded
claim is stated in §7.

---

## 1. The ten questions

| # | question | answer | evidence |
|---|---|---|---|
| 1 | Does a real protein-specific affinity interaction exist above noise? | **Yes, in single-laboratory profiling panels.** 59.6% of Metz panel variance; 38% of the residual reproducible; protein-side geometry replicates at `r=0.885` across disjoint compound halves and `r=0.565` across an independent platform (`p<5e-4`) | XP1-A |
| 2 | Which biological inputs are needed to observe it? | **Unknown, and not any tested.** Every deployment-observable representation is indistinguishable from a capacity-matched random control under proper closure | XP1-D, XP2-D/E, XP4, XP5 |
| 3 | Can it be represented by a fixed deployable basis? | **No, at the pose-free rung.** Ten named physicochemical complementarity channels give `R2_gamma = -0.0015` and derangement specificity of exactly `+0.00000 [-0.00069,+0.00066]` | XP5 |
| 4 | Can `k <= 5` support identify its task-specific section? | **No.** Identifiable dimension is exactly `min(k-1,d)` — measured `0,1,2,3,3` at `k=1..5`; `R2_gamma` peaks at `+0.025` against a frozen `0.05` floor | XP2-C |
| 5 | Does correct support beat ligand-only, foreign support and random features? | **No** under double closure. Derangement `+0.00185 [-0.00477,+0.00552]`; random ligand features reproduce the entire gain (`+0.0154` vs `+0.0199`) | XP2-D |
| 6 | Does it transfer across both protein groups and ligand scaffolds? | **No.** Specificity present under protein-only closure vanishes entirely when scaffold closure is added | XP2-C vs XP2-D |
| 7 | Does it replicate externally? | **No.** Klaeger direction transfer `+0.00346 [-0.00234,+0.00863]` | XP2-F |
| 8 | Deep biological meta-learning, or a panel-local statistical learner? | **Panel-local at best, and even that is generous** — random ligand features reproduce the gain | XP2-D §6.2 |
| 9 | Does the statistic satisfy the frozen theory? | Interface-legal only conditionally, and moot: no statistic passed an empirical Gate | XP2-G |
| 10 | Admitted into `z` and the probability law? | **No.** `K(B(z)F(z))` was never scored, deliberately — there was never an admitted statistic to score it with | all stages |

---

## 2. The governing trade-off in public data

This is the central finding of the census, and it explains every negative above.

| design | independent protein-side units | interaction vs noise | verdict |
|---|---|---|---|
| single-lab kinase profiling (Metz) | **8** — capped by Manning kinase taxonomy | real: 38% of residual reproducible | signal, too few units, fails deployability |
| single-lab kinobeads (Klaeger) | ~8 | 93.6% of cells at the floor — a hit matrix | not a continuous crossed panel |
| literature aggregation (PDSP) | many | per-report `sigma = 0.714`; radioligand is nearly a function of the target | assay confounded with the estimand |
| literature panels (BindingDB curated) | **70 protein clusters / 85 panels** | per-report `sigma = 0.777`; `gamma` sd `0.406` is **below** the `0.650` that noise alone transmits | **below the noise floor** |

**Public data offers either low noise with few independent protein units, or many
independent units with noise exceeding the interaction. No accessible source
offers both.** That is the stopping boundary, and it is a property of how the
data were generated, not of any model.

The BindingDB result is quantitative and independent of any model: the most
favourable non-parametric use of each panel's own data — similarity-weighted
neighbour transfer within the same panel and the same target — scores
`R2_gamma = -0.539` on 2,068 cells. There is nothing there to fit.

---

## 3. What each stage falsified

| stage | hypothesis tested | verdict |
|---|---|---|
| XP1 | a support-identified latent section transfers to unseen protein groups | **partly true** at `k=16` with ligands reused; zero-shot protein landing fails |
| XP2-A | XP1's evidence is reproducible | **`XP1_EVIDENCE_REPRODUCED`**, 18/18 checks |
| XP2-B | the ligand loading is a lookup table | **falsified** — `LIGAND_LOADING_RECOVERABILITY_OBSERVED`, ECFP `R2 = +0.199 [+0.133,+0.261]` vs random `+0.025` |
| XP2-C | the section survives `k <= 5` | **failed** — below the frozen floor at every `k` |
| XP2-D | it survives double closure | **failed** — specificity exactly zero; random ligand features reproduce the gain |
| XP2-E | protein representations contribute | **failed** — all within `0.002` of a random protein embedding |
| XP2-F | it replicates externally | **failed** |
| XP4 | a deployment-observable bilinear basis works given many independent panels | **failed, because the data is below its own noise floor** |
| XP5 | a fixed named typed basis, computed rather than predicted, carries the interaction | **failed** — specificity exactly `+0.00000` |

XP5 deserves emphasis because it removed XP2's diagnosed failure mechanism. XP2
failed partly because `u(L)` had to be *predicted*, injecting noise. XP5's channels
are *computed* analytically from both partners — no prediction, no gauge freedom,
no lookup. It still returns nothing. The failure is therefore not the estimator.

---

## 4. What was deliberately not attempted, and why

**Rung `B2` (experimental or predicted 3D structure with docked poses)** is the
only rung of the ladder left untested. It was not attempted for three stated
reasons, not for convenience:

1. **Compute.** The estimand needs pose-derived features for the crossed cells:
   928 ligands x 147 kinases = 136,416 dockings. At a realistic 10 s per pose
   that is ~380 GPU-hours, far outside a lightweight staged programme.
2. **Subsampling defeats the purpose.** Any tractable subsample reduces
   protein-group units below the 8 that already proved marginal.
3. **The project has already tested the geometry rung.** `P1B` validated
   correct-partner contact/distance geometry (correct AUPRC `0.43885`,
   wrong-protein `0.05149`) and `P1C` then found that this validated geometry
   yields no affinity direction. XP1 independently found that the
   structure-aligned KLIFS pocket and KLIFS conformational-state availability are
   both indistinguishable from a random null under group closure. Pose error
   would inject exactly the prediction noise that XP2 showed destroys the
   target-specific component preferentially.

`B2` is therefore **recorded as untested with a stated feasibility barrier**, not
as refuted. §7 states the missing experiment.

---

## 5. Theory position

`THEORY_INTERFACE_AUDIT.md` stands: the candidate statistic is interface-legal
only conditional on a declared gauge, a two-term outer radius whose second term
bounds the unidentified component, and placement of the discrete coordinates in
the finite context map `kappa`. Abstention is representable as the existing
`p = e_0` simplex vertex and needs no new operator.

Two theory-side results are worth retaining independently of the negatives:

- **The identifiability ledger.** With an unpenalised support intercept, the
  identifiable section dimension is exactly `rank(U_S - mean U_S) = min(k-1, d)`.
  A frozen `k <= 5` therefore caps any support-identified section at four
  dimensions and gives **zero** at `k = 1`, where ridge returns `v = 0` exactly.
  Any future design must budget for this.
- **The gauge result.** Coordinate-wise loading `R2` is ~0 or negative in the same
  runs where the gauge-invariant reconstruction `R2` is `+0.199`. Latent
  coordinates are not stable objects and must never be given biological names.

`model/`, production `scripts/`, `contracts/`, `theory/`, CSMO, Band, `K` and the
fixed mesh were not modified at any point across all five stages.

---

## 6. Honest limitations of this conclusion

1. **Metz is now a development panel** and was used across XP1, XP2 and XP5. XP5
   is a single registered configuration with no sweep, but a *positive* result
   there would have required external replication; the negative is less
   vulnerable to reuse, since reuse inflates rather than deflates apparent signal.
2. **Eight protein-group units** underpin every Metz interval. The intervals are
   correspondingly wide, and a small true effect could hide inside them — but the
   Gate floor was set at `0.05` precisely to make "too small to matter" a
   failure rather than an ambiguity.
3. **BindingDB's noise estimate (`sigma = 0.777`) rests on 298 replicated cells**
   and mixes within-paper duplicate curation with genuine replication. It is an
   upper bound on within-assay noise. The chemistry-ceiling result
   (`R2 = -0.539`) does not depend on it and reaches the same conclusion.
4. **`B2` is untested**, as stated in §4.
5. The endpoint throughout is `Ki`/`Kd`-type log-affinity. Nothing here speaks to
   kinetic, thermodynamic or cellular endpoints.

---

## 7. The exact missing experimental design

The stopping boundary identifies precisely what public data does not provide, and
therefore what would be needed to answer the question:

> **A single-laboratory, single-protocol crossed panel spanning at least 30
> homology-disjoint protein clusters from more than one protein family, with at
> least 200 scaffold-diverse ligands measured against every target, quantitative
> `Ki` or `Kd`, replicate measurements sufficient to put the per-cell `sigma`
> below ~0.3 log units, and uncensored or explicitly interval-censored values.**

Order of magnitude: ~30 x 200 = 6,000 cells with replicates. That is one
well-designed profiling campaign. It does not exist in the open literature in the
form required, which is exactly why the census returns insufficient.

Failing that, the next tractable *research* step is rung `B2` on a deliberately
reduced design — accepting fewer protein units in exchange for pose-derived
features — registered in advance as a feasibility probe rather than a Gate.

---

## 8. Artifacts

| stage | preregistration | artifacts |
|---|---|---|
| XP1 | `research/crossed_panel_identification/PREREG_XP1.md` | `report/crossed_panel_identification/` |
| XP2 | `research/crossed_panel_deployability/PREREG_XP2.md` | `report/crossed_panel_deployability/` |
| XP3 | census, label-blind | `report/observability_census/` |
| XP4 | `research/multipanel_interaction/PREREG_XP4.md` | `report/multipanel_interaction/` |
| XP5 | `research/typed_basis/PREREG_XP5.md` | `report/typed_basis/` |

Recovery commit for the complete reproducible tree: `3281780` on branch
`research/xp1-xp2-crossed-panel`. All release checksums, dependency versions,
seeds, licences and label-read counts are recorded in `history.md`.
