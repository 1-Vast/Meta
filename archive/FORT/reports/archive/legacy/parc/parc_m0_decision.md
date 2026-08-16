# PARC M0 decision — `PARC_M0_COORDINATE_NOT_LOAD_BEARING_STOP`

**Run:** 2026-07-27, seed 1729, CUDA (RTX 4060), 223.2 s.
**Preregistration:** `reports/active/parc_m0_preregistration.md` (frozen before any statistic below).
**Runner:** `research/parc_m0.py`. **Result:** `reports/active/parc_m0.json`. **Tests:** `tests/test_parc_m0.py` (8 passed).
**Firewall:** `panel_development_labels_read=false`, `panel_confirmation_labels_read=false`,
`davis_panel_touched=false`, `chembl_confirmation_labels_read=false`, `sealed_test_consumed=false`.
Train cells only.

---

## 1. What was asked

Six structurally independent parameterisations of the interaction have failed on the **same** protein
coordinate, and HQ-GBMA Stage D found the ESM-conditioned Grassmann map worse than a protein-free shared
basis. M0 re-ran Stage D with **exactly one changed input** — the protein coordinate — to separate *"the
protein is uninformative"* from *"the coordinate was wrong"*.

PARC = one-hot of the 85 aligned KLIFS pocket residues (1700 → 32 by PCA whitening), against pooled
ESM-2 (1280 → 32) and the protein-free shared global subspace. Rank 6, ridge 1.0, 400 map steps, 5 frozen
component folds — all carried over from Stage D unchanged, with **no rank re-selection for the new arm**.

Substrate: Metz dense pKi panel, train cells only. 111 of 112 targets eligible (CHEMBL5026 / O00418
eEF2K has no 85-residue KLIFS pocket and was dropped from **every** arm, so all contrasts are matched);
12,560 cells, 100 components, 77 with evaluable containment. All 111 pockets are distinct
(mean pairwise pocket identity 0.373). Projection converged; orthonormality error 4.8e-7; positive-signal
fraction 0.775; derangement fixed points 0.

## 2. G0 — the estimator is sensitive, and specific to the coordinate

Before any real contrast was read, synthetic coefficients were generated whose subspace is a **known
function of the PARC coordinate**, with noise matched to the real per-target estimation noise.

| synthetic contrast | mean | 95% CI |
| --- | --- | --- |
| `parc_pocket − shared_global` | **+0.0431** | [+0.0185, +0.0696] |
| `esm_pooled − shared_global` (specificity) | +0.0206 | [−0.0037, +0.0468] |

**G0 PASS.** The estimator recovers coordinate-dependence when it exists (≥ the 0.02 threshold, LCB > 0),
and it does *not* attribute that dependence to the wrong coordinate. **The null below is therefore a real
null, not a broken estimator.** No prior gate in this program carried this control.

## 3. Result — all five substantive gates fail

Preregistered gates use the component-bootstrap **mean**, carried over from Stage D unchanged.

| gate | contrast | mean | LCB95 | requirement | result |
| --- | --- | ---: | ---: | --- | --- |
| G1 | `parc − esm` | −0.6608 | −2.1929 | LCB > 0 | **FAIL** |
| G2 | `parc − shared_global` | −0.0956 | −0.2887 | mean ≥ 0.02 and LCB > 0 | **FAIL** |
| G3 | `parc − random_positions` | +0.2712 | −0.0053 | LCB > 0 | **FAIL** (by 0.0053) |
| G4 | `parc − wrong_target` | +0.1370 | −0.1614 | LCB > 0 | **FAIL** |
| G5 | `parc − random_features` | −0.2120 | −0.9460 | LCB > 0 | **FAIL** |

**Verdict `PARC_M0_COORDINATE_NOT_LOAD_BEARING_STOP`** under the frozen rule. **M1 is not authorised.**
No threshold, rank, dimension, step count, seed or fold was changed after seeing any result, and no arm
was added after the fact.

## 4. A metric defect was found, and it is inherited from Stage D — not introduced here

`containment_fraction = inside / signal` with `signal = beta_hat^T beta_hat − tr(V_t)` is an **unbounded
ratio**: the denominator is bounded below only by zero. On this substrate the minimum positive signal is
**0.00309** against a median of **0.1473** — a factor of 48. A single component reached containment
**56.09** in the ESM arm and moved that arm's mean from 0.151 to 0.878 on its own.

Robust statistics are therefore reported **as non-gating diagnostics** (no gate reads them):

| arm | mean | median | 5% trimmed | frac. positive | max abs |
| --- | ---: | ---: | ---: | ---: | ---: |
| `shared_global` | +0.3125 | **+0.2112** | **+0.2871** | 0.87 | 4.31 |
| `parc_pocket` | +0.2169 | +0.1350 | +0.1741 | 0.82 | 3.15 |
| `parc_pocket_composition` | +0.6064 | +0.1257 | +0.1808 | 0.86 | 30.28 |
| `parc_wrong_target` | +0.0799 | +0.1256 | +0.1215 | 0.77 | 8.13 |
| `esm_pooled` | +0.8777 | +0.0937 | +0.1464 | 0.83 | **56.09** |
| `esm_wrong_target` | −0.0691 | +0.0958 | +0.1228 | 0.77 | 16.48 |
| `random_features` | +0.4289 | +0.0908 | +0.1332 | 0.86 | 24.79 |
| `parc_random_positions` | −0.0543 | +0.1086 | +0.1064 | 0.81 | 14.28 |

| contrast | mean | LCB95 | median | 5% trimmed | frac. positive |
| --- | ---: | ---: | ---: | ---: | ---: |
| `parc − esm` | −0.6608 | −2.1929 | +0.0253 | +0.0148 | 0.55 |
| `parc − global` | −0.0956 | −0.2887 | −0.0243 | −0.0738 | 0.39 |
| `parc − random_positions` | +0.2712 | −0.0053 | +0.0317 | +0.1113 | 0.57 |
| `parc − wrong_target` | +0.1370 | −0.1614 | +0.0269 | +0.0886 | 0.58 |
| `parc − random_features` | −0.2120 | −0.9460 | +0.0246 | +0.0487 | 0.58 |
| `esm − global` | +0.5652 | −0.3202 | −0.1165 | −0.1551 | 0.31 |

**The verdict is unchanged under either reading.** Robustly the ordering is
`shared_global (0.211) > parc_pocket (0.135) > esm_pooled (0.094) ~ random_features (0.091)`: PARC is a
slightly better coordinate than pooled ESM (median contrast +0.025, positive in 55% of components) and
**neither beats the protein-free shared basis** (`parc − global` positive in only 39% of components,
`esm − global` in only 31%). G1 and G2 fail on the mean and on every robust statistic alike.

### 4.1 Correction to the Stage D record

Stage D reported `true − global = −0.232 [LCB −0.083]`. The identical estimator with **one of 112 targets
removed** gives a mean of **+0.565**, i.e. the *sign of the mean flips* under a one-target perturbation,
entirely because of the 56× outlier. On robust statistics Stage D's qualitative finding reproduces
cleanly and in the same direction (`esm − global` median −0.117, positive in 31% of components), so
`HQGBMA_STAGE_D_FAIL_STOP` stands as a verdict. **Its point-estimate magnitudes should not be cited as
effect sizes.** This applies to every containment number in `reports/active/hqgbma_stage_d.json`.

**Registered blueprint correction (not applied retroactively to any verdict):** any future use of the
containment estimand must aggregate as a **ratio of sums** (`sum_t inside_t / sum_t signal_t`, a
ratio-of-means estimator) rather than as a mean of per-target ratios, and must report the
signal-denominator distribution. M0 did not re-decide anything on that estimand, because the robust
statistics above already give the same ordering.

## 5. What this shows, and what it does not

**Shows.** On the one substrate in this program where the interaction is demonstrably real, low-rank and
identified (`PA5 p = 0.000488`, `PA4 LCB95 +0.2431`), replacing pooled ESM-2 with a positional,
structure-aligned pocket coordinate moves the protein arm in the predicted direction by a **negligible
and unresolvable amount**, and does not close the gap to a protein-free shared global basis. The
subspace-misspecification diagnosis (`task.md` §2.4) is **not repaired by the pocket restriction** here.
The coordinate was not the whole problem.

**Does not show.**
1. **Not** that pocket information is useless. G3 (`parc − random_positions`, the decisive control that
   holds protein identity and sequence source fixed and destroys only the pocket) is directionally
   positive — median +0.032, trimmed +0.111, positive in 57% of components, mean +0.271 with LCB
   −0.0053, missing its bar by 0.0053. The pocket restriction is neither shown load-bearing nor shown
   inert. It is unresolved at 77 components.
2. **Not** a test of the structure-token half of PARC. The Foldseek/3Di arm (A3) was **not run**: no 3Di
   tooling exists in the `drug` env and `dataset/structure/alphafold` was removed in the 2026-07-25
   cleanup. M0 tested the positional/pocket half only. The KLIFS alignment is structure-derived, but
   discrete local-structure tokens were never evaluated.
3. **Not** a multi-family result. KLIFS is kinase-only; every eligible target is a kinase. The
   within-kinome limitation that closed TR-0 and PFSC-0 applies unchanged.
4. **Not** a predictive result. The estimand is subspace containment of per-target empirical
   coefficients, train-only. No affinity prediction was made, no development/confirmation label was read,
   and nothing here may be cited as evidence about dual-cold ranking accuracy.
5. **Not** a refutation of the declared expected outcome in the optimistic direction — the declared
   expectation was "G1 passes, G2 fails". **G1 also failed.** The result is weaker than predicted, and
   the prediction is recorded as wrong.

## 6. Consequences

* `PARC M1` (the `gamma`-interpolation estimator) is **not authorised**. Part 8 stops at M0.
* No predictive claim, no F1–F4 unlock, no Mamba comparison, under any reading.
* Reopening the PARC coordinate requires **new information**, not a re-run: specifically (a) the 3Di /
  structure-token arm on an audited dependency, **and** (b) a multi-family substrate on which a pocket
  coordinate can be tested outside the kinome — the same measurement-design condition the program has
  already reduced to. Increasing width, rank, steps, dimension or seeds is not admissible.
* The model-side innovation budget is **not** consumed by this stop: P1 was tested and failed its
  identifying gate; P2 (`gamma`-interpolation) was never built.
