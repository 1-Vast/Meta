# Stage C: where the error lives, and what MSE ≤ 1.00 pK² would require

**Finding: the ≤ 1.00 pK² target is not reachable at k=0 with the current
inputs, and the missing information is named precisely — zero-shot target-level
affinity calibration for an unseen homology component.**

This is a measurement result, not a training result. Nothing was trained here.
`meta_test` was not read; its audited status is *logical exclusion after parsing
with an open process-isolation incident*.

Authority: `FEASIBILITY.json`, `LEVEL_CEILING.json`, `REPRESENTATION.json`.
Baseline: the leak-free Stage B `T` checkpoint (checkpoints selected on
`meta_train` internal-validation components only).

---

## 1. The error decomposes exactly, and it is mostly level

Every squared error splits with no residual:

```
MSE = (mean(p) − mean(y))²  +  mean(((p − mean(p)) − (y − mean(y)))²)
    =      level²           +            centered_MSE
```

| k | MSE | level² | centered (shape) | level share | **MSE with a perfect level** |
|---|---:|---:|---:|---:|---:|
| 0 | 2.7425 | 1.8664 | 0.8761 | **68%** | **0.8761** |
| 1 | 1.8549 | 0.9788 | 0.8761 | 53% | **0.8761** |
| 2 | 1.3628 | 0.5560 | 0.8068 | 41% | **0.8068** |
| 3 | 1.2481 | 0.4500 | 0.7981 | 36% | **0.7981** |
| 5 | 1.0096 | 0.2754 | 0.7342 | 27% | **0.7342** |

Two things follow immediately.

**A perfect target-level predictor would put every k below 1.00.** The target is
arithmetically reachable, and reachable through level calibration alone — no
ordering improvement required.

**Support labels are already the level mechanism.** level² falls 1.87 → 0.28 as
k goes 0 → 5, while the shape term barely moves (0.876 → 0.734). The few-shot
gain this project has been measuring for many stages is target-level
calibration, arriving through the labels.

### The model orders no better than a constant

k=0 centered MSE is **0.8761** against a within-target label variance of
**0.8525**. The ratio is **1.0277** — the model's within-target ordering is
*worse* than predicting each target's own mean. That is consistent with Stage R
(inert operator, query spread 0.0027 pK), Stage P (protein response with +0.022
truth alignment) and Stage L2 (weak directional signal), measured here at the
endpoint.

## 2. Target level is not predictable from any representation tested

Predicting `mean_ligands(pK)` for an unseen protein, `meta_val`, 41 targets.
Weight decay chosen on `meta_train` component folds; `meta_val` read once.

| method | level MSE | RMSE | vs calibrated constant |
|---|---:|---:|---:|
| ESM linear probe | 6.5368 | 2.557 | **4.85× worse** |
| ESM nearest neighbour | 5.1292 | 2.265 | 3.81× worse |
| sequence-length probe | 2.5283 | 1.590 | 1.88× worse |
| `meta_train` grand mean | 2.1703 | 1.473 | 1.61× worse |
| incumbent model | 1.7078 | 1.307 | 1.27× worse |
| **ESM MLP probe** (best) | **1.6357** | **1.279** | **1.21× worse** |
| *calibrated constant — REFERENCE, reads meta_val labels* | *1.3471* | *1.161* | *1.00* |
| oracle | 0 | 0 | 0 |

**Not one method beats a constant.** The weight-decay sweep selected the
*largest* value (1.0) — maximal regularization, i.e. the fit is being driven
toward a constant because the features carry no cross-component level signal.
Five independent approaches agree, including the full 1.8M-parameter model
trained end-to-end.

There is also a genuine covariate shift: the `meta_train` grand mean scores
2.1703 against the calibrated constant's 1.3471, so the two splits' target-level
distributions differ by **≈ 0.91 pK**.

### The arithmetic that closes k=0

`k=0 MSE = level² + centered`, and `centered ≥ 0`. Therefore

> **k=0 MSE ≥ 1.6357 using the best legitimate level predictor, even with
> perfect within-target ordering.**

Reaching 1.00 needs level MSE ≤ `1.00 − 0.8761 = 0.1239`, a **13.2× reduction**
from the best measured predictor, i.e. explaining ~91% of the between-target
variance (1.3471) that nothing currently explains at all.

### Why this is plausibly irreducible rather than a modelling failure

A BindingDB target's mean pKi depends heavily on *which ligands were tested
against it*: a mature medicinal-chemistry programme contributes optimized,
high-affinity compounds, an early screening campaign contributes weak ones. Part
of "target level" is therefore a property of the assay history and library
composition, not of the protein sequence — and no sequence-derived
representation can recover it. This is a hypothesis consistent with the
measurements, not something these measurements prove.

## 3. The ligand-varying signal is not absent — and one representation carries it

Stage B reported within-target cosine 0.997 and inferred collapse. Cosine alone
was insufficient evidence, and the fuller measurement changes the picture.

| representation | width | protein-constant var | ligand-varying var | ligand share | eff. rank | rel. separation |
|---|---:|---:|---:|---:|---:|---:|
| `embed` | 96 | 0.5816 | 0.1352 | 0.189 | 4.11 | 0.0327 |
| `readout_hidden` | 96 | 0.0371 | 0.0076 | 0.171 | 3.85 | 0.0317 |
| **`occupancy`** | **24** | 0.1435 | 0.0694 | **0.326** | 2.06 | **0.0533** |
| `section` | 48 | 0.2003 | 0.0441 | 0.180 | 4.11 | 0.0267 |

17–33% of each representation's variance is ligand-varying, spread over roughly
2–4 effective directions. So there *is* a ligand-varying subspace.

The direct test — a frozen linear probe fitted on `meta_train` panels, weight
decay chosen on component folds, predicting within-target centered affinity:

| representation | train-fold r | **meta_val r** |
|---|---:|---|
| `embed` | +0.2383 | +0.0074 [−0.1550, +0.1718] |
| `readout_hidden` | +0.2301 | +0.0256 [−0.1375, +0.1936] |
| **`occupancy`** | +0.2029 | **+0.2182 [+0.0751, +0.3670]** |
| `section` | +0.2174 | +0.0603 [−0.1026, +0.2250] |

**`occupancy` is the only representation whose within-target ordering signal
survives to unseen components, and its interval excludes zero.** The three wide
representations overfit — ~0.23 in-fold collapsing to ~0.01–0.06 out of
component — while the 24-dimensional contact-type vector holds its value.

### The model has this signal and does not use it

`occupancy` reaches the endpoint through exactly one path: `contact_weight`, a
single `Linear(24 → 1)`. The probe above has the *same* capacity and extracts
r = +0.218, while the model's own endpoint orders at ratio 1.0277 — no better
than a constant.

The likely reason is that `contact_weight` is trained jointly against a total
MSE that is 68% level error at k=0, so the optimizer spends those 24 parameters
on level rather than on ordering. That is a mechanistic hypothesis, and it is
exactly the separation the Stage 2 plan proposes: a protein-conditioned level
representation kept apart from a ligand-dependent shape representation.

**But it cannot close the gap.** Fully exploiting r = 0.218 would reduce the
centered term by `1 − r² ≈ 4.8%`, from 0.8761 to about 0.834 — roughly 0.04 pK².
Combined with the best level predictor, k=0 would still be ≈ 2.47.

## 4. The boundary, stated precisely

| claim | status |
|---|---|
| k=0 MSE ≤ 1.00 with current inputs | **not reachable** — bounded below by 1.6357 |
| k=1…5 MSE ≤ 1.00 in principle | reachable *only* if ordering becomes near-perfect; k=5 is already 1.0096 |
| what is missing | **zero-shot target-level affinity calibration for unseen homology components** |
| what is present but unused | a resolved within-target ordering signal in `occupancy` (r = +0.218), worth ≈ 0.04 pK² |

Everything the project has attributed to few-shot mechanism — Stage A's inner
loop, Stage B's hybrid, the incumbent's Tanimoto transport — is target-level
calibration arriving through the support labels. The ordering half has never
moved, and the measurement above says the endpoint currently extracts none of
the ordering signal that its own `occupancy` representation contains.

## 5. What would change the verdict

Not a training schedule, and not an adapter. Only new information about the
recipient protein's affinity level:

1. **external protein representations reported as external data** — MSA depth
   and conservation (the preregistered but never-run M0 lane), or structure-derived
   pocket descriptors;
2. **assay/library metadata**, if the level really is partly a property of the
   testing history rather than the protein — this is directly testable by
   regressing target level on assay covariates within `meta_train`;
3. **accepting a different target**: `centered` MSE, CI and Spearman are the
   metrics the current inputs can actually move, and the `occupancy` lever is
   the one measured, resolved, unexploited opportunity.
