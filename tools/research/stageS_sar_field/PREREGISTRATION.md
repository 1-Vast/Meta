# Stage S preregistration — cross-target protein-conditioned SAR field

Frozen **before any arm was trained**. Phase 0 (`PHASE0_AUDIT.json`) and Phase 1
(`tests/test_structural.py`, 24 tests) were complete when this file was written;
no evaluation metric from any trained arm had been computed. Nothing in this
file may be changed after the first arm's internal-validation metric is read.
A post-hoc threshold change voids the stage.

## 1. The hypothesis under test

A ligand encoder learns a **target-independent** shared coordinate `phi(L)`. The
directed chemical change is `u_ab = phi(L_b) - phi(L_a)`; no protein information
enters `phi` or `u`. A protein encoder produces response coefficients
`alpha(P)`. A nonlinear SAR potential `V(P,L)` is FiLM-modulated by `alpha`, and
the prediction is a potential difference:

    dy_hat(P,a,b) = V(P,L_b) - V(P,L_a)

**Question.** Does the correct protein carry information about the *signed*
within-target affinity difference between two ligands, over and above what a
protein-blind field extracts from chemistry alone?

The repository's prior is that it does not. Stage P drove a centered
protein-conditioning objective at 4x weight and got a reproducible but
truth-unaligned protein response (alignment +0.022 against a +0.10 threshold);
the A2 exact operator's wrong-protein control failed inverted; Stage L2 found a
directional SAR signal that is **protein-independent** (`embed` +0.2119 ± 0.0112)
and strongest exactly where fixed Morgan/Tanimoto transport is already strong.
This stage is not an attempt to rescue any of those. It is a different function
class — a conservative potential field trained directly on signed pair
differences, with the between-target level term algebraically removed — tested
against a frozen gate that the prior predicts it will fail.

## 2. What the construction guarantees, and what is excluded

Because the prediction is the difference of a scalar potential, it is a
curl-free field and the three identities hold for **every** parameter setting,
before and after training:

| identity | status |
|---|---|
| `dy_hat(P,a,b) = -dy_hat(P,b,a)` | bitwise exact in IEEE-754 (negation is exact) |
| `dy_hat(P,a,a) = 0` | bitwise exact (self-subtraction is exact) |
| `dy_hat_ab + dy_hat_bc + dy_hat_ca = 0` | exact in exact arithmetic; measured `< 1e-12` in float64 |

Excluded by construction, each with a test:

* **no explicit `e^T H e` quadratic term.** It is even under `e -> -e`, so it
  cannot appear in a signed prediction. The potential is a function of one
  ligand coordinate and never receives a difference at all.
* **no double protein conditioning.** `phi` has no protein argument; `alpha` has
  no ligand argument; they meet exactly once, inside the potential.
* **no 3D geometry claim.** `phi` is a 2D graph + Morgan-fingerprint coordinate
  and `alpha` is a frozen ESM-2 150M sequence readout. Nothing here is a bound
  pose, a docking result, or a protein-ligand atomic contact, and no artifact of
  this stage may describe it as one.
* **no closed-form adaptation.** Ordinary AdamW forward/backward only: no ridge,
  no pseudoinverse, no analytic solver, no test-time query-label gradient, no
  second dataset, no support set of any kind.

## 3. Data and governance

* Corpus: governed BindingDB-Ki `main_v0`, double-cold split `v1`, mounted
  through the **physically isolated split view**
  (`bindingdb_ki_double_cold_v1_views`). The `meta_test` label artifact is not
  present on the mounted surface, so no sealed label is decompressed or parsed
  by any process in this stage. `seal_record()` is copied verbatim into every
  artifact.
* The development-validation split (`meta_val`) is **not opened anywhere in this
  stage** — not for training, not for selection, not for reporting, not for a
  census. A parsed-AST test fails the suite if its name appears as a string
  constant in any stage module.
* Partition: the frozen `scripts/internal_validation.py` component partition of
  `meta_train` (`PARTITION_SEED = 20260818`, 12% held out): **227 fit
  components / 31 internal-validation components**. Training reads fit only;
  reporting reads internal-validation only.
* **Checkpoint selection: none.** Every arm trains for exactly the same number of
  steps and the final parameters are evaluated. Selecting on the
  internal-validation components would leak into the reported population;
  selecting on `meta_val` is the defect this repository measured at ~0.62 pK^2.

### Disclosure that travels with every number in this stage

The internal-validation partition is **protein-component-hard but not ligand-cold**.
Phase 0 measured that **183 of 638 (28.7%)** internal-validation ligands are
identity-shared with the fit components and **51.4%** share a Murcko scaffold.
That is a property of partitioning `meta_train` by protein homology, which is
what the task specifies; the double-cold `meta_val`/`meta_test` populations do
not have it. Consequences: (a) absolute numbers here are **not** comparable with
any k=0/k=5 figure on the double-cold protocol; (b) the chemical-novelty strata
below are mandatory, not optional, and the most-novel tercile is the honest
read of transfer.

## 4. Phase 0 result that shapes the design

| population | components | targets | pairs | same-panel | cross-panel | cliffs |
|---|---:|---:|---:|---:|---:|---:|
| fit | 227 | 290 | 138,695 | 42,195 | 96,500 | 3,372 |
| internal validation | 31 | 42 | 13,279 | 6,775 | 6,504 | 563 |

* Pairs per target are extremely skewed (median 36, max 26,335), so **global
  pair sampling would make the loss a report about three targets**. All banks
  are target-balanced and stratum-balanced.
* Same-panel `|dpK|` is 0.943 mean against cross-panel 1.286: cross-panel
  differences carry inter-assay offset on top of chemistry. They are therefore
  **excluded from the training loss** and reported as their own stratum.

## 5. Arms (matched budget, one code path, one seed)

| arm | protein input | counterfactual loss | labels |
|---|---|---|---|
| **A** `A_ligand_only` | none (learned constant response) | no | true |
| **B** `B_protein` | correct protein | no | true |
| **C** `C_protein_cf` | correct protein | yes | true |
| **D** `D_protein_shuffled` | stable cross-component permutation | no | true |
| **E** `E_label_shuffled` | correct protein | no | within-target permuted |

Frozen hyperparameters, identical in every arm: seed 20260819; **4,000 steps**;
16 targets x 8 pairs = 128 pairs per step; AdamW lr 3e-4 with cosine decay,
weight decay 1e-4, gradient clip 5.0; training bank = same-panel fit pairs,
target-balanced, <= 512 per target (24,125 pairs over 285 targets);
loss = Huber(delta=1.0) on signed dpK + 0.2 x softplus sign term;
counterfactual weight 0.5 with a 0.25 pK^2 hinge margin (arm C only).
Model: `phi` 64-d, `alpha` 64-d, potential width 128 depth 3; 523,524 parameters
protein-conditioned, less the protein encoder in arm A (reported per arm).

The step budget was fixed from a **timing-only** smoke (50 steps, wall clock and
parameter count only; no evaluation metric was computed). B and C are the two
declared candidate arms; reporting two doubles the multiplicity and that is
stated rather than hidden — both are reported in full whatever they show, and
each must clear every gate on its own.

### Hard wrong proteins

Selected by a frozen rule (`features.hard_wrong_protein_map`): a different
CD-HIT40 component; drawn from `meta_train` and chosen using `meta_train`-only
protein features (frozen pooled PLM vectors and the governed component map);
the **most similar admissible protein by PLM cosine**, so it is the hardest
available negative rather than a random one; and with **zero DOI overlap**
between the recipient's and the donor's panel documents, so the swap cannot be
solved by recognising a shared testing programme. Only the protein input is
replaced — the recipient's ligands and its `delta_y` are bitwise unchanged.

## 6. Evaluation

Population: the withheld internal-validation components. Primary bank:
**same-panel** pairs, target-capped at 512 (4,270 pairs / 42 targets /
29 components). Secondary bank: cross-panel pairs (2,028). Diagnostic bank:
same-panel fit pairs the training draw did not take — a training-health monitor
that never selects anything.

Reported for every arm: equal-component target-mean MSE; pooled MSE; Pearson;
Spearman; CI (Kendall tau-b rescaled); sign accuracy; activity-cliff sign
accuracy; same-panel and cross-panel strata; chemical-novelty terciles cut at
the Phase 0 quantiles **0.3014 / 0.5607** of a pair's mean
max-Tanimoto-to-fit-ligands; correct vs ligand-only; correct vs matched-wrong
protein; correct vs shuffled-protein training; correct vs shuffled-label
training; component-paired bootstrap intervals (2,000 draws, resampling the 29
components, each drawn copy counted separately).

Primary Pearson is **pooled over the internal same-panel bank**; the
per-target-then-equal-component version is reported as a robustness check.

## 7. Decision thresholds — FROZEN

Let `X` be a candidate arm (B or C). All contrasts are component-paired
bootstraps on the internal same-panel bank unless stated.

**G1 — protein beats ligand-only.**
`Pearson(X, correct) - Pearson(A) >= +0.05` **and** the 95% interval lower
bound `> 0`.

**G2 — protein beats a matched wrong protein.**
`Pearson(X, correct) - Pearson(X, hard-wrong) >= +0.05` **and** lower bound `> 0`.

**G3 — calibration and ranking improve together.**
Against `A`: the equal-component target-mean MSE improves (delta `< 0`) with a
resolved interval, **and** at least one of Spearman / CI / sign accuracy
improves with a resolved interval, **and** none of Spearman / CI / sign accuracy
shows a *resolved degradation*.

**G4 — the controls cannot reproduce the gain.**
(a) `Pearson(D) - Pearson(A) <= 0.5 x [Pearson(X) - Pearson(A)]`;
(b) `Pearson(E) <= 0.10` and `Pearson(X) - Pearson(E) >= +0.05` with lower bound `> 0`.

**G5 — the high-confidence stratum does not reverse the effect.**
The `X - A` Pearson gain on the pooled (same-panel + cross-panel) bank and on
the same-panel bank have the **same sign**, and the same-panel gain is the one
that meets G1. If cross-panel carries the effect and same-panel reverses it, the
route fails regardless of the pooled number.

**G6 — activity cliffs are not sacrificed and easy pairs do not hide the result.**
(a) cliff-stratum sign accuracy of `X` `>= 0.50`;
(b) `sign_accuracy_cliff(X) - sign_accuracy_cliff(A) >= 0` and its interval is
not resolved-negative;
(c) the `X - A` Pearson gain computed on the cliff stratum alone `>= 0`.

**The route passes only if a single arm clears G1-G6 simultaneously.**

## 8. Stop rules

* **S1.** G1 or G2 fails -> the SAR-field family is **stopped**. It is not
  rescued by Cartesian tensors, conformer ensembles, extra attention layers, a
  few-shot adapter, a larger budget, more seeds, or any threshold change. The
  report must name the precise failure mechanism.
* **S2.** G4 fails (a control reproduces the gain) -> the measured gain is
  recorded as an artifact of capacity or of the training signal, not of protein
  information, and the family is stopped.
* **S3.** No arm may be retrained after its metrics are read. A bug found after
  reading metrics voids the affected numbers and requires a new preregistration.
* **S4.** Nothing is promoted to `model/` or `scripts/` unless the gate passes
  **and** a subsequent multi-seed confirmation is run. A single seed can reject
  the hypothesis; it cannot establish performance.

## 9. If and only if the gate passes

The next isolated experiment, in this order and not before:

1. deterministic offline RDKit conformer ensembles with stable seeds and
   recorded provenance;
2. a **ligand-intrinsic geometry** O(3) scalar/vector/rank-2 encoder against the
   same 2D SAR field — called ligand-intrinsic geometry, never bound pose and
   never atomic protein-ligand recognition;
3. only after geometry provides incremental evidence, a k=0/1/2/3/5 low-capacity
   few-shot field adapter.

## 10. Commands

```bash
python -m tools.research.stageS_sar_field.phase0_audit
python -m pytest tools/research/stageS_sar_field/tests -q
python -m tools.research.stageS_sar_field.train --arm A_ligand_only
python -m tools.research.stageS_sar_field.train --arm B_protein
python -m tools.research.stageS_sar_field.train --arm C_protein_cf
python -m tools.research.stageS_sar_field.train --arm D_protein_shuffled
python -m tools.research.stageS_sar_field.train --arm E_label_shuffled
python -m tools.research.stageS_sar_field.analyse
```
