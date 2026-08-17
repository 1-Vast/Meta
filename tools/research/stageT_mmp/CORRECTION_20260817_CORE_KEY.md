# Scoped correction to Stage T — the transformation key omitted the shared core

**This document does not alter Stage T's `PREREGISTRATION.md`, its
machine-readable results (`T0_RELIABILITY.json`, `T1_CENSUS.json`,
`T1_COVERAGE.json`, `T2_RESULT.json`), its prediction rows or its `REPORT.md`.
Those remain exactly as recorded.** It scopes what they support.

Forensic authority: `tools/research/stageV_core_mmp/STAGE0_FORENSICS.json`,
produced by `tools/research/stageV_core_mmp/stage0_forensics.py`, which reads
Stage T's own construction code unmodified.

## 1. The defect

The requested estimand was

    tau = (shared core, R_a -> R_b, attachment context, stereochemistry, charge)

Stage T implemented

    exact_key = f"{attachment_context}|{R_a}>>{R_b}"          # mmp.py:189

with **no shared core**. It then

* aggregated every observation of one `(exact_key, target)` with a median
  (`t2_dataset.py:93`), pooling different scaffolds into one "target effect";
* built the model's transformation descriptor from `R_a`, `R_b`, the attachment
  context, charge and a stereo flag — **the core is never read**
  (`mmp.py`, `descriptor`).

Stage T's `REPORT.md` states that the crossed double difference "cancels the
target-level offset **and** the generic chemical effect `mu_tau` exactly". That
statement is true only when both targets realise the *same* tau. Under a
core-blind key,

    D = [mu_core1 - mu_core2] + [delta(t1,tau) - delta(t2,tau)] + noise

and the first bracket is a generic chemical residual with no protein content.

## 2. How large the residual actually is

The decisive measurement holds the **protein fixed**: for one target realising
one Stage T key on two or more cores, the spread of `delta_y` across cores is
pure generic chemical context.

| statistic | value |
|---|---:|
| within-target across-core absolute gap, median | **0.269 pK** |
| mean | 0.415 |
| p95 | **1.268** |
| max | 3.401 |
| within-target across-core range, median | 0.454 |
| range p95 | 1.549 |

For scale, the truth `D` on Stage T's internal evaluation bank has sd
**0.804 pK**. The nuisance term is the same order as the entire quantity being
modelled.

How often it bites:

| | fit same-panel | internal same-panel |
|---|---:|---:|
| Stage T keys spanning >1 core | 1,848 / 26,695 (6.9%) | 290 / 3,402 (8.5%) |
| max cores under one key | **57** | 22 |
| observations under a multi-core key | 7,130 (**18.8%**) | 990 (**21.6%**) |
| target effects median-pooling >1 core | 1,801 (5.1%) | 279 (6.9%) |
| **D rows whose two targets have disjoint core sets** | **8,010 / 19,851 (40.4%)** | **219 / 759 (28.9%)** |

**40.4% of Stage T's training rows and 28.9% of its evaluation rows compared two
targets with no core in common.** For those rows the exact-cancellation claim is
false, and the contamination has median 0.27 pK and p95 1.27 pK.

## 3. What Stage T still validly establishes

1. **The concrete Stage T discriminator failed, and that measurement stands.**
   A pooled-protein MLP over a core-blind edit descriptor, trained on `D`, was
   beaten 2.4x on MSE by the trivial zero predictor (0.6603 vs 1.5775), its
   shuffled-protein substitution *improved* it (+0.1490 vs +0.0588 Pearson), and
   a label-shuffled arm generalised better (+0.1836). Those are real properties
   of that implementation.
2. **The defect supplies a mechanism for that failure rather than excusing it.**
   If 40% of the training target is contaminated by a core-mismatch nuisance of
   the same magnitude as the signal, then fitting `D` is substantially fitting
   noise — which is precisely why a model trained on destroyed labels
   generalises better and why the zero predictor wins. Stage T's own numbers are
   consistent with the defect.
3. **T0 stands unchanged.** Provenance recovery, the aggregation rule, the
   reliability levels and the non-identifiability statements do not depend on
   the key.
4. **T1's coverage direction stands.** MMP transformations exist in quantity;
   deployment coverage `C_k` = 0.226/0.362/0.442/0.526 is computed from the MMP
   *relation*, not from the key, and is unaffected.

## 4. What Stage T does NOT establish, and is hereby withdrawn

* The claim that `D` cancels `mu_tau` **exactly** in Stage T. It does not, for
  40% of training rows.
* **The global closure claim.** Stage T's `REPORT.md` section 8 and the
  corresponding lines in `task.md`, `history.md` and `report/EVIDENCE_LEDGER.md`
  stated that "protein-conditioned SAR latent space is formally closed under the
  current BindingDB protocol". **That is withdrawn.** Stage T tested one
  coarsened-key, pooled-protein discriminator. It cannot close a representation
  family it did not test.
* The reusability of the figure "1,112 exact transformation keys spanning >=3
  targets and >=3 components". Under the core-inclusive key the corresponding
  count is **1,001** and every downstream count changes; the 1,112 figure must
  not be quoted for the corrected estimand.

## 5. What remains genuinely unmeasured

* The requested core-inclusive estimand, with any protein operator.
* A **local** (non-pooled, non-target-ID) protein-region operator of the kind
  the active contract's Phase 1 requires. Stage T's discriminator concatenated a
  pooled ESM vector with an edit vector; it never had residue-region tokens.
* Stereochemical edits (1 fit / 0 internal observations) and charge-changing
  edits (326 fit / 2 internal). Unmeasured, not falsified.

## 6. Effect of the corrected key on the census

| | Stage T core-blind key | requested core-inclusive key |
|---|---:|---:|
| fit same-panel observations | 37,945 | 37,945 |
| fit keys | 26,695 | 30,463 |
| fit keys with >=3 targets and >=3 components | 1,112 | **1,001** |
| fit `D` rows | 19,851 | **12,740** |
| fit `D` components / EIU | 162 / 162 | 99 / **99** |
| internal `D` rows | 759 | **546** |
| internal `D` components / EIU | 20 / 20 | 10 / **10** |
| internal `D` rows whose key is repeated in fit | 270 | **32** (4 components) |
| internal `D` rows whose key is absent from fit | 489 | **514** (7 components) |

The fit side remains rich. The **evaluation** side does not: requiring a
complete chemical context that repeats across targets *inside the withheld
protein components* leaves 32 rows over 4 components on the repeated-key
surface.

## 7. Consequence

Stage T is reclassified from "protein-conditioned SAR latent space closed" to
**"the coarsened-key pooled-protein discriminator is rejected"**. Phase 1 of the
active contract returns to **pending**, and the corrected test is carried out
under a new preregistration in `tools/research/stageV_core_mmp/`.

`task.md`, `history.md` and `report/EVIDENCE_LEDGER.md` are corrected
accordingly; Stage T's own frozen artifacts are left untouched and this document
is linked from each of them.
