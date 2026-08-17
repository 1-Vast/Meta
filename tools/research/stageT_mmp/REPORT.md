# Stage T — true-MMP transformation space: T0 measured, **T1 PASS**, **T2 FAIL**

Authorities: `PREREGISTRATION.md` (SHA-256 recorded in `T2_RESULT.json`),
`T0_RELIABILITY.json`, `T1_CENSUS.json`, `T1_COVERAGE.json`, `T2_RESULT.json`,
`runs/<arm>/RUN.json`, `runs/<arm>/*.rows.json`.
Structural verification: `tests/test_structural.py`, **29 passed** (`RUN_SLOW=1`).

**Verdict: the preregistered T2 gate fails, 3 of 10.** Per the frozen stop rule,
**protein-conditioned SAR latent space is formally stopped under the current
BindingDB protocol.** Nothing promoted; `model/` and `scripts/` unmodified; the
sealed split never mounted; the development-validation split never read.

---

## 1. What was actually tested

`delta_y(t, tau) = mu_tau + delta(t, tau) + noise`, and the estimand is the
**crossed double difference**

    D(tau, t1, t2) = delta_y(t1, tau) - delta_y(t2, tau)

which cancels the target-level offset **and** the generic chemical effect
`mu_tau` exactly. This is a new estimand, not a rescue of Stage S: it removes by
construction all three shortcuts that explained every earlier apparent positive
in this repository (target level, generic medicinal chemistry, target-identity
key).

`method_ladder/CLOSURE_MAP.md` family 8 recorded that pairs had never been
matched-molecular-pair identified, so the MMP constraint had never been
instantiated. It now has been, with RDKit's supported Hussain–Rea machinery.

## 2. Stage T0 — measurement reliability (not an MSE floor)

Both pre-aggregation artifacts were located and **hash-verified against the
corpus manifest** (`exact_labels.jsonl.gz` = `labels_sha256`,
`metadata_projection.jsonl.gz` = `projection_sha256`). **100% of `meta_train`
source rows recovered** (5,983).

Aggregation rule actually used, quoted from the manifest: `within_panel_
aggregation: median`, then `cross_panel_pair_aggregation: equal-panel median`,
endpoint `exact positive uncensored Ki`, `pKi = 9 - log10(Ki[nM])`.

| level | groups | median range | p95 | max | residual sd |
|---|---:|---:|---:|---:|---:|
| L1 same panel, same protocol | 99 (98 disagreeing) | **1.019** | 2.11 | 3.06 | 0.657 |
| L2 same panel, different protocol | **0** | — | — | — | not identifiable |
| L3 across panels | 133 (54 disagreeing) | 0.429 | 3.00 | 4.99 | 0.426 |

Derived difference-label variance under the actual aggregation rule:
**same-panel ≈ 0.858 pK², cross-panel ≈ 1.221 pK²**. That ordering is T0's job —
it is why cross-panel observations are the S3 weak stratum and are excluded from
the primary bank.

**What T0 does not establish.** Three things are stated in the artifact and
repeated here because they bound the claim:

* **L2 is not identifiable.** The projection's assay `protocol_sha256` never
  splits a `meta_train` panel, so technical replication cannot be separated from
  a paper reporting two conditions. L1 is *same-document repeat disagreement*,
  not pipetting noise, and no technical-replicate variance is reported.
* **59.4% of apparent cross-panel repeats have exactly zero range** — one
  physical measurement curated under two DOIs. The pooled L3 figure is deflated
  by them; the disagreeing-only figures are reported beside it.
* **96% of cells carry one source row**, so their error is unobservable, and the
  4.1% repeated subset is selected toward reference compounds. The corpus also
  dropped 333 conflicting ligands at admission, so this is a **lower bound** on
  raw disagreement.

**None of this is a universal benchmark MSE floor and it is not quoted as one.**

## 3. Stage T1 — the census **PASSES** all five thresholds

| threshold | required | measured | |
|---|---:|---:|---|
| same-panel fit observations | 2,000 | **37,945** | PASS |
| fit targets | 50 | **243** | PASS |
| exact keys spanning ≥3 targets **and** ≥3 components | 30 | **1,112** | PASS |
| internal observations | 300 | **5,035** | PASS |
| internal protein components | 10 | **25** | PASS |

Built from `rdMMPA.FragmentMol` single-cut fragmentation over 5,643 ligand
slots; 118 ligands admit no cut; 0 unparsable. Degree concentration is healthy
(top key 0.2% of observations, top-10 1.4%), the evidence graph is well
connected (largest component 97.6% of nodes), and 78.6% of exact keys are
singletons — so the identifiable core is the 1,112 repeated keys, not the bulk.

Two dimensions of the hypothesis are **not testable on this corpus**: stereo
edits (1 in fit, 0 in internal) and charge-changing edits (326 fit, 2 internal).

Fit↔internal reuse: 16.0% of internal exact keys also occur in fit.

**Double-difference sets:** fit **19,851** D-pairs over 4,621 keys / 162
components; internal **759** over 559 keys / **20 components**. Effective
independent units on the evaluation surface are therefore **20** — the binding
constraint on all power below.

## 4. Deployment-coverage audit

`C_k = P(at least one support-query pair forms a valid MMP)` on the frozen
nested banks, labels not read:

| k | C_k exact | component-equal | novelty-high | internal |
|---|---:|---:|---:|---:|
| 1 | 0.226 | 0.242 | 0.213 | 0.339 |
| 2 | 0.362 | 0.388 | 0.351 | 0.437 |
| 3 | 0.442 | 0.473 | 0.435 | 0.498 |
| 5 | **0.526** | 0.564 | 0.517 | 0.610 |

Coarse-key coverage is identical to exact on this corpus. So MMP reaches at most
**52.6%** of governed queries at k≤5: usable as a training signal and as a
partial inference mechanism, **not** a universal reference-based one, and no
artifact may present it as one.

## 5. Stage T2 — the gate **FAILS**, 3 of 10

Six matched arms, one code path, seed 20260820, 3,000 steps, identical batch
sequence, **no checkpoint selection**. 126,721 parameters each.

### The estimand verified itself

Arms **A (zero)** and **B (transformation-only)** produce **bitwise identical**
predictions (`D_hat ≡ 0`) and identical training loss (0.29818). A protein-free
model is structurally incapable of expressing `D` — confirming that `mu_tau` and
the target level cancel exactly, which is the whole point of the estimand.

### Headline, internal-validation components (759 rows, 20 components, EIU 20)

| arm | MSE | MAE | Pearson | Spearman | sign |
|---|---:|---:|---:|---:|---:|
| **A/B (predict D = 0)** | **0.6603** | **0.6232** | undefined (constant) | — | — |
| C correct protein | 1.5775 | 0.9580 | +0.0588 | +0.0414 | 0.5284 |
| D shuffled protein | 1.5650 | 0.9042 | +0.0402 | **+0.1434** | **0.5812** |
| E matched-wrong protein | 1.1759 | 0.8514 | −0.0822 | −0.0678 | 0.4148 |
| F labels shuffled in key | 0.8006 | 0.6743 | **+0.1836** | +0.1561 | 0.5694 |

**The trivial zero predictor is 2.4× more accurate than the protein-conditioned
model.** Predicting "no target-specific transformation response at all" beats
every trained arm on both MSE and MAE.

### The gate

| gate | result | |
|---|---|---|
| 1 correct − shuffled Pearson ≥ +0.05 | **+0.0186** | FAIL |
| 2 lower bound > 0 | [−0.2142, +0.2512] | FAIL |
| 3 correct − matched-wrong Pearson ≥ +0.05 | +0.1441 | pass |
| 4 lower bound > 0 | [−0.0840, +0.2908] | FAIL |
| 5 error **and** ranking both improve | MSE +0.0125 (worse), Spearman −0.1020, sign −0.0528 | FAIL |
| 6 label shuffle destroys the effect | C − F Pearson **−0.1249**; C − F MSE **+0.7770 [+0.2714, +1.3543] RESOLVED worse** | FAIL |
| 7 not confined to one transformation/component | one component = **4.66×** the whole effect; one key = 1.96× | FAIL |
| 8 protein shift aligned with truth | alignment +0.139, shift/truth sd **1.27** | pass |
| 9 no target-key shortcut contrast on fit-unsampled | −0.0071 | pass |
| 10 transformation-cold does not reverse | repeated +0.0733 vs disjoint **−0.0094** | FAIL |

`route_passes = false`.

## 6. The failure mechanism, precisely

**6.1 A wrong protein makes the model better.** Arm C evaluated with a
**shuffled** protein input scores Pearson **+0.1490** against **+0.0588** with
the correct one (paired, identical rows, only the protein input replaced). The
falsification control is not merely unresolved — it is **inverted**, the same
signature the A2 exact operator produced.

**6.2 Destroying the labels improves generalisation.** Arm F, trained on `D`
permuted *inside each transformation key*, reaches Pearson **+0.1836** on the
internal components against C's +0.0588, and C is **resolvedly worse** than F on
MSE (+0.7770 [+0.2714, +1.3543]) and MAE (+0.2838 [+0.1267, +0.4620]). A model
that learned nothing real generalises better than one that fitted the real
double differences. The fitted structure is therefore anti-transferable across
protein components.

**6.3 Both arms memorise the fit double differences equally.** On the
fit-unsampled bank (same targets, same keys, unseen rows) arm C reaches Pearson
**0.9121** and arm D — with protein identity **permuted across components** —
reaches **0.9192**, both at MSE ≈ 0.16. Out of component they fall to 0.0588 and
0.0402. The model learns per-(target, key) responses by rote, and a permuted
protein serves that purpose slightly better than the real one, because it is
still a unique per-target vector. Gate 9 passes only because the *contrast* is
≈ 0 there; the shortcut itself is enormous and shared.

**6.4 The +0.019 "effect" is one component.** Removing a single protein
component changes the effect by **4.66× its own magnitude**; removing one
transformation key changes it by 1.96×. There is no effect to attribute.

**6.5 The protein path is loud, not informative.** The protein-induced shift has
sd 1.025 against a truth sd of 0.804 — **127%** — with alignment +0.139. Stage P
and Stage S produced the same signature in two other function classes.

## 7. What this does and does not establish

* **T0** measured supervision reliability on a small, selected, partly
  duplicate-inflated subset. It is **not** an MSE floor.
* **T1** established that the transformation graph **is** identifiable: 1,112
  exact keys span ≥3 targets and ≥3 components, 37,945 same-panel fit
  observations, a well-connected evidence graph. **The hypothesis got a fair
  test — it was not blocked by coverage.** That is the most important thing T1
  contributes, because it removes "not enough MMP data" as an explanation.
* **T2** measured protein × transformation interaction on that identifiable
  graph and found none: below threshold, unresolved, with an inverted
  wrong-protein control, a label-shuffle control that beats the real model, and
  an effect smaller than the influence of any single component.
* None of the three is a zero-shot or few-shot DTA performance claim.

### Limitations, stated plainly

1. **One seed.** A single seed may reject; it cannot establish performance. The
   rejection here does not rest on one number: the zero predictor beating every
   arm, the inverted shuffled-protein control and the label-shuffle control all
   point the same way independently.
2. **Effective independent units = 20** (14 on the transformation-disjoint
   surface, 17 on the repeated-key surface). Intervals are wide, and a true
   effect below roughly ±0.20 Pearson could not be resolved here. What the data
   exclude is a *large* effect; several point estimates additionally sit in the
   wrong direction, which is not what an underpowered null looks like.
3. **Stereochemistry and charge changes are untestable** on this corpus (1 and 2
   internal observations). The stereochemical part of the hypothesis is
   **unmeasured**, not falsified.
4. The internal partition is protein-component-hard but not ligand-cold; MMP
   pairs are by construction chemically close, so the novelty caveat is stronger
   here than elsewhere.
5. T0's provenance step reads a single all-label artifact under a `meta_train`
   allow-list applied in the same parsing pass — **logical exclusion after
   parsing**, weaker than the model path's physical isolation, and disclosed as
   such in `T0_RELIABILITY.json`.

## 8. Decision

**The frozen stop rule fires. Protein-conditioned SAR latent space is stopped
under the current BindingDB protocol.** It is not rescued by DrugBAN,
PSICHIC-style attention, Cartesian tensors, conformers, MSA, more capacity, more
seeds, meta-learning or threshold changes.

The conditional next stage (transformation/edit tokens querying residue-level
protein regions) was gated on a pass and is **not authorized**; no part of it was
designed or run.

Two results survive and are worth carrying forward, neither protein-conditioned:

* the **MMP transformation graph is identifiable** on this corpus (T1), and
  MMP relations reach 52.6% of governed queries at k=5 (coverage audit);
* the **generic transformation effect `mu_tau` is the whole of the signal** — a
  protein-blind model is structurally zero on `D`, and predicting zero beats
  every protein-conditioned arm by 2.4× in MSE.

## 9. Verification

* `tests/test_structural.py`: **29 passed** with `RUN_SLOW=1`, including
  deterministic MMP decomposition, a positive control against silent empty
  fragmentation, inverse/sign consistency, attachment and stereochemistry
  preservation, no cross-target or cross-panel contamination, the physical
  meta-test seal, bank stability across `PYTHONHASHSEED` 0/1/12345, identity /
  antisymmetry / protein-cycle consistency, no dead trainable parameters, and
  parsed-AST checks that no label, `hash()`, or development-validation split
  name appears on any construction path.
* Two defects were found and fixed **before** any arm trained, both of which
  would have produced a false negative: `rdMMPA.FragmentMol` mis-binds four
  positional integers (and segfaults when fully named), and returning an RDKit
  `Atom` from a helper outlived its parent `Mol` and crashed on access. Both are
  now pinned by tests.
* One methodological repair was made **after** training and **before** reading
  any evaluation metric, and is disclosed: the per-arm batch seed included the
  arm name, so arms drew different training rows. It was changed to a shared
  sequence and all six arms retrained. Arms A and B then produced identical
  training loss, which is the expected consequence.
* Environment: conda env `drug`, Python 3.11.15, torch 2.6.0+cu124, CUDA
  available, RDKit 2023.09.6, numpy 1.26.4.

## 10. Commands

```bash
python -m tools.research.stageT_mmp.t0_reliability
python -m tools.research.stageT_mmp.t1_census
python -m tools.research.stageT_mmp.t1_coverage
RUN_SLOW=1 python -m pytest tools/research/stageT_mmp/tests -q
python -m tools.research.stageT_mmp.t2_train --arm A_zero
python -m tools.research.stageT_mmp.t2_train --arm B_transformation_only
python -m tools.research.stageT_mmp.t2_train --arm C_protein
python -m tools.research.stageT_mmp.t2_train --arm D_protein_shuffled
python -m tools.research.stageT_mmp.t2_train --arm E_protein_matched_wrong
python -m tools.research.stageT_mmp.t2_train --arm F_label_shuffled
python -m tools.research.stageT_mmp.t2_analyse
```
