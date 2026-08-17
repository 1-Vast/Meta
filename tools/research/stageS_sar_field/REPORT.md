# Stage S — cross-target protein-conditioned SAR field: **FAIL**

Authorities: `PREREGISTRATION.md` (SHA-256
`45eb15a318543b5bdb5529f5bde21a4ee791e1962991f6509f351fa67073b2bd`, recorded
inside `RESULT.json`), `PHASE0_AUDIT.json`, `RESULT.json`, `ATTRIBUTION.json`,
`runs/<arm>/RUN.json` and `runs/<arm>/*.rows.json`.
Structural gates: `tests/test_structural.py`, **24 passed**.

**Verdict: the preregistered gate FAILS.** Arm `B_protein` clears 1 of 6 gates,
arm `C_protein_cf` clears 0 of 6. Stop rules S1 and S2 both fire. The
SAR-field family is stopped. Nothing was promoted to `model/` or `scripts/`; the
sealed confirmation split was not mounted; the development-validation split was
not read at any point.

---

## 1. What was built, and what it guaranteed

`dy_hat(P,a,b) = V(P,L_b) - V(P,L_a)`, where `V` is a FiLM-modulated nonlinear
potential over a protein-free ligand coordinate `phi(L)` and protein response
coefficients `alpha(P)`. Because the prediction is the difference of a scalar
potential it is a curl-free field, so the three identities hold for **every**
parameter setting, before and after training. Measured, not assumed:

| identity | measured |
|---|---|
| `dy_hat(P,a,b) = -dy_hat(P,b,a)` | `torch.equal` — bitwise exact |
| `dy_hat(P,a,a) = 0` | `torch.equal` against zeros — bitwise exact |
| `dy_hat_ab + dy_hat_bc + dy_hat_ca = 0` | `< 1e-12` in float64 (rounding residual) |
| identities survive a parameter update | bitwise exact after an SGD step |

The rejected `e^T H e` term is absent structurally: the potential is a function
of **one** ligand coordinate and never receives a difference, so no even term
can exist. Protein conditioning happens exactly once — `phi` has no protein
argument, `alpha` has no ligand argument. No 3D coordinate, pose, docking result
or protein-ligand contact was constructed, and nothing in this stage may be
described as one.

## 2. Phase 0 — the eligible pair population

Governed BindingDB-Ki `main_v0` on the double-cold `v1` split, mounted through
the **physically isolated split view** (`seal_record()["isolation"]["level"] =
"physically_isolated"`, `evaluated = false`, 768 sealed cells withheld and
absent from the mounted surface). Frozen `scripts/internal_validation.py`
partition: **227 fit / 31 internal-validation** components.

| population | components | targets | pairs | same-panel | cross-panel | cliffs |
|---|---:|---:|---:|---:|---:|---:|
| fit | 227 | 290 | 138,695 | 42,195 (30.4%) | 96,500 | 3,372 |
| internal validation | 31 | 42 | 13,279 | 6,775 (51.0%) | 6,504 | 563 |

* **Assay context matters and is not pooled.** Same-panel `|dpK|` averages
  **0.943**; cross-panel **1.286**. Cross-panel differences carry inter-assay
  offset on top of chemistry, so they were **excluded from the training loss**
  and reported as their own stratum. Primary supervision is same-panel only:
  24,125 target-balanced pairs over 285 fit targets (`local` 8,865 / `medium`
  6,704 / `distant` 6,095 / `cliff` 2,461).
* **Target-balanced sampling was mandatory, not stylistic.** Pairs per target
  run median 36, max 26,335. Global pair sampling would have made the loss a
  report about three targets.
* Zero duplicated (target, ligand) rows, so `delta_y` is unambiguous.
* **Disclosure carried by every number here:** the internal partition is
  protein-component-hard but **not ligand-cold** — 183 of 638 (28.7%)
  internal-validation ligands are identity-shared with fit components and 51.4%
  share a Murcko scaffold. Absolute values in this stage are therefore **not**
  comparable with any k=0/k=5 figure on the double-cold protocol, and the
  chemical-novelty strata below are the honest read of transfer.

## 3. Arms and budget

One code path, seed 20260819, 4,000 steps, 128 pairs/step (16 targets x 8),
AdamW 3e-4 cosine, weight decay 1e-4, Huber(1.0) on signed `dpK` + 0.2 x
softplus sign term. **No checkpoint selection of any kind** — every arm trains
the same number of steps and the final parameters are evaluated, because the
reporting population is the internal-validation components and selecting on
them would leak.

| arm | params | final regression loss | elapsed |
|---|---:|---:|---:|
| `A_ligand_only` | 366,084 | 0.0364 | 342 s |
| `B_protein` | 523,524 | 0.0413 | 349 s |
| `C_protein_cf` | 523,524 | 0.0650 (+ 0.226 counterfactual) | 450 s |
| `D_protein_shuffled` | 523,524 | 0.0389 | 436 s |
| `E_label_shuffled` | 523,524 | 0.0732 | 342 s |

**The protein input does not even help in-sample.** `B` has 157,440 more
parameters than `A` and its final training regression loss is *higher*
(0.0413 vs 0.0364).

## 4. Result on the withheld internal-validation components

Primary bank: same-panel pairs, 4,270 pairs / 42 targets / 29 components.

| arm | MSE | equal-component target-mean MSE | Pearson | Spearman | CI | sign acc |
|---|---:|---:|---:|---:|---:|---:|
| `A_ligand_only` | 2.2541 | 2.6078 | **0.1599** | 0.1633 | 0.5553 | 0.5648 |
| `B_protein` | 2.2995 | 2.7102 | **0.1664** | 0.1765 | 0.5595 | 0.5615 |
| `C_protein_cf` | 4.0894 | 7.7976 | 0.0739 | 0.0897 | 0.5298 | 0.5420 |
| `D_protein_shuffled` | 2.2583 | 3.2794 | **0.1684** | 0.1788 | 0.5603 | 0.5581 |
| `E_label_shuffled` | 2.8772 | 3.7549 | −0.0355 | −0.0324 | 0.4891 | 0.4882 |

Matched-wrong protein (same rows, same arm, **only the protein input replaced**):

| arm | Pearson correct | Pearson hard-wrong | MSE correct | MSE hard-wrong |
|---|---:|---:|---:|---:|
| `B_protein` | 0.1664 | 0.1226 | 2.2995 | 2.2303 |
| `C_protein_cf` | 0.0739 | −0.0181 | 4.0894 | **12.6204** |
| `D_protein_shuffled` | 0.1684 | 0.1933 | 2.2583 | 2.1508 |

### The gates

| gate | `B_protein` | `C_protein_cf` |
|---|---|---|
| **G1** protein beats ligand-only (Pearson >= +0.05, lo > 0) | **FAIL** +0.0065 [−0.0761, +0.0865] | **FAIL** −0.0860 [−0.2258, +0.0292] |
| **G2** protein beats matched-wrong protein | **FAIL** +0.0438 [−0.0047, +0.1124] | **FAIL** +0.0920 [−0.0270, +0.1990] |
| **G3** MSE and ranking improve together | **FAIL** eq-comp MSE +0.1025 (worse, unresolved); 0 resolved ranking gains | **FAIL** eq-comp MSE **+5.1898 [+2.0973, +9.1305] RESOLVED degradation** |
| **G4** controls cannot reproduce the gain | **FAIL** shuffled protein +0.0085 >= candidate +0.0065 | **FAIL** |
| **G5** same-panel does not reverse | **FAIL** same-panel +0.0065 vs pooled −0.0037 (opposite signs; both below threshold) | **FAIL** |
| **G6** cliffs not sacrificed | PASS (cliff sign 0.5827 vs 0.5680, +0.0147 [−0.0884, +0.0820]) | **FAIL** cliff sign 0.5294 |

`route_passes = false`.

### Strata (Pearson on the primary bank)

| stratum | n | A | B | C | D | E |
|---|---:|---:|---:|---:|---:|---:|
| all | 4,270 | +0.1599 | +0.1664 | +0.0739 | +0.1684 | −0.0355 |
| cliff | 544 | +0.1023 | +0.2339 | −0.0467 | +0.0980 | −0.1022 |
| local | 1,726 | +0.0958 | +0.0853 | +0.0368 | +0.1355 | −0.0186 |
| medium | 1,314 | +0.1747 | +0.1384 | +0.0865 | +0.1633 | −0.0702 |
| distant | 686 | +0.2565 | +0.2779 | +0.1756 | +0.2754 | +0.0503 |
| novelty low | 1,696 | +0.0767 | +0.2332 | +0.0759 | +0.1374 | −0.1355 |
| novelty mid | 1,253 | +0.0594 | −0.0464 | −0.0283 | +0.0556 | −0.0786 |
| **novelty high** | 1,321 | **+0.4003** | +0.3335 | +0.1596 | +0.3655 | +0.1310 |

The `B − A` cliff gain (+0.1316 [−0.0616, +0.2631]) and the `B − A`
high-novelty loss (−0.0668 [−0.1673, +0.0626]) are both unresolved. The one
consistent reading is that **the ligand-only field is the best arm on the most
chemically novel pairs**, which is where a transferable mechanism would have to
earn its keep.

## 5. The precise failure mechanism

This is not the "inert operator" failure of the A2 exact probe, and it is not a
dead gradient. It is the Stage P failure mode — **loud and misaligned** —
reproduced in a completely different function class, plus one new mechanism
that the controls isolate.

### 5.1 The protein path is fully alive and points nowhere

`ATTRIBUTION.json`, measured on the trained fields over the 4,270 primary pairs
(label spread 1.1288 pK):

| arm | protein-induced spread (pK) | / label spread | shift alignment with truth | response-vector pairwise cosine |
|---|---:|---:|---:|---:|
| `A_ligand_only` | 3.6e-08 | 3.2e-08 | n/a (constant response) | 1.000 |
| `B_protein` | **0.7044** | **0.624** | **+0.0926** | 0.245 |
| `C_protein_cf` | **1.9843** | **1.758** | +0.0593 | 0.223 |
| `D_protein_shuffled` | 0.3437 | 0.304 | −0.0747 | 0.246 |
| `E_label_shuffled` | 0.5922 | 0.525 | −0.0918 | 0.282 |

Swapping the protein moves `B`'s prediction by **62% of the entire label
spread** — against the A2 exact operator's 0.0028 pK of query-specific content.
The mechanism works. Its alignment with the signed truth is **+0.093**, and its
alignment with the part of the truth the ligand-only field leaves unexplained is
**−0.051**. Arm `A`'s constant response is genuinely constant (spread 3.6e-08,
pairwise cosine 1.000), so the contrast is clean.

**Conclusion: the field learns a large, protein-specific response that carries
essentially no signed-affinity information.**

### 5.2 The protein encoder is used as a target key, not as a response coordinate

The training-health bank (same-panel fit pairs the training draw did not take —
same targets, unseen pairs, diagnostic only, never used for selection):

| arm | fit-unsampled Pearson | internal-validation Pearson |
|---|---:|---:|
| `A_ligand_only` | 0.7190 | 0.1599 |
| `B_protein` | **0.8282** | 0.1664 |
| `D_protein_shuffled` | **0.8205** | 0.1684 |
| `E_label_shuffled` | −0.0693 | −0.0355 |

The protein input buys **+0.109 Pearson in-distribution** and **+0.0065 out of
component**. And a protein whose identity has been **permuted across
components** buys the same +0.108 in-distribution. A permuted protein is still a
unique per-target vector, so what the potential learned to use is the *key*, not
the *biology*. That is the whole of the measured protein contribution, and it is
by construction non-transferable to an unseen homology component.

This also explains the one contrast that does resolve in `B`'s favour: correct
minus hard-wrong equal-component target-mean MSE **−0.4123 [−0.8708, −0.0104]**.
The hard-wrong donor is drawn from the **fit** components, so it injects a
strongly fitted, target-specific response into a target it does not belong to;
the correct (unseen) protein lands in a region of response space the model never
fitted and behaves closer to a mean. A resolved wrong-protein degradation is
therefore **not** evidence of transferable protein information here, and this
stage recommends that any future wrong-protein control be paired with a
shuffled-protein arm before it is read as one.

### 5.3 The counterfactual loss found the degenerate solution

Arm `C` satisfies "the correct protein must fit better than a hard wrong one" by
making the **wrong** side explode, not the correct side better: hard-wrong MSE
**12.6204** against its own correct MSE 4.0894, and a correct-vs-wrong
equal-component MSE gap of **−22.93 [−50.82, −3.17] RESOLVED**. The
preregistered hinge (saturating at 0.25 pK²) was intended to bound exactly this
and did not: bounding the *reward* does not remove the *direction*, because the
term only stops pushing once the gap is met and any parameter change that
inflates the wrong branch reaches the margin fastest.

The cost is a resolved regression on the real task: same-panel MSE
**+1.8353 [+0.3116, +4.5365]**, equal-component MSE **+5.1898 [+2.0973,
+9.1305]**, and on the high-novelty stratum Pearson **−0.2407 [−0.4433,
−0.0121]**, Spearman **−0.1683**, CI **−0.0619**, sign **−0.0666**, every one
resolved. Protein-induced spread rose to 1.98 pK — **176% of the label
spread** — while alignment with truth stayed at +0.059. The counterfactual
objective maximised protein *sensitivity* without producing protein
*informativeness*.

The B/C split does what it was designed to do: representation gain and
training-objective gain are separable, and both are measured at approximately
zero and clearly negative respectively.

### 5.4 The label control is clean

`E_label_shuffled` (within-target permuted labels, target level preserved
exactly) scores Pearson **−0.0355** on the primary bank against `A`'s +0.1599 —
resolved degradation at −0.1954 [−0.2992, −0.0960], and −0.3032 [−0.4077,
−0.1719] pooled. It also fails on the fit-unsampled bank (−0.0693). There is no
label leakage and no route by which the geometry of the bank alone produces the
signal.

## 6. Answer to the hypothesis

> Does the correct protein carry information about the signed within-target
> affinity difference between two ligands, beyond what a protein-blind field
> extracts from chemistry alone?

**Not measurably, on this protocol, with these inputs, at this budget.** The
protein-conditioned field is +0.0065 Pearson over the ligand-only field, a
permuted protein reproduces that difference, and the protein-induced movement is
large and truth-unaligned.

What the signed-difference formulation *did* deliver is the ligand-only result:
a conservative potential field trained directly on signed same-panel `dpK`
transfers to unseen protein components at Pearson **+0.1599**, Spearman
**+0.1633**, sign accuracy **0.5648**, and **+0.4003 Pearson on the most
chemically novel tercile**. That is consistent with, and independent of, the
Stage L2 frozen-probe finding of a **protein-independent** directional SAR
signal — reached here through end-to-end training of a different function class
rather than a linear probe on frozen features. It is a ligand-side result and
must not be reported as protein-conditioned meta-learning.

## 7. Limitations, stated plainly

1. **One seed.** The preregistration says a single seed may reject but cannot
   establish performance, and that is the correct reading of G1 in isolation:
   the interval [−0.0761, +0.0865] does **not** exclude a true gain as large as
   +0.087. What the single seed does establish is (a) the preregistered gate
   fails, (b) the shuffled-protein arm reproduces the observed difference in a
   paired comparison on identical rows, and (c) the protein path moves the
   prediction by 0.70 pK with +0.093 alignment. (b) and (c) are the load-bearing
   evidence, not the G1 point estimate.
2. **Heavy in-distribution fit.** Fit-unsampled Pearson 0.72–0.83 against
   internal-validation 0.16 is a large generalisation gap. A different budget or
   stronger regularisation would change the absolute numbers. It would have to
   change the *relative* protein contribution to change the verdict, and the
   shuffled-protein control is evidence that it would not.
3. **The internal partition is not ligand-cold** (28.7% shared ligand identity,
   51.4% shared scaffold). Absolute values are not comparable with the
   double-cold protocol; the novelty strata are the transfer read.
4. **Scope.** This is one architecture (FiLM-modulated potential over a 2D-graph
   + Morgan coordinate and a frozen ESM-2 150M sequence readout), one training
   objective family, one seed, one budget, on BindingDB-Ki. It is not a theorem
   about protein-conditioned SAR models, and it does not extend to other
   datasets, other protein representations or untested architectures.
5. **Two candidate arms were gated**, which doubles the multiplicity. Both are
   reported in full; neither passed anything close to threshold, so the
   multiplicity does not change the verdict.
6. `C`'s degenerate solution is a property of *this* counterfactual loss form.
   A different relational loss might avoid it — but under the stop rules that is
   not a rescue this stage is permitted to attempt.

## 8. Decision

**Stop rules S1 (G1/G2 fail) and S2 (a control reproduces the gain) both fire.
The SAR-field family is stopped.**

It is **not** rescued by Cartesian tensors, conformer ensembles, extra attention
layers, a few-shot adapter, a larger budget, more seeds, or any threshold
change. The conditional geometry experiment in section 9 of the preregistration
(offline RDKit conformer ensembles, a ligand-intrinsic O(3) encoder, then a
low-capacity few-shot field adapter) was gated on a pass and is **not
authorized**; no part of it was designed or run.

Nothing was promoted. `model/` and `scripts/` are unmodified. The sealed
confirmation split was not mounted (physically isolated view, 0 evaluations).
The development-validation split was not read at any point — a parsed-AST test
fails the suite if its name appears as a string constant anywhere in this stage.

## 9. Verification

* `tests/test_structural.py`: **24 passed** (`RUN_SLOW=1`, 52.6 s), including
  the cross-process bank-stability check under `PYTHONHASHSEED` 0/1/12345.
* Maintained suite `python main.py verify tests`: **310 passed / 6 skipped**,
  unchanged from the repair commit.
* Environment: conda env `drug`, Python 3.11.15, torch 2.6.0+cu124, CUDA
  available (RTX 4060 Laptop), rdkit 2023.09.6, numpy 1.26.4.

## 10. Commands

```bash
python -m tools.research.stageS_sar_field.phase0_audit
RUN_SLOW=1 python -m pytest tools/research/stageS_sar_field/tests -q
python -m tools.research.stageS_sar_field.train --arm A_ligand_only
python -m tools.research.stageS_sar_field.train --arm B_protein
python -m tools.research.stageS_sar_field.train --arm C_protein_cf
python -m tools.research.stageS_sar_field.train --arm D_protein_shuffled
python -m tools.research.stageS_sar_field.train --arm E_label_shuffled
python -m tools.research.stageS_sar_field.analyse
python -m tools.research.stageS_sar_field.diagnose
```

Per the repository's stage-artifact rule the five failed checkpoints
(`runs/<arm>/field.pt`) were deleted after `ATTRIBUTION.json` was written. The
prediction rows, run records and both analyses are retained; every checkpoint is
reproducible from the recorded seed, budget and commands above.
