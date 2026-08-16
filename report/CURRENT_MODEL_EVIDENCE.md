# Current Model and Experimental Evidence

Last consolidated: 2026-08-16 (R5-R13 cycle complete; record audited)

This is the single narrative authority for the active MetaSieve model. Numeric
authority remains the corresponding leaf `RESULT.json`. Historical failures are
retained as machine-readable artifacts and summarized here rather than in
multiple overlapping reports.

The consolidated end-state of the R0-R13 chain — the corrected k=0 Pareto
frontier, the scoped activity-cliff record, and the `meta_test` separation —
is `report/BOUNDARY_20260816.md`.
`python -m scripts.audit_research_record` recomputes all of it from the leaf
artifacts, and `tests/test_research_record.py` fails if this narrative drifts
from them.

## Reading the word `meta_test` in this document

Two different populations carry that name here.

* Stages 4, 6 and 7 report `bindingdb_ki_main_v0` `meta_test` — 42 targets /
  6-7 components, under the **older, protein-only-cold** protocol. It was
  legitimately used and is **consumed**.
* From Stage R1 onward, `meta_test` means
  `bindingdb_ki_double_cold_v1` — 22 targets / 10 components. It is
  **physically sealed and has never been opened**
  (`QPSMPData(include_meta_test=False)`).

A `meta_test` result in a Stage 4/6/7 section is therefore *not* evidence of
generalization under the double-cold protocol.

## Evidence grades

Every claim below carries one of four grades. They are not interchangeable.

| grade | meaning |
|---|---|
| **exploratory** | measured once; no preregistration; may use oracle or transductive quantities |
| **development** | preregistered gates, but hyperparameters selected on the same population used for inference |
| **conditional** | resolved, but conditional on a stated unresampled factor (trained seeds) or an uncontrolled population property (ligand overlap) |
| **confirmed** | selection separated from inference by an outer fold, preregistered gates, positive component-level lower bound on a population never used for selection |

**No result in this project is currently `confirmed`, and the Stage 10 k=0
result is withdrawn.** Stage R0 falsified it: under leave-one-component-out
selection the exact-ligand-free effect is **-0.217 [-0.785, +0.261]**, i.e. the
prior makes genuinely new ligands *worse*, and the whole benefit sits in the
near-duplicate and exact-overlap strata (1.3581 -> 1.0114). See
`meta_fewshot/stageR0_retrieval_falsification_20260815/REPORT.md`.

The Stage R0 audit
(`meta_fewshot/stageR0_retrieval_falsification_20260815/AUDIT_VERIFICATION.md`)
verified eleven binding corrections by recomputation; all eleven hold. The two
that change the standing conclusions:

* the Stage 10 k=0 result is **development, conditional on ligand overlap**.
  48.9% of its query cells contain a ligand present verbatim in `meta_train`, and
  on the 12 exact-ligand-free targets the reduction falls from +0.198
  [+0.016,+0.405] to **+0.050 [-0.074,+0.175]**, i.e. unresolved. 6 of 10
  components improve, not all;
* Stage 9's "protein retrieval is weak" is **scope-limited to raw pooled ESM
  cosine**, whose similarities occupy a band of width 0.21 around 0.90 with a
  0.024 spread across the nearest 16 training targets — near-uniform softmax
  weights by construction. Train-only centring widens that spread to 0.238.
  Protein representation for k=0 is reopened; Mac-Diff locality, conformer
  routing, PBCNet2.0 and Cartesian equivariance are **not** reopened, having been
  rejected on structural-input coverage and multi-seed training evidence.

## Task and protocol

- BindingDB Ki, strict cold target, CD-HIT40 component-hard split.
- Nested support k=0/1/2/3/5 with unique, disjoint support/query ligands.
- Fixed evaluation bank (`evaluation_seed=73101`) and equal-component evidence.
- Query labels are metrics/loss targets only; target IDs are lookup keys only.
- Normal single-stage episodic training; no ridge, closed-form adaptation,
  inner loop, or test-time gradients.
- The active bank has no common-frame complex coordinates. Current performance
  is sequence+2D evidence, not atomic 3D recognition.

### Two banks are now reported

The frozen protocol bank uses `eval_targets_per_component=1` and contains **6
episodes per k**. It remains the retained comparator. Because the Stage 0 audit
showed it cannot resolve differences below about 0.05 MSE, every result is now
also reported on a **wide bank** built identically over all 42 eligible
meta-test targets across 7 components.

## Active model

`model/interaction_grammar.py` (`--arch grammar`) is the active candidate:

1. residue encoder over pooled ESM slots with a chemistry bias;
2. atom-to-residue cross attention feeding a globally shared contact-type
   dictionary;
3. one interaction embedding that serves as both the zero-shot readout and the
   few-shot kernel key;
4. `f = f0(P,Lq) + shrink(n) * sum_k softmax_k(sim) * rho(q,k) * r_k` with
   `r_k = y_k - f0(P,L_k)` and `rho in (0,2)`.

Because `rho` depends on the query, a single support observation already
produces a query-specific correction; the retained model's k=1 correction was
identically a scalar. `rho == 1` with flat weights recovers the shrunken support
mean, so level-only abstention remains inside the hypothesis class.

`model/qpsmp_meta.py` (`--arch bpsf`) is retained unchanged as the control arm.

## Headline result

Frozen protocol bank, three seeds each, MSE in pK squared:

| configuration | k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---:|---:|---:|---:|---:|
| retained no-warmup baseline | 2.115 | 1.704 | 1.341 | 1.270 | 1.253 |
| **interaction grammar** | **1.786** | **1.493** | **1.155** | **1.077** | **1.026** |
| reduction | -15.6% | -12.4% | -13.9% | -15.2% | -18.1% |

Every seed of the new model is below the best retained-baseline seed at every k.
Wide bank, three-seed mean: 3.246 / 1.977 / 1.642 / 1.493 / 1.271 against
3.589 / 2.386 / 1.993 / 1.831 / 1.581 for the retained baseline (one seed).

**The governed admission gate did not pass.** See "Stage 4" below. The
improvement is attributed to the zero-shot trunk plus target-level calibration,
not to query-specific SAR transfer.

## Stage 0 audit of the retained baseline

Artifacts: `meta_fewshot/stage0_audit_20260815/`.

1. **Confirmed bug — dead trainable branches.** 303,362 of 3,788,937 parameters
   (8.0%) received exactly zero gradient from the query loss at every k:
   `meta.ligand_baseline.*`, `meta.protein_level.*`, `meta.term.state.weight`,
   `pair_section.latent.interaction.response_weight`,
   `pair_section.latent.section.weight`, `protein_encoder.proj.*`,
   `protein_encoder.norm.*`. `qpsmp_meta.py` also bound `query_section` and
   never used it.
2. **Confirmed — degenerate endpoint.** Zero-shot spread across the queries of
   one episode was 0.065 pK against a 0.93 pK label spread; the mechanism
   dictionary was a constant 0.4913 with cross-episode std 0.00029.
3. **Confirmed — protein blindness.** Swapping in a cross-component donor
   protein moved the zero-shot output by 0.0093 pK and cost 0.012 MSE.
4. **Architectural — k=1 is scalar by construction**, verified at gradient
   level: `local_scale_logit`, `log_temperature` and `interaction_key.*` have
   exactly zero gradient at k=1.
5. **Architectural — the level path is gradient-inert.** In
   `QPSMPMetaLearner.infer` the two `level_adjustment` terms cancel identically
   and the surviving term is built from a detached residual, so the few-shot
   loss delivers no gradient to the endpoint on support ligands.
6. **Optimization — not merely undertrained.** A 400-step probe (4x budget)
   improved the validation score to 1.595 while worsening every frozen-bank test
   metric (k=0 2.047 -> 2.308) and reaching 8,399 MB on an 8,188 MB device.
7. **Falsified: a capacity explanation.** On a synthetic protein-by-ligand
   contact-type bilinear task the retained trunk reaches relative held-out MSE
   0.003-0.008 at lr 1e-3 and 3e-4, collapsing to a constant only at 3e-3 and
   above. Both trunks can express the interaction; the retained one is 4.4x more
   expensive per step with a narrower stable learning-rate band.
8. **Confirmed performance bug — episode materialization.** The compact ligand
   bank kept a one-shard LRU and reloaded 0.62 MB npz archives dozens of times
   per episode: 1,158 ms per episode, roughly 80% of every training step. With
   all shards resident it is 15.2 ms. Training-step cost fell from 5.57 s to
   0.75 s (`bpsf`) and 0.205 s (`grammar`). Numerics are unchanged.
9. **Missing diagnostic — evaluation power.** On the wide bank the retained
   baseline scores k=0 3.589, which is **worse than the meta-train global-mean
   constant (3.441)**, and k=5 1.581, which is **worse than the plain support
   mean (1.523)**. The frozen 6-episode bank is a favourable subsample.
10. **Missing diagnostic — reproducibility.** Re-running the identical seed and
    configuration diverged from step 40 onward; CUDA training here is not
    bitwise deterministic.

Label-only references, wide bank: global mean 3.441, ligand-average prior 3.119,
support mean 2.346 / 2.180 / 1.918 / 1.523 at k=1/2/3/5, oracle target-mean
level ceiling 1.100.

## Stage 1: held-out synthetic gates

`tests/test_interaction_grammar_synthetic.py`, 13 gates, all passing: exact k=0
identity, query-specific **and gradient-trainable** k=1 label effect, support
permutation invariance, query permutation equivariance, level-only abstention
(zero residual gives exactly zero transport), shared-mechanism recovery,
single-support shared-mechanism recovery, private-mechanism rejection, linearity
in the support labels, no query-label input, no dead trainable branch, and a
protein-conditioned zero-shot trunk gate.

## Stage 2: matched-budget architecture discriminator

`meta_fewshot/stage2_grammar_discriminator_20260815/`. Seed 20260812, 60 steps,
one changed variable. Wide bank: 3.996 / 2.204 / 1.829 / 1.616 / 1.372 against
4.085 / 2.631 / 2.146 / 1.866 / 1.611. Peak memory 1,698 MB against 6,053 MB.
Two preregistered gates failed at this budget: zero-shot spread 0.107 pK
(threshold 0.20) and wide-bank wrong-protein zero-shot gap -0.019 (threshold
0.05). Both were carried forward as blocking gates for zero-shot claims.

## Stage 3: budget, schedule and diagnosed capacity

`meta_fewshot/stage3_scaled_budget_20260815/`. 2000 steps, cosine schedule,
lr 6e-4, `backbone_lr_scale` 1.0, validation bank widened to 2 targets per
component (the frozen test bank is untouched; `val_targets_per_component` is a
separate field).

| arm | architecture | parameters | wide k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---|---:|---:|---:|---:|---:|---:|
| A | `bpsf` | 3.79M | 3.720 | 2.407 | 2.023 | 1.839 | 1.597 |
| B | `grammar` | 1.82M | 3.080 | 1.875 | 1.583 | 1.469 | 1.244 |
| C | `grammar` | 7.29M | 3.021 | 1.955 | 1.641 | 1.461 | 1.241 |

* **The budget and schedule alone bought nothing.** Arm A received the identical
  2000-step cosine schedule and the identical data-path fix and matched its own
  100-step version; its zero-shot spread *fell* to 0.0065 pK. Peak allocation
  9,110 MB exceeded the device.
* **Capacity 1.8M -> 7.3M bought no consistent MSE gain** (B wins k=1, k=2;
  C wins k=0, k=3, k=5). Its one reproducible benefit is control sign: arm B's
  permutation control is inverted at k=2, 3 and 5; arm C's is not. Recorded as a
  negative result for width.
* Carried-forward R2 now passes: wrong-protein zero-shot gap 0.674.
  R1-secondary still fails: spread 0.186-0.196 pK against 0.20.

## Stage 4: three-seed governed admission run

`meta_fewshot/stage4_grammar_admission_20260815/`. Seeds 20260812/13/14,
7,294,171 parameters, 8,000 episodes each, 6,507-6,863 MB peak, 42 meta-test
targets over 6 components, 9,999 paired-component bootstrap draws.

| # | gate | outcome |
|---|---|---|
| 1 | k=0 does not regress | **pass** |
| 2 | k=1 gain query-specific | point estimate only |
| 3 | k=2,3,5 improve in every seed | **pass** |
| 4 | `full` beats `sar_cut` with positive bootstrap lower bound | **fail** at every k |
| 5 | correct support beats permuted and matched wrong support | wrong support **pass**; permutation **fail** at k=2,3,5 |
| 6 | CI and Spearman improve with MSE | **fail** (CI 0.647 -> 0.571-0.610; Spearman 0.372 -> 0.169-0.257) |
| 7 | checkpoint re-evaluation reproduces | **pass**, exactly |
| 8 | no dead trainable branch | **pass** (0 at k>=2; k=0/k=1 zeros are semantic) |

**Admission is refused.** The supported conclusion is that the improvement comes
from a genuinely protein-conditioned zero-shot trunk plus target-level
calibration. The transferability gate produces a small squared-error gain that
does not survive a six-component bootstrap, is net negative at k=5, and reduces
within-target ranking quality at every k. Permuting support labels leaves
`mean(r)` exactly unchanged, so the permutation contrast isolates the
query-specific channel — and it is negative at k=2, 3 and 5.

## Stage 5: geometry audit and signed relative transport

`meta_fewshot/stage5_relative_transport_20260815/`.

### Geometry is closed on data

`scripts/audit_geometry_coverage.py`: **0 of 17,717** BindingDB deployment cells
have a common-frame protein-ligand complex. 15/499 targets have an exact holo
sequence but always bound to a different ligand; 84/9,880 ligands share a holo
SMILES; 110/499 targets match at containment level. The processed structure
assets store contact maps and distance bins, not coordinates.

Every Cartesian/equivariant family (PBCNet2.0, TensorNet, PaiNN, MACE,
Equiformer, SE(3)-EGNN) fails on this single input constraint, so none is
applicable to the deployment path. `model/cartesian.py` stays verified by
`tests/test_cartesian.py` and unused; the active models raise on coordinate
inputs. Reopening this requires new data, not a new architecture.

### The relative transport was rejected

`model/relative_grammar.py` implemented PBCNet2.0's Siamese *relative*
formulation without geometry, plus an AdaMBind-style leave-one-out
label-consistency credit, with no MAML, inner loop or test-time gradient. All 17
Stage 1 algebraic and synthetic gates pass.

On `meta_val` at matched seed and budget the operator is **inert**: `full` minus
`level_only` is +0.001 / -0.002 / -0.046 / -0.019 at k=1/2/3/5, and the
permutation gap is identically zero at k>=2 — the algebraic signature of flat
weights and a null operator. Neither explicit difference supervision nor removing
the reliability credit changed this. Rejected under its preregistered gates.

### Why both transports failed

Ranking on the same `meta_val` episodes, `full` against `level_only`:

| checkpoint | budget | k=2 CI | k=3 CI | k=5 CI |
|---|---:|---|---|---|
| Stage 5 grammar arm | 800 steps | 0.576 / 0.570 | 0.591 / 0.570 | 0.583 / 0.570 |
| Stage 4 grammar, three seeds | 2000 steps | worse in 2/3 | worse in 3/3 | worse in 3/3 |

At 800 steps the query-specific gate **improves** the concordance index; at 2000
steps, same split and same architecture, it **degrades** it in 9 of 12 (seed, k)
cells. The Stage 4 ranking failure is therefore a property of the mechanism
under optimization, not of the consumed `meta_test` split.

Mechanistic explanation, supported by both stages: the transport is trained on
squared error; the MSE-optimal use of k noisy support residuals is shrinkage
toward their mean; a level shift is constant across queries and cannot change
ranking; so gradient descent buys MSE with shrinkage and pays in within-target
discrimination, and buys more of that trade the longer it trains. The signed
difference operator does not escape it — its MSE optimum is `delta = 0`, the same
solution in a different parameterisation, which is why it went inert rather than
wrong.

**The blocker is the objective, not the operator.** Any query-specific channel
trained to minimise squared error on cold-target episodes converges to level
calibration whatever its functional form.

## Stage 6: chemistry-grounded support weighting (accepted at k>=2)

`meta_fewshot/stage6_bottleneck_20260815/`. Corrects the Stage 5 conclusion.

A label-and-chemistry-only audit shows the MSE optimum is **not** the support
mean: a fixed Morgan/Tanimoto kernel (production 1024-bit contract) beats it by
0.19 / 0.21 / 0.25 MSE at k=2/3/5 on `meta_val`, and the within-target
correlation between Tanimoto similarity and absolute affinity gap is -0.35. A
frozen-checkpoint audit shows the same weighting helps the *residuals* of a
`grammar` checkpoint that never saw the mechanism, so the bottleneck was the
similarity representation, not the residual decomposition and not the objective.

`model/similarity_grammar.py` adds `w_qk = softmax_k(gamma * Tanimoto)`. It is a
**fixed kernel with learned scalar calibration** (`gamma` settles at 7.98-7.99
from an initialisation of 8.0), inactive at k=0 and degenerate at k=1.

Primary evidence is **within-checkpoint**: the same trained model with and
without the weighting, three seeds, complete banks.

| split | k | dMSE | dCI | dSpearman | component LB>0 |
|---|---:|---:|---:|---:|---|
| meta_val | 2 | +0.119 | +0.038 | +0.106 | MSE/CI/rho |
| meta_val | 3 | +0.154 | +0.053 | +0.156 | MSE only |
| meta_val | 5 | +0.235 | +0.104 | +0.275 | MSE/CI/rho |
| meta_test | 2 | +0.097 | +0.047 | +0.086 | MSE/CI/rho |
| meta_test | 3 | +0.141 | +0.106 | +0.250 | MSE/CI/rho |
| meta_test | 5 | +0.118 | +0.143 | +0.326 | MSE/CI/rho |

18/18 positive point estimates and permutation gaps +0.40 to +0.51 in every
seed-k cell against the incumbent's +0.06 to +0.20. **This is the first
mechanism in this project to improve squared error and ranking together,
reproducibly, on both splits.**

Two required qualifications: component-level lower bounds are **not** all
positive (16/18 exclude zero; `meta_val` k=3 CI and Spearman cross zero, so that
ranking cell is unestablished), and every interval is **conditional on the three
trained seeds** because the analysis averages seeds per (component, target)
before resampling components, so seed variance is not resampled.

Not established, and explicitly not claimed:

* **superiority over the incumbent `grammar` transport.** Cross-arm results
  contradict across splits — F wins on `meta_val`, the incumbent wins at every
  k>=1 on `meta_test` — and no lower bound excludes zero on `meta_test`;
* any k=0 or k=1 effect: the mechanism is inactive at k=0 and degenerate at k=1,
  with `full - level_only` exactly 0.000000 in all seeds;
* protein-specific transport: Tanimoto consumes no protein information, so the
  wrong-protein control is a full-system perturbation here, not mechanism
  evidence.

Audit corrections applied: unused `key`/`log_temperature` frozen when
`use_learned_key=False` (they had `grad=None`; AdamW skips such parameters, so no
completed run was invalidated), zero-fingerprint and "zero-learning" wording
corrected, signal audit rerun at the production fingerprint width, and every
decision re-derived on the complete 44-episode `meta_val` bank rather than the
6-episode automatic bank. Details in `AUDIT_ADDENDUM.md`.

## Stages 9-10: k=0 diagnosed, and the first accepted k=0 improvement

`meta_fewshot/stage9_k0_decomposition_20260815/`,
`meta_fewshot/stage10_retrieval_prior_20260815/`.

### The k=0 diagnosis (training-free, meta_train-only indices)

`MSE = calibration + shape` on `meta_val` k=0, 50 episodes:

| estimator | MSE | calibration | shape | CI | Spearman |
|---|---:|---:|---:|---:|---:|
| `global_mean` | 2.420 | 1.677 | 0.743 | 0.500 | — |
| `protein_neighbor_esm` | 2.079 | 1.336 | 0.743 | 0.500 | — |
| accepted `model_f0` | 1.821 | 1.081 | 0.740 | 0.525 | 0.075 |
| `ligand_neighbor_b24` | 1.625 | **0.697** | 0.928 | 0.635 | 0.313 |
| composed retrieval *(transductive)* | 1.356 | 0.697 | 0.659 | 0.660 | 0.392 |
| `target_oracle` *(oracle)* | 0.743 | 0 | 0.743 | 0.500 | — |

The composed row re-centres on the query panel, so it is a **transductive
diagnostic upper bound** (grade: exploratory), not a deployable predictor. The
best per-query train-only estimator is `ligand_neighbor_b24` at 1.625, a 10.8%
reduction against the model.

**k=0 error is 59% target-level calibration.** Re-centring the model on the true
target mean gives 0.7403 against a flat constant's 0.7430: **the trained
zero-shot endpoint has essentially no within-target ligand discrimination**
(CI 0.525 against a 0.500 coin flip). Protein retrieval is weak; ligand
retrieval beats the protein-conditioned model at calibration using no protein.

### The accepted intervention

A training-free blend of a `meta_train`-only retrieval prior into the endpoint,
`f0' = 0.5*f0 + 0.5*retrieval`, transport unchanged:

| k | accepted | with prior | reduction |
|---|---:|---:|---:|
| 0 | 1.612 | **1.415** | **12.3%** |
| 2 | 1.064 | 0.884 | 16.9% |
| 5 | 0.896 | 0.687 | 23.3% |

**Grade: WITHDRAWN by Stage R0.** As originally measured: component bootstrap
+0.198 [+0.012, +0.416], positive in 3/3 seeds, CI +0.026, Spearman +0.081. All
five Stage R0 gates then failed on the identical population. The corrected
statement is that the prior helps query ligands already present in `meta_train`
(1.3581 -> 1.0114) and harms those that are not (2.8019 -> 3.0193), and that
tuning it on the population it is reported on is worth 0.468 MSE by itself. It
is retained as a named baseline only.

**Mandatory caveats, all verified by recomputation.**

* The CD-HIT40 split is component-hard on proteins only. **305 of the 624 query
  cells (48.9%) contain a ligand that appears verbatim in `meta_train`.** On the
  12 targets with no exact overlap the reduction is **+0.050 [-0.074, +0.175]**,
  about a quarter of the headline effect and statistically unresolved. The 12.3%
  therefore describes a population that is half exact recall.
* **6 of 10 components improve**; four regress by -0.155 to -0.012. The bank has
  10 components, not 11.
* `beta=24`, the retrieval source and `w=0.5` were selected on `meta_val` — the
  same population the bootstrap resamples. Selection and inference are not
  separated.
* `w=0` matches the checkpoint only approximately: the script fixes `beta=8`
  while the checkpoints learned 7.9743/7.9849/7.9897.
* The prior is an **offline evaluator**. It is not in `ARCHITECTURES`, not in any
  checkpoint, and not on the standard evaluation path.

## Stage R cycle: double-cold protocol, factorization and routed training

`meta_fewshot/stageR0_retrieval_falsification_20260815/`,
`stageR1_double_cold_split_20260815/`,
`stageR2_representation_discriminator_20260815/`,
`stageR3R4_level_shape_20260815/`.

### The protocol changed, and it invalidates earlier development decisions

`dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1` is a governed
two-axis split (CD-HIT40 protein components x Bemis-Murcko scaffold clusters,
then document closure, then eligibility). Verified overlap with the training
block is **zero** on exact ligand, scaffold, protein component and document, and
**81.6%** of development ligands are below Tanimoto 0.40 to every training
ligand. Development: 41 targets / 19 components. Confirmation: 22 targets /
10 components, **never opened**.

The cost is real: 44.1% of the corpus is retained and the training block falls
from 12,633 to 5,643 cells. Absolute numbers on this protocol are therefore not
comparable with anything above; only within-protocol contrasts are.

### Two core innovations, one admitted and one rejected

**Innovation A — level-shape factorized predictor** (`model/level_shape.py`):
`f = ligand_prior(L) + target_level(P) + centered_interaction(P,L)`, where the
centered branch subtracts a per-protein constant built from **learned anchor
ligand embeddings**. It therefore has exactly zero mean in the anchor basis and
cannot carry a target-level offset, while every prediction stays inductive.
19 structural gates pass.

**Rejected on performance.** Under ordinary training it is significantly *worse*
than the incumbent trunk at k=0 (**-0.268 [-0.482, -0.059]**), k=2 and k=5. With
the full training method it recovers to +0.094 [-0.142, +0.314] against the
incumbent — a 4.3% point gain that does not clear a component bootstrap.

**Innovation B — counterfactual level-shape gradient-routed training**
(`scripts/train_level_shape.py`). The squared error decomposes exactly into
`mean(p-y)^2 + var(p-y)`; the method routes each term only to the components
that own it, keeping one joint scalar prediction and one backward pass. Three
counterfactual contrasts run in the same step, each routed so it cannot be
satisfied by the wrong module.

**Admitted as a real and major performance source, relative to conventional
training of the identical architecture:**

| contrast (same architecture) | k=0 MSE | k=0, `< 0.40` tier |
|---|---|---|
| **full method vs ordinary training** | **+0.3618 [+0.0275, +0.7002]** | **+0.4174 [+0.1248, +0.7308]** |
| routing alone | +0.0944 [-0.1556, +0.3597] | +0.1203 [-0.1075, +0.3703] |
| counterfactual alone | +0.2674 [-0.0695, +0.6535] | +0.2971 [-0.0077, +0.6698] |

15.0% at k=0 and 17.5% on the low-similarity tier, three seeds, component lower
bounds above zero. **Neither half is separately resolved** — only the whole
method clears the bootstrap, so no claim is made for routing or counterfactual
supervision on its own.

The premise is measured, not assumed: under ordinary training the level and
shape objectives have gradient cosine **-0.334** on the interaction trunk and
**-0.532** on the ligand encoder. Under routing the interaction's level-gradient
is exactly **0.000e+00**.

### First resolved protein specificity, and what falsified the hypothesis

Against a **similarity-matched donor from the same evaluation split**, the
routed arm's wrong-protein gap is +0.4216 [+0.1180, +0.7436] at k=2 — **4.2x**
the incumbent's +0.1016 [-0.0060, +0.2234], which crosses zero. At k=0 every
arm's gap crosses zero, so the claim is limited to k >= 2.

But the cycle's central hypothesis is **falsified**. Removing the gradient
conflict improved calibration (1.4897 -> 1.1309) and left within-target shape
flat (0.9274 -> 0.9244, against the incumbent's 0.9130). Deleting the
interaction branch entirely changes k=0 MSE by **-0.0241 [-0.0587, +0.0024]**.
The protein specificity lives in `target_level`, not in `centered_interaction`.
**The level/shape conflict was real and fixable, and it was not the binding
constraint on within-target ordering.**

### Admission decision

**No arm meets the zero-shot admission target** (Z1 4.3% against 10%; Z2 crosses
zero; **Z3 fails with a resolved CI regression of 0.049 [0.010, 0.090]**; Z4 2 of
3 seeds; Z5 fails at k=0; Z6 unresolved). Under the preregistered rule
`meta_test` is therefore **not opened**, and the double-cold confirmation
population remains pristine for future work.

## Stage R5-R9 cycle: contract repairs, the relative-transport closure, and the pair-level diagnosis

Full cycle narrative: `meta_fewshot/README.md` (unified R5-R8 summary);
artifact inventory and per-stage verdicts: `EVIDENCE_LEDGER.md`; gates and
reports: `meta_fewshot/stageR{5..9}_*/`.

### Standing conclusions of the cycle (all measured, all retained)

1. **The experimental contract is repaired** (R5): same-split wrong-protein
   donors with meta_train-only whitening, aggregated gradient-cosine
   diagnostics, a physical `meta_test` seal (`QPSMPData
   include_meta_test=False`), and complete artifact records (config, split
   hash, seed, checkpoint sha256, per-target predictions, donors, activation
   statistics, gradient coverage). Covered by
   `tests/test_stage0_contract_fixes.py`.
2. **The relative-transport/gate model family is closed.** Three transport
   mechanisms were each falsified under their preregistered gates (R6a
   saturating gate, R6b additive correction, R7 linear rho gate), all with
   the same signature — deployment-inert (`nogate` gap ~0.000) while their
   training gradients disturb calibration. It is the seventh query-specific
   channel in the project with that signature, now under ranking-primary
   objectives, so the objective is no longer the explanation.
3. **The shape-first training method is the project's first real
   within-target shape source** (measured, R7/R8): same-architecture shape
   term 0.943 -> 0.896, k=0 cliff sign 0.512 -> 0.598, k=5 cliff sign 0.675
   -> 0.768 (best on record). It did not convert into a global CI or k=0
   MSE gain: R8's best arm ties A0 at k=0 (2.167 vs 2.149) while CI
   regresses (0.535 vs 0.580), so the family was closed for the double-cold
   zero-shot target under its preregistered rule. `meta_test` was never
   opened.
4. **The pair-level audit localizes the CI regression** (R9,
   `stageR9_cliffweight_20260816/PAIR_AUDIT_meta_val.json`): at k=0 the only
   component-resolved stratum is the **mid-similarity band
   (0.4 <= Tanimoto < 0.6): +0.119 [+0.022, +0.220]** per-target sign
   accuracy against A0 — the band immediately below the activity-cliff
   weight's discontinuity at 0.6. Cliff pairs themselves improve
   (-0.049, unresolved). Low-similarity pairs are unresolved at the target
   level (-0.022); mid-gap pairs (0.5-1.0 pK) are near-resolved
   (+0.120 [-0.008, +0.263]). The R9 cliff-weight dose response (1x/2x/4x)
   tests the discontinuity hypothesis directly.
5. **The routed level readout converges to incumbent calibration at the
   full budget** (attention pooling; 1.271 vs A0's 1.236), and the
   counterfactual/identifiability machinery is verified at the gradient
   level — transferable components for the next design.

### Current experiments (R9-R11, preregistered and executed)

- R9 (`stageR9_cliffweight_20260816/`): pair-level audit + cliff-weight dose
  response. **The x4 activity-cliff weight is a net negative for ranking
  itself** — C1 (w=1) beats B1 (w=4) on global CI (0.562 vs 0.535) *and* on
  cliff pairs (0.606 vs 0.577); C2 (w=2) gives the first 3-seed k=0 below A0
  (2.119, unresolved) with calib 1.218. The only component-resolved CI-loss
  stratum was the mid-similarity band (0.4-0.6), which the weight removal
  closed.
- R10 (`stageR10_variance_20260816/`): `shape_variance_weight 1.5 -> 0.5` on
  the C1 base — **falsified** (CI 0.552 vs 0.562; cliff gain lost). The
  variance term is not the margin-compression cause.
- R11 (`stageR11_grammar_shape_20260816/`): shape-first routing applied to
  the incumbent trunk itself, zero architecture change — **falsified**: the
  incumbent's calibration lives in the interaction branch; routing the level
  away from it degrades calibration 1.236 -> 1.488 and CI 0.580 -> 0.525.
  The level/shape routing trades calibration for shape on **every**
  architecture tested.
- R12 (`stageR12_margin_20260816/`): margin-ranking (hinge) shape objective
  on the C2 base — **falsified as the actionable lever**: CI +0.003 only
  (0.548 -> 0.551). The margin compression is a symptom of the shape
  branch's expressivity, not the loss form.
- R12 (`stageR12_margin_20260816/`): `REPORT.md` and `RESULT.json` were
  backfilled on 2026-08-16 from the retained comparison artifact — the stage
  had run without them. Gate **M5 is recorded as not evaluable**: the
  artifact carries `D2_vs_A0` but no `D2_vs_C2` contrast, so the
  preregistered bootstrap against its stated control was never computed.
- R13 (`stageR13_shape_direct_20260816/`): direct interaction-head shape
  with difference supervision (supervision-leak fix) — **gate-blocked at
  Stage 1**: the MLP shape branch collapses under the shape variance term
  on the synthetic interaction task (mean CI 0.60, gap 0.14 against gates
  0.70/0.20). The suite collects **18 gates: 16 pass, 2 recorded `xfail`**.
  An earlier revision reported "16 gates / 15 of 16 pass", which is
  inconsistent with itself and with the suite; corrected 2026-08-16.

### Standing conclusion of the R0-R13 chain

The shape-first training method produces the project's first real
within-target shape and activity-cliff gains, and **no tested configuration
has converted them into a joint k=0 MSE + CI improvement over the
incumbent.**

The k=0 Pareto frontier over the two preregistered primary metrics (MSE
down, CI up) is exactly three whole configurations, on double-cold
`meta_val`, three seeds:

| arm | k=0 MSE | k=0 CI |
|---|---:|---:|
| B3 (R3R4 full method) | **2.055** | 0.531 |
| C2 (R9, cliff w=2) | 2.119 | 0.548 |
| A0 (incumbent) | 2.149 | **0.580** |

No model achieves both ends; none of the MSE differences against A0 is
resolved by a component-level paired bootstrap. Earlier phrasing of this as
"MSE frontier 2.055-2.119, CI frontier 0.548-0.580" mixed metrics from
different models and described a configuration that does not exist.

The best activity-cliff ordering, **k=5 cliff sign 0.782**, is a
double-cold **`meta_val` development record** on arm C1 (R9) — and C1 is
Pareto-dominated on both primary metrics (MSE 2.235, CI 0.562). It has never
been measured on the sealed `meta_test` and is not a generalization claim.
What the R9 dose response does establish is mechanism, not leaderboard: the
cliff-ordering ability comes from the shape-first training itself, not from
the activity-cliff pair weight (C1 at weight 1 beats B1 at weight 4 on
global CI *and* on cliff pairs).

The consolidated reachable-boundary statement is
`report/BOUNDARY_20260816.md`. The double-cold `meta_test` remains sealed
and unopened.

## Controlled model-stage decisions

| Stage/family | Main observation | Decision |
|---|---|---|
| B/C slot-gated readout | Three seeds: no stable MSE gain; binding gaps worsened | Reverted |
| D atom-aware localization | k0 and wrong-protein specificity weakened | Reverted |
| F/G primitive matching | k0-k2 regressed; learned primitive gate ~0.0012 | Reverted |
| H/I residual context | k3/k5 improved in two seeds, k0 regressed strongly | Reverted |
| J/K detached context | k0 remained worse | Reverted |
| L 0.05 pK label noise | Short run worsened k0/k1/k2 and wrong-protein control | Reverted |
| BPSF v2 relevance/shared latent | Support state remained near zero and lost to level | Rejected |
| HyperSAR/PBC-inspired matching | Permutation gap near zero | Rejected |
| D-MEMT | Permutation/foreign controls failed | Rejected |
| CIPF+TERM | Most gain was scalar level | Rejected |
| L-CIPF+ELMT | Controlled development failed its admission boundary | Rejected |
| E no-warmup baseline | Retained comparator | Retained |
| Stage 0 budget probe (4x steps, retained arch) | Validation improved, every test metric worsened | Rejected |
| Stage 3 arm A (retained arch, 2000 steps + cosine) | No improvement; endpoint collapsed further | Rejected |
| Stage 3 capacity 1.8M -> 7.3M | No consistent MSE gain; fixed control sign only | Weak, retained for control validity |
| **Interaction grammar trunk** | **12-18% MSE reduction at every k, three seeds, both banks** | **New development baseline; not admitted** |
| **Transferability gate `rho`** | No positive bootstrap lower bound; ranking degrades at every k, and the degradation grows with training budget | **Claim rejected; retained in code as the tested mechanism** |
| Cartesian / equivariant interaction encoder | 0 of 17,717 deployment cells have a common-frame complex | Rejected on data; `model/cartesian.py` verified and unused |
| **Signed relative-difference transport** | Algebraically sound (17 gates) but inert on real data: `full == level_only`, permutation gap identically zero at k>=2 | **Rejected; opt-in `--arch relative`, retained as evidence** |
| Explicit per-(query, support) difference supervision | Worse than the unsupervised relative arm at k=0,1,2,3 | Rejected |
| Leave-one-out label-consistency credit | No effect while the operator it weights is null | Not evaluable; rejected with its parent |

## Failure localization (double-cold era)

Confirmed limitations:

1. No trained model orders within a cold target better than the incumbent
   at k=0 (best-ever CI 0.575 at 300 steps, A3; incumbent 0.580 at 1200
   steps) — and every arm that improves the shape term pays for it in
   global CI, concentrated in the mid-similarity band (R9 audit).
2. Seven query-specific few-shot channels (across three model families)
   measured deployment-inert while disturbing calibration under their
   respective objectives. The retained useful few-shot mechanism remains
   the Stage 6/7 fixed Morgan/Tanimoto residual weighting.
3. The double-cold k=0 MSE target (-10% vs the incumbent's 2.149) has never
   been approached: best recorded is B3's 2.0554 (-4.3%), best of the R5-R8
   cycle is R8 B1's 2.167 (-0.8%, unresolved tie).
4. `meta_test` (22 targets / 10 components) is pristine and physically
   sealed; it may be opened once, only after every preregistered meta_val
   gate passes.
5. The active bank has no complex coordinates; no atomic 3D claim is
   possible, and the Cartesian family stays closed on data.

Open hypotheses (each falsifiable, preregistered):

- the R9 discontinuity hypothesis: the activity-cliff pair weight's
  threshold at Tanimoto 0.6 starves the 0.4-0.6 band (audit-resolved
  regression; dose response in R9);
- the shape-variance hypothesis: the variance term compresses shape margins
  on uncertain pairs, hurting their ordering (B1 margins 0.090 vs A0 0.121);
- the trunk hypothesis: shape-first training applied to the incumbent
  grammar trunk itself (best calibration 1.236) rather than to factorized
  trunks (`scripts/train_grammar_shape.py`, preregistered as R10/R11).

## Engineering evidence (double-cold, 2026-08-16)

| quantity | A0 incumbent (R3R4) | reltransport (R7/R8) |
|---|---:|---:|
| trainable parameters | 7,294,171 | 1,885,113 |
| steps per run | 1200 | 1200 |
| peak CUDA memory | ~6,500 MB | 550-650 MB |
| wall time per 1200-step run | ~150 s | ~300-330 s |
| checkpoint selection | val admission score | mean-over-k MSE |
| zero-gradient tensors | 0 | 0 (gate suite; k=0/k=1 semantic exceptions documented) |

Complete maintained test suite (82 pytest modules), 2026-08-16. The pre-audit
suite was 410 passed / 3 skipped / 2 xfailed in 410 s. After adding
`tests/test_research_record.py` and splitting the regression-suite tiers:

| tier | command | outcome | wall |
|---|---|---|---:|
| default | `pytest tests -q` | 413 passed, 9 skipped | 105 s |
| research gates | `RUN_RESEARCH_GATES=1 pytest tests -q` | 416 passed, 3 skipped, 2 xfailed | 410 s |

The six deferred tests are the synthetic **training** gates of the two closed
families; their verdicts are retained as immutable evidence. See
`docs/PROJECT_FILE_ORGANIZATION.md` "Regression-suite tiers".
`scripts/evaluate_qpsmp.py` releases CUDA memory between model seeds.

## Admission requirements for the next model

Unchanged, plus two additions forced by this stage:

1. no material k=0 regression and preferably a zero-shot improvement;
2. genuine query-specific k=1 improvement beyond scalar level calibration;
3. consistent k=2/3/5 gains;
4. full beating level/adaptation-cut with a positive component-bootstrap lower
   bound on the preregistered primary metric;
5. correct support beating permuted and magnitude-matched wrong support —
   **note that permutation leaves `mean(r)` invariant, so this is the sharpest
   available test of the query-specific channel**;
6. CI, Spearman and pairwise sign/ranking improving with MSE — **an MSE gain
   bought by shrinking toward the target level will fail here**;
7. reproducible checkpoint evaluation and complete gradient/activation logs;
8. measured parameter count, wall time, utilization and safe peak GPU memory.
