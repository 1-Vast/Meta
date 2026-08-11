# Meta-learning Transfer & Robustness Audit

**Status:** static research audit completed 2026-08-11; no model was trained and
no confirmation label was opened. The frozen operator, positive-ridge section
and theory were not modified.

## 1. Current problem localization

The question is not whether another optimizer can lower ordinary DTA error. It
is whether source-task selection and realistic uncertainty handling can make a
`k=1/2/3/5` support-conditioned section transfer more stably to unseen targets
without preserving the same ligand, intercept or wrong-protein shortcut.

The current evidence imposes three hard boundaries.

1. The uncentered positive analytic ridge is retained. The targeted V1 repair
   reduced RMSE by `7.6-8.9%`, but absolute quality, support specificity and
   all-k partner-identity Gates still failed.
2. The existing ATS and fixed support-noise arms already failed their matched
   controls. ATS did not beat its equal-capacity permuted-statistics null at
   every k, and support noise did not beat clean training. These are negative
   results, not untested ideas.
3. F-132 is stronger negative biological evidence: on 219 fresh components the
   exact-pair residual was worse than the capacity-matched additive atom/residue
   model (`RPS 0.040522` versus `0.039302`). Therefore R1 affinity-field
   training and V1 exact-pair integration are closed. A scheduler cannot turn
   this into an admitted partner-specific observable.

The operative failure class is consequently **information/identifiability plus
support instability**, not merely under-training. An optimization improvement
that survives protein/support corruption must be recorded as
`TRAINING_REGULARIZATION_USEFUL_BUT_NOT_META_IDENTIFIABILITY_REPAIR`.

## 2. AdaMBind mechanism audit

Primary sources are the [Nature Communications article](https://www.nature.com/articles/s41467-026-70554-5),
its complete 19-page Supplementary Information, and official repository commit
[`01a169a6`](https://github.com/Moohyun-w/AdaMBind/tree/01a169a6d62fba0d6c003f47bfba539e55f5b344).
The supplement was checked as OOXML (84 paragraphs, 22 tables), exported to PDF
and visually inspected on all 19 pages.

### 2.1 Paper-level computation

Each protein is a task. Drugs measured on that protein are split into support
and query examples. A two-layer drug GNN, three-layer protein CNN and fusion
regressor are adapted with MAML. The adaptive module receives, for each
candidate task:

- query loss after support adaptation;
- support/query gradient agreement at the shared initialization; and
- training progress.

The supplement describes two bidirectional LSTMs over the loss and layer-wise
gradient features, a set-context mapping and a scalar sampling score. It then
describes an approximate bilevel cycle: adapt candidate tasks, score them,
sample a meta-batch, take a provisional meta update, measure separate
meta-validation tasks, update the scheduler, rescore/resample, and update the
meta-learner. This is derived from [Yao et al. ATS](https://proceedings.neurips.cc/paper/2021/hash/3dc4876f3f08201c7c76cb71fa1da439-Abstract.html).

High weight has only a conditional meaning: the task has low current query
loss and/or support/query gradients compatible with the current model, in a
mixture that happened to improve validation loss. It does **not** identify
label reliability, biological information or transfer to a remote protein.
The paper's “easy-to-hard” account is a post-hoc trajectory interpretation;
neither the score definition nor the code imposes monotone easy-to-hard order.

### 2.2 Label perturbation

The paper and supplement use zero-mean uniform perturbations with a noise-scale
grid from `0.1` to `0.6`. The described training intervention touches support
and query labels; clean testing is separate. The reported hyperparameter sweep
is not an empirical assay-noise model: it contains no replicate-derived
variance, endpoint-conditioned scale, clean/noisy robustness slope, or
cross-assay calibration test. The ablation tables report five repeats of
aggregate DTA metrics, chiefly in the majority setting, and do not establish
that noise specifically stabilizes the target-conditioned section.

### 2.3 Task and data construction

Supplementary Table 1 reports:

| Dataset | Proteins | Drugs | Pairs |
| --- | ---: | ---: | ---: |
| BindingDB | 1,088 | 9,862 | 42,203 |
| KIBA | 229 | 2,068 | 118,254 |
| Davis | 379 | 68 | 30,056 |

Random tasks use an 8:1:1 split. The novel-task protocol clusters protein
sequences with CD-HIT at 40% identity and then allocates clusters 8:1:1.
Few-shot uses five supports and “majority” uses forty; reported comparisons use
five repetitions. This target-as-task organization and explicit support/query
adaptation are transferable. The public CSV snapshot, however, exposes only
SMILES, target sequence and scalar affinity. It does not carry assay identity,
document provenance, replicate groups or the published fixed split manifests.
KIBA is also a composite score rather than an assay-matched Ki/Kd label. These
files cannot support MetaSieve's assay-aware uncertainty or dependency closure.

### 2.4 Official-code audit and paper/code boundary

The official code is useful evidence about reproducibility, but is not safe to
reuse directly.

- [`train.py`](https://github.com/Moohyun-w/AdaMBind/blob/01a169a6d62fba0d6c003f47bfba539e55f5b344/train.py)
  concatenates `train_idxs + val_idxs`, then samples both scheduler-training and
  scheduler-validation tasks from that same pool. The claimed separate
  validation-task bilevel role is therefore not preserved, and task-role
  overlap is possible.
- The scheduler update is a REINFORCE loss on sampled indices scaled by
  validation loss. This is related to ATS, but it is not the differentiable
  one-step update depicted in the AdaMBind supplementary algorithm.
- [`scheduler.py`](https://github.com/Moohyun-w/AdaMBind/blob/01a169a6d62fba0d6c003f47bfba539e55f5b344/model/scheduler.py)
  uses detached layer-wise cosine similarities and samples with replacement by
  default; duplicate tasks can occupy one meta-batch.
- `adaptive_tasks` is parsed but never gates the training path. In
  [`Trainer.py`](https://github.com/Moohyun-w/AdaMBind/blob/01a169a6d62fba0d6c003f47bfba539e55f5b344/model/Trainer.py),
  `noise=0` leaves local `y` undefined in training, while query noise is drawn
  but the reported query MSE is computed from clean labels. Thus the public
  flags do not provide a reliable scheduler/no-noise ablation contract.
- [`DataSplit.py`](https://github.com/Moohyun-w/AdaMBind/blob/01a169a6d62fba0d6c003f47bfba539e55f5b344/utils/DataSplit.py)
  shuffles each target and takes a prefix as support. No scaffold, document,
  assay or recipient closure is enforced there.

The article's benchmark and screening results can support “promising
engineering method,” but not the stronger MetaSieve claim that the selected
tasks or noise expose partner-specific affinity information. The official
implementation receives **`REFERENCE_ONLY`**; direct code reuse receives
**`REJECT`**.

## 3. AdaMBind versus MetaSieve

| Axis | AdaMBind | Current MetaSieve | Decision |
| --- | --- | --- | --- |
| Task | one target and its ligands | same | `ADOPT` |
| Adaptation | full-model MAML fine-tuning | `d<=5` analytic positive ridge | keep MetaSieve; MAML `REJECT` |
| Scheduler inputs | query loss, support/query gradient similarity, progress | log query loss, cosine agreement, progress | already implemented in smaller form |
| Scheduler supervision | paper: validation-task bilevel; code: overlapping train+val pool and REINFORCE | disjoint `meta_val` target clusters, first-order gradient utility | MetaSieve isolation is stricter |
| Sampling | code defaults to replacement | without replacement | keep MetaSieve |
| Noise | generic uniform, paper describes support+query | source support only, clean query | current fixed noise already failed |
| Novel target | CD-HIT40 cluster split | complete CD-HIT40 clusters plus controls; full protocol also requires scaffold/document/recipient closure | keep stricter protocol |
| Biology controls | aggregate DTA metrics/ablations | ligand-only, wrong protein/support, permuted labels, intercept controls | keep MetaSieve |
| Confirmation | benchmark and application evaluations | fresh dependency-closed confirmation required | AdaMBind data are consumed engineering data only |

## 4. Current implementation audit (Gate 0)

### 4.1 What the scheduler optimizes

`model/metasieve_v1.py:152-183` defines a three-input MLP. During training,
`research/meta_fewshot/train_main_v1.py:328-372` supplies:

1. `log1p(query episode MSE)`: current post-support predictive difficulty,
   which mixes true difficulty, label noise, representation failure and model
   misspecification;
2. cosine between gradients of support population residual and clean query
   episode loss: local directional compatibility, not the unnormalized inner
   product in the ATS derivation and not proof of biological transfer; and
3. scalar training progress, identical for every candidate in an iteration.

The scheduler is fitted to approximate a soft distribution over candidate
query-gradient cosine with one clean, disjoint meta-validation task. Selected
source tasks then update the model. Inputs and the validation gradient are
detached, so there is no second-order gradient through the provisional model
update. This is a source-only first-order task-ranking surrogate, not full
AdaMBind/ATS bilevel optimization. The equal-capacity `ats_null` arm permutes
the three-statistic rows, and sampling is without replacement. The scheduler
is not loaded at test time.

Therefore, high scheduler weight means **predicted current-gradient utility to
one sampled meta-validation episode**. It may favor an easy/familiar target, a
low-noise assay, a dense ligand series or the same shortcut present in both
support and query. The current implementation never conditions directly on
protein homology, scaffold overlap, sample count, assay type, affinity spread,
replicate reliability, target family or measured interaction strength.

### 4.2 Curriculum-shortcut and corruption evidence

No logged scheduler-weight audit currently separates those covariates. More
importantly, V1's outcome already supplies a destructive warning: ATS failed
its all-k matched-null Gate while correct-versus-wrong protein and support
specificity also failed. Gradient agreement therefore has no current license
to mean biological transferability. A new scheduler may proceed only after an
offline ranking audit shows that its score changes under correct-to-wrong
protein/support corruption and is not explained by homology, scaffold
familiarity, task size, assay or baseline loss.

### 4.3 Current label noise

`model/metasieve_v1.py:186-195` correctly converts a declared standard
deviation to `Uniform(-sqrt(3)sigma,+sqrt(3)sigma)`. The V1 runner converts a
fixed pK half-width of `0.2` into standardized units and perturbs only source
support labels with a separate deterministic CUDA generator. Query labels,
sealed evaluation labels and test-time model state remain clean; clean and
noisy-support arms are paired.

This implementation is reproducible, but the distribution is arbitrary with
respect to assay and provenance. At `k<=5`, support noise directly changes the
ridge right-hand side and hence the identifiable section; it is not analogous
to harmless large-batch augmentation. The formal run found no all-k advantage
over clean training. Fixed uniform support noise is therefore **`REJECT` as a
candidate repair** and retained only as a negative baseline.

## 5. Cross-domain literature matrix

The ten required questions are abbreviated as: **Q1** original root cause;
**Q2** same MetaSieve root; **Q3** mathematical/statistical isomorphism; **Q4**
source of new information; **Q5** optimization or identifiability; **Q6** clean
meta-validation required; **Q7** leakage risk; **Q8** suitability for `k<=5`;
**Q9** classification dependence; **Q10** cheapest falsifier.

| Work | Q1-Q5 | Q6-Q10 | Verdict |
| --- | --- | --- | --- |
| [Yao et al., ATS, NeurIPS 2021](https://proceedings.neurips.cc/paper/2021/hash/3dc4876f3f08201c7c76cb71fa1da439-Abstract.html) | Q1 noisy/imbalanced tasks under budget. Q2 partly: source tasks differ. Q3 only at task-selection level; gradient agreement is not biological information. Q4 clean validation-task loss. Q5 optimization/generalization, not pair identifiability. | Q6 yes. Q7 high unless validation tasks are dependency-isolated. Q8 computationally yes, statistically unstable. Q9 experiments are classification plus a drug-activity regression benchmark; mechanism is loss-agnostic. Q10 compare ranking before/after protein/support corruption and against a permuted-statistics null. | `ADOPT_WITH_MODIFICATION` |
| [Ren et al., ICML 2018](https://proceedings.mlr.press/v80/ren18a.html) | Q1 biased/noisy examples. Q2 partly: heterogeneous assay cells. Q3 gradient reweighting isomorphic only if a trusted validation estimand exists. Q4 small clean unbiased validation set. Q5 optimization/robustness; no new interaction information. | Q6 yes. Q7 severe if validation documents/scaffolds leak. Q8 weights may have extreme variance at very small k; use source cells, not test supports. Q9 no essential classification assumption. Q10 one-step nonnegative weights versus uniform on a sealed meta-val component, with partner corruptions. | `ADOPT_WITH_MODIFICATION` |
| [Shu et al., Meta-Weight-Net, NeurIPS 2019](https://proceedings.neurips.cc/paper/2019/hash/e58cc5ca94270acaceed13bc82dfedf7-Abstract.html) | Q1 unknown loss-to-weight rule under noisy/biased data. Q2 possibly at assay-cell level. Q3 only optimization-isomorphic; loss alone cannot distinguish hard biology from noise. Q4 clean meta-data learns the mapping. Q5 optimization, not identifiability. | Q6 yes. Q7 same meta-data reuse risk. Q8 poor direct fit for five support labels; plausible only across many source cells. Q9 not inherently, although evidence is mostly classification. Q10 test whether learned weights are predicted by assay reliability after controlling for loss and whether benefit survives clean data. | `REFERENCE_ONLY` |
| [Chen et al., Eigen-Reptile, ICML 2022](https://proceedings.mlr.press/v162/chen22aa.html) | Q1 sampling and label noise destabilize gradient-based meta-learning. Q2 support sampling instability exists. Q3 not algorithmically isomorphic: MetaSieve adapts by closed-form ridge, not parameter trajectories. Q4 historical inner-loop parameter directions and ISPL-selected samples. Q5 optimization stability. | Q6 ISPL requires a reliable selection signal/training partition. Q7 cross-task history can leak if components overlap. Q8 motivation fits few-shot, implementation does not fit the analytic section. Q9 empirical method is classification-oriented but trajectory idea is generic. Q10 compare ridge-section coefficient dispersion over resampled supports; no Reptile implementation needed. | `REFERENCE_ONLY` |
| [Bertinetto et al., R2-D2, ICLR 2019](https://openreview.net/forum?id=HyxnZh0ct7) | Q1 expensive/unstable neural fine-tuning. Q2 yes: few-shot adaptation should be small and identifiable. Q3 strong at representation-to-differentiable-ridge level. Q4 source query losses shape the representation; supports solve the task section. Q5 improves representation/conditioning, but identifiability still comes from support geometry. | Q6 ordinary isolated meta-validation for selection. Q7 standard episodic leakage risk. Q8 explicitly designed for few-shot; dual/Woodbury solve is favorable. Q9 its reported loss is classification, but ridge differentiation transfers to regression. Q10 freeze current features and measure coefficient stability, rank, conditioning and query coverage under support resampling. | `ADOPT` for the existing solver principle; no new module |
| [Oreshkin et al., TADAM, NeurIPS 2018](https://proceedings.neurips.cc/paper/2018/hash/66808e327dc79d135ba18e051673d906-Abstract.html) | Q1 a fixed metric is suboptimal across tasks. Q2 only if a protein-conditioned coordinate is already informative. Q3 weak: class-centroid metric conditioning is not affinity regression or pair geometry. Q4 support-set summary. Q5 representation/optimization; can create a support shortcut. | Q6 yes for model selection. Q7 support/query and class leakage risks. Q8 conditioning from five points is high variance. Q9 strongly classification/metric based. Q10 compare conditioning to intercept-only and wrong-support controls with the pair feature frozen. | `REFERENCE_ONLY` |
| [Agarwal et al., NeurIPS 2021](https://proceedings.neurips.cc/paper/2021/hash/ab73f542b6d60c4de151800b8abc0a6c-Abstract.html) | Q1 extreme dependence on which natural supports are selected. Q2 yes. Q3 strong diagnostic isomorphism: low-k support geometry controls the adapted solution. Q4 none; it exposes instability through support/margin analysis. Q5 diagnostic of identifiability/robustness. | Q6 no clean meta labels beyond ordinary evaluation. Q7 low if resampling remains within sealed episodes. Q8 directly applicable. Q9 margin theory is classification-specific; sensitivity audit is not. Q10 enumerate/resample support subsets and report worst/median coefficient and prediction spread. | `ADOPT` as a diagnostic |
| [Bronskill et al., TaskNorm, ICML 2020](https://proceedings.mlr.press/v119/bronskill20a.html) | Q1 batch normalization mismatches hierarchical tasks. Q2 no: current V1 section has no batch-normalization state. Q3 absent. Q4 per-task support statistics. Q5 normalization only. | Q6 ordinary. Q7 query-statistic use would leak. Q8 proposed for few-shot, but irrelevant here. Q9 tied to deep image feature normalization, not strictly classification. Q10 confirm no episode-coupled normalization exists in the retained model. | `REJECT` |
| [Minot & Reddy, Cell Systems 2024](https://doi.org/10.1016/j.cels.2023.12.003) | Q1 noisy and under-labelled biological sequence-function data. Q2 partly: affinity sources are heterogeneous and sparse. Q3 trusted-data bilevel reweighting is analogous, but antibody display labels and affinity assays are not exchangeable. Q4 experimentally trusted labels and positive/unlabelled task construction. Q5 robustness; it cannot create missing partner interaction information. | Q6 yes, genuinely trusted biological data. Q7 high if replicate/doc/scaffold families overlap. Q8 conceptually, although source experiments provide more repeated observations. Q9 no essential classification assumption, but their tasks differ from continuous affinity. Q10 use replicate-derived uncertainty on source-only assay groups and test section stability against fixed-noise and clean arms. | `ADOPT_WITH_MODIFICATION` |

## 6. Transferable mechanisms

1. **Representation/solver separation.** R2-D2 independently supports the
   retained pattern “source episodes learn coordinates; a small closed-form
   ridge solves the new task.” This is already present and should not be
   replaced by full-parameter MAML.
2. **Trusted, dependency-isolated validation controls weighting.** ATS, Ren,
   Meta-Weight-Net and the antibody study agree that weighting obtains its new
   information from a cleaner or more representative validation signal. In
   MetaSieve, “clean” must also mean document/scaffold/protein-component
   isolated and assay-semantic compatible.
3. **Support sensitivity is a primary outcome.** The candidate mechanism must
   reduce coefficient/prediction dispersion across natural support subsets,
   not merely average RMSE, while preserving correct over wrong/permuted
   support and correct over wrong protein.
4. **Reliability is not difficulty.** Loss and gradient alignment may enter a
   scheduler only beside explicit nuisance audits. Weight cannot be called
   biological or transferable until corruption and familiarity controls pass.

## 7. Mechanisms not transferable

- Full-model MAML, TaskNorm and TADAM-style task-conditioned feature
  modulation add high-dimensional support freedom without evidence that the
  retained low-dimensional solver is the bottleneck. They are rejected now.
- AdaMBind's public data split and CSV schema are weaker than MetaSieve's
  protein-homology, scaffold, document, support/query and recipient isolation.
- Generic support/query uniform noise is an augmentation baseline, not an
  assay-noise model. Query-noise training is especially uninformative under
  squared loss in expectation and contaminates the clean meta-objective.
- Low loss or positive support/query gradient cosine is not a proxy for
  protein-ligand mechanism. It can reward intercept, ligand familiarity or a
  mutually shared wrong-protein coordinate.
- Neither a scheduler nor a robust objective supplies the new observable that
  F-132 says is missing. Biological R0/R1 remains independently gated.

## 8. Candidate mechanisms (maximum three)

### Candidate A: cross-fitted task reliability/transferability audit

**Verdict: `ADOPT_WITH_MODIFICATION`, research-only.** Do not begin with a
larger scheduler. First produce offline, cross-fitted task scores from current
loss/alignment and audit them against homology, scaffold familiarity, sample
count, assay/source, affinity variance, family and baseline difficulty. Add
reliability/provenance only where measured. The score must change when correct
protein/support information is destroyed, and any eventual learner must use
dependency-isolated meta-validation components, a matched-capacity
permuted-statistics null and sampling without replacement.

This can improve which source gradients are trusted; it cannot admit biology.

### Candidate B: assay-aware support uncertainty in positive ridge

**Verdict: `ADOPT_WITH_MODIFICATION`, research-only.** Estimate source-only
uncertainty from assay type, document provenance and replicate disagreement.
Use it as precision weighting or a covariance/sensitivity analysis inside the
existing positive-ridge support solve. Do not perturb sealed queries and do
not pool Ki, Kd, Kdapp, IC50 or KIBA semantics. Compare clean, fixed-uniform,
Gaussian, heteroscedastic and replicate-derived arms only where each variance
has evidence.

The target outcome is lower section-coefficient and prediction instability
under realistic noisy support, with correct support/protein advantages intact.

### Candidate C: robust episodic section objective

**Verdict: `ADOPT_WITH_MODIFICATION`, conditional on A/B diagnostics.** Keep
the current `d<=5` representation and analytic solver. During source-only
episodic learning, penalize instability of the same query predictions or
section coefficients across two valid support resamples from the same target,
weighted by candidate B's uncertainty. This imports the support-sensitivity
lesson and the robust-direction idea from Eigen-Reptile without importing
Reptile parameter trajectories.

It must be compared with an equal-compute resampling-only control; otherwise a
gain may be ordinary data augmentation.

## 9. Cheapest falsification experiments

No full training is authorized by this report.

### Gate 0: static audit -- completed

- Gradient paths, detach points, source/meta-validation roles, clean query,
  without-replacement sampling, noise scale and deterministic streams were
  inspected.
- Current MetaSieve code implements its narrow first-order claims more cleanly
  than the official AdaMBind code, but its statistics do not identify biology.
- The existing fixed support-noise arm and current ATS are already negative.

### Gate 1: scheduler identifiability -- not authorized for confirmation

Run source/meta-validation only. Cache one candidate pool and compare uniform,
current ATS, loss-only, alignment-only, random and, only if justified, a
cross-fitted reliability learner. Re-score the identical pool after correct
protein to homology-matched wrong protein, correct to wrong-target support,
label permutation, ligand-only, intercept-only and destroyed-pair transforms.

Fail if ranking/benefit is preserved after biological destruction; if score is
mostly explained by homology, scaffold overlap, task size or baseline loss; or
if no learned arm beats its matched null with component-bootstrap uncertainty.

### Gate 2: uncertainty robustness -- not run

Before training, use the analytic ridge to propagate recorded support
uncertainty and measure coefficient/prediction spread across support resamples.
Only then compare source training arms. Report
`Delta_clean`, `Delta_noisy`, robustness slope and negative-transfer rate.
Fail if artificial-noise performance improves while clean unseen-target
performance or correct-support/protein contrasts decline.

### Gate 3: short end-to-end -- closed

Open only after Gates 1 and 2 pass on source/meta-validation components. Use
multiple seeds or component bootstrap and report target-macro CI, Spearman,
NDCG@10, RMSE, target-centered RMSE and negative-transfer rate, with all
support/protein/ligand/no-scheduler/no-noise controls. Confirmation requires a
fresh dependency-closed target cohort. No current result is fresh
confirmation for these candidates.

## 10. Recommended order

1. Freeze the current model and perform Candidate A's offline covariate and
   corruption audit. This is the cheapest way to falsify the scheduler story.
2. Inventory assay/replicate provenance. If uncertainty is not estimable,
   stop Candidate B rather than inventing a distribution.
3. Run analytic support-resampling and uncertainty propagation without source
   retraining. Stop if correct support/protein separation does not survive.
4. Only after those checks, preregister one short source/meta-validation run of
   A or B; Candidate C is conditional on a demonstrated residual instability.
5. Keep the F-132 biological replacement programme separate. Scheduler/noise
   evidence never authorizes R1 or V1 exact-pair integration.

## 11. Documentation changes

- `README.md` now lists AdaMBind as a core methodological reference, names the
  accepted cross-domain principles, and records the no-copy/no-admission
  boundary.
- `task.md` now records this audit, the current ATS/noise negative evidence,
  the three research-only candidates and Gate authorization state.

No code, theory, model checkpoint, dataset or experiment ledger was changed.

## 12. Final verdict

| Object | Verdict | Reason |
| --- | --- | --- |
| AdaMBind target-as-task and explicit support/query principle | `ADOPT` | structurally aligned and already retained |
| AdaMBind adaptive/noise ideas | `ADOPT_WITH_MODIFICATION` | useful hypotheses, but require isolated validation, realistic uncertainty and biological falsifiers |
| AdaMBind MAML architecture and public split as MetaSieve protocol | `REJECT` | excessive adaptation freedom and weaker dependency/assay closure |
| Official AdaMBind implementation | `REFERENCE_ONLY` | paper/code mismatches and ablation-path defects preclude direct reuse |
| Current MetaSieve ATS as biological transfer mechanism | `REJECT` | matched-null and all-k biological Gates failed |
| Fixed uniform support noise as robustness repair | `REJECT` | failed clean comparison and lacks assay semantics |
| Positive-ridge low-dimensional section | `ADOPT` | retained mathematical core; not itself proof of partner biology |
| Candidates A-C | `ADOPT_WITH_MODIFICATION` | research-only, sequentially gated, no fresh confirmation |

```text
CURRENT_ATS_NOISE_REPAIR_FAILED
TRAINING_REGULARIZATION_IS_NOT_BIOLOGICAL_IDENTIFIABILITY
NO_FRESH_CONFIRMATION_IN_THIS_AUDIT
GATE1_GATE2_GATE3_NOT_PASSED
R1_AFFINITY_FIELD_TRAINING_AUTHORIZED=false
V1_EXACT_PAIR_INTEGRATION_AUTHORIZED=false
PRODUCTION_MIGRATION_AUTHORIZED=false
```
