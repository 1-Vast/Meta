# MetaSieve-DTA

The canonical current-state overview is
[`PROJECT_SUMMARY.md`](PROJECT_SUMMARY.md). This README describes the project;
the summary records the binding result ladder and active/archive boundary.

MetaSieve is a trainable bioinformatics model for **few-shot drug-target
affinity prediction on unseen targets**. Its central learning problem is not
generic pocket detection: it must learn transferable target-ligand knowledge
from large open source datasets and adapt that knowledge from `k` measured
support ligands for a new target.

## Current status

MetaSieve-main v0 established useful support adaptation, but its biological
specificity failed the cluster-level wrong-protein audit.  K1 then tested the
proposed calibration repair across BindingDB, Davis and KIBA:

```text
KEEP_UNCENTERED_POSITIVE_RIDGE
CENTERED_SECTION_CROSS_DATASET_FAIL
BIOLOGICAL_PAIR_REPRESENTATION_IS_THE_NEXT_REPLACEABLE_MODULE
```

The production candidate therefore remains the original uncentered,
strictly-positive ridge section. A fresh 219-component R0-C confirmation found
that the exact residue--atom residual strongly improves the frozen distance
prior but is worse than a capacity-matched additive atom/residue null. Its
terminal verdict is `MARGINAL_OR_SLOT_RECALIBRATION_ONLY`; R1 affinity training
and migration into V1, `model/` or `z` are not authorized.

## Core task

Each protein target is a meta-learning task. Source training uses target-wise
support/query episodes; evaluation holds out the target family, ligand
scaffold and source document. The primary support sizes are `k=1/2/3/5`.

The minimal intended predictor is

```text
m(P,L) = U^T phi(P,L),                         d <= 5
y_hat  = f_L(L, endpoint) + w0^T m(P,L) + a_t^T m(P,L)
```

where `phi` is a frozen, audited biological interaction measurement and `U`,
`w0` are learned across source targets. The target section `a_t` is estimated
from the support set with a strictly positive ridge penalty and is restricted
to support-observable directions. This makes the meta-learner trainable while
keeping continuous task freedom no larger than the support rank.

The core methodological reference is Wan et al., [A meta learning and task
adaptive approach for drug target affinity prediction](https://www.nature.com/articles/s41467-026-70554-5)
(Nature Communications, 2026; [official AdaMBind
code](https://github.com/Moohyun-w/AdaMBind)). It establishes three useful
principles: targets are tasks, transfer is trained through support/query
episodes, and source-task scheduling and training-label perturbation can be
studied as explicit robustness interventions. MetaSieve retains its stricter
closed-form identifiable section in place of full MAML adaptation.

AdaMBind is a core reference, not an unchecked dependency. The completed
mechanism/code audit found that its transferable principles are task-adaptive
training, explicit target-as-task support/query episodes and robust-label
hypothesis testing, not its full MAML architecture or public split. MetaSieve's
protein-family, ligand-scaffold, source-document, assay, support/query and
recipient isolation remain stricter. The existing MetaSieve ATS and fixed
support-noise arms failed their matched controls and are retained as negative
baselines, not production components.

The current biological-bridge reference also includes Li et al.,
[UniPert-G2CP bridges genetic and chemical screens from molecular
representation to phenotype modeling](https://doi.org/10.1016/j.cell.2026.06.005)
(Cell, 2026). MetaSieve does not import UniPert-G2CP as a target-ligand
similarity model. The borrowed architecture is narrower: a two-stage bridge in
which external cross-modal supervision first makes a perturbagen or interaction
coordinate identifiable, then downstream transfer predicts the task-specific
effect. In MetaSieve terms this is implemented as a target-conditioned
ligand-pair SAR-delta bridge before any few-shot section update:

```text
UniPert-G2CP:
genetic / chemical perturbagen representation -> genetic-to-chemical phenotype transfer

MetaSieve translation:
source target ligand-pair SAR-delta supervision -> cold-target affinity/ranking transfer
```

The registered BindingDB retest is
`report/crossed_interaction/unipert_g2cp_sardelta_bridge_gate1_20260812`.
It passed its local bridge Gate: 20,423 train pairs, 1,033 development pairs,
8 development dependency components, correct MSE `0.233714` versus zero-delta
MSE `0.580072`, component reduction `+0.376740`, and LCB95 `+0.270172`.
The stricter attribution Gate
`report/crossed_interaction/bindingdb_sardelta_attribution_gate1_20260812`
then showed that this PASS must not be attributed to a true UniPert-style
target x chemical-transformation bridge yet. The additive concat arm
`[protein; ligand_delta]` beat zero, but it violated SAR-delta antisymmetry
(`max |f(i,j)+f(j,i)| = 1.654930`). The bilinear target x ligand-delta arm was
antisymmetric, but failed against ligand-delta, wrong-target and shuffled-target
controls. The correct current conclusion is therefore: SAR-delta transfer
signal is identified; target-specific conditioning and end-to-end Cold Target
utility remain open; scalar and unordered panel-edge lifts are rejected.

The screened methodological references also include [adaptive task
scheduling](https://proceedings.neurips.cc/paper/2021/hash/3dc4876f3f08201c7c76cb71fa1da439-Abstract.html),
[differentiable closed-form solvers](https://openreview.net/forum?id=HyxnZh0ct7),
[support-set sensitivity](https://proceedings.neurips.cc/paper/2021/hash/ab73f542b6d60c4de151800b8abc0a6c-Abstract.html)
and [trusted-data reweighting in biological
engineering](https://doi.org/10.1016/j.cels.2023.12.003). Their accepted roles
are scheduler falsification, representation/solver separation,
support-stability diagnostics and assay-aware uncertainty. None is an admitted
production component. The completed audit and its decisions are consolidated
in [`PROJECT_SUMMARY.md`](PROJECT_SUMMARY.md); the full report is archived.

Candidate A's formal follow-up is now complete. A hash-bound, three-seed CUDA
audit fitted a two-input difficulty/alignment scorer out of whole CD-HIT40
components and compared it with an equal-capacity stratified-permutation null.
Nuisance covariates never entered the scorer; empty scaffolds were treated as
missing. Clean score-to-transfer-utility Spearman was
`-0.060/-0.163/0.010/-0.089` at `k=1/2/3/5`, and matched-null, correct-protein,
wrong-support and shortcut criteria did not pass all k. The scheduler lane is
therefore `REJECT_TASK_SCHEDULER_GATE1_FAIL_CLOSED`; `train_main_v1.py` was not
changed and no short training Gate was authorized. The full report and failed
implementation are archived; the binding verdict is in `PROJECT_SUMMARY.md`.

The executable Cold Target V1 development runner is exposed through
`python main.py v1 train-evaluate --output report/meta_fewshot/<new-run>` and
is CUDA-only. It uses a
paired factorial design (uniform/ATS task selection x clean/support-noisy
source episodes), an equal-capacity permuted-statistics scheduler null, and
fixed `k=1/2/3/5` evaluation episodes. The current BindingDB development
corpus has complete CD-HIT-40 target-cluster separation, but not ligand-
scaffold or source-document closure; its result is therefore development
evidence only, not a fresh confirmation result.

The 2026-08-11 formal CUDA development run used three seeds, produced 24
hash-bound checkpoints and 1,251,600 label-free prediction rows, and selected
`uniform_clean` on meta-validation. On the shared 33-target Cold Target cohort,
target-macro RMSE/CI were `1.582/0.529`, `1.429/0.534`, `1.363/0.538`, and
`1.313/0.536` for `k=1/2/3/5`. ATS did not beat its matched null across all k,
support-noise did not beat clean training across all k, and the support-
specificity and wrong-protein Gates failed. The governed verdict is therefore
`COLD_TARGET_FEWSHOT_V1_NOT_YET_GOOD`; V1 is trainable and connected, but its
current representation is not admitted as a good Cold Target predictor.

The subsequent CUDA vectorization batches candidate episodes, dual ridge
solves, per-task ATS gradients, validation and inference without changing the
Meta-Section equation. A full three-seed reproduction completed in 388 seconds
instead of 2,066 seconds (`5.32x` faster); same-config per-model speedups were
`11.52x` for uniform training and `5.09x` for ATS. Frozen-checkpoint predictions
matched within `5.25e-6`. Float32 optimization was not trajectory-bitwise—the
selected arm changed to `ats_clean`—so the vectorized implementation owns new
checkpoints. The terminal scientific verdict and failed biological Gates did
not improve.

Targeted V1 repairs then changed one axis at a time. A support-only
Meta-Section removed the pair-dependent zero-shot anchor, a 64-wide ligand
population MLP improved global calibration, and meta-validation selected
`d=4, ridge=1`. The resulting CUDA run reached RMSE
`1.495/1.357/1.290/1.230` and CI `0.549/0.553/0.558/0.567` at `k=1/2/3/5`, a
consistent `7.6-8.9%` RMSE reduction over vectorized V1. A nonlinear pair
encoder regressed, wider population networks were not selected, and source
population pretraining selected zero steps. Absolute quality,
support-specificity and partner-identity Gates still fail, so this is a
retained development improvement, not a good or admitted model.

R0-B exact geometry is complete for 2,845 governed complexes with zero mapping
exclusions and 26,044,068 exact atom-residue cells. Its pre-fit audit found
structural headroom (`S_prior=0.12322`, `S_add*=0.02296`), but 53 heldout
protein components give `MDE80=0.00746`, above the frozen
`delta*=0.00616`. R0-B therefore stopped before training with
`R0B_NOT_RUN_FAIL_CLOSED`; no Gate was relaxed.

R0-C then supplied a genuinely new, geometry-blind and chemistry/protein-closed
confirmation lane. Its 219 independent systems passed pre-fit power
(`MDE80=0.00568 < delta*=0.00637`) and enabled the registered three-seed CUDA
run. Full RPS improved from the frozen prior's `0.12749` to `0.04052`, but the
capacity-matched additive null reached `0.03930`; the paired Full advantage was
negative (`-0.00122`). The field therefore learned useful marginal distance
recalibration, not an admitted exact pair interaction.

## Biology and mathematics

- Frozen ESM2, ligand graph states and P1B geometry provide biological inputs.
- The 288D radial chemistry T-BASIS passed its structural Gate but failed
  affinity/selectivity admission; it is retained as a legacy baseline.
- The replacement candidate keeps exact residue and ligand-atom identities,
  uses P1B only as a distance prior, and exposes six typed local interaction
  channels. It remains research-only until real correct-partner Gates pass.
- Open datasets are used according to measurement semantics: Ki/Kd/Kdapp are
  not pooled, and inhibition/displacement panels provide ordinal rather than
  absolute-affinity supervision.
- No raw pair map or arbitrary neural latent may enter the mathematical state
  `z`. A biological statistic must first beat ligand-only and wrong-partner
  controls, replicate independently, and remain identifiable from support.
- The authoritative downstream operator remains unchanged:

```text
A(F,z) = K(B(z)F(z))
```

The rank bound is a linear-algebra property of the section design; it is not
retroactively claimed as a theorem of `FINAL_FROZEN_THEORY`.

## Current evidence

```text
REAL_META_SECTION_META_EFFECT_IDENTIFIED
WRONG_PROTEIN_SPECIFICITY_NOT_IDENTIFIED_ACROSS_CDHIT_CLUSTERS
TBASIS_SELECTIVITY_SIGNAL_NOT_IDENTIFIED
KEEP_UNCENTERED_POSITIVE_RIDGE
CENTERED_SECTION_CROSS_DATASET_FAIL
R0_R1_EXACT_PAIR_SOFTWARE_AND_SYNTHETIC_CONTRACT_PASS
R0C_PREFIT_ADMISSION_PASS
R0C_EXACT_PAIR_INCREMENTAL_FAIL
MARGINAL_OR_SLOT_RECALIBRATION_ONLY
SAR_DELTA_TRANSFER_SIGNAL_IDENTIFIED
TARGET_CONDITIONING_NOT_YET_IDENTIFIED
UNIPERT_STYLE_INTERACTION_BRIDGE_NOT_YET_PROPERLY_TESTED
SCALAR_AND_UNORDERED_EDGE_LIFTS_REJECTED
END_TO_END_COLD_TARGET_UTILITY_OPEN
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_FEWSHOT_DTA_MODEL
```

## Repository boundaries

- `theory/FINAL_FROZEN_THEORY/`: authoritative probability-law mathematics.
- `model/`: verified operator, encoder and geometry primitives; no assembled
  validated few-shot DTA model yet.
- `scripts/`: governed data, sealing, structure and training utilities.
- `research/crossed_interaction/`: current open-data training programme.
- `report/`: active experiment outputs only; historical evidence is archived.
- `history.md`: chronological decisions and failure lessons.

Read [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md), then [task.md](task.md) and
[history.md](history.md).

## Verification

```powershell
conda run -n drug python main.py verify tests
```

Large third-party releases, embedding banks and caches are not redistributed;
see `DATA_AVAILABILITY.md`.
