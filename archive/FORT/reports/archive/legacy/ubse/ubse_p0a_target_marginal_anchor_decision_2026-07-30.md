# UBSE-P0A Target-Marginal Anchor Decision

**Decision date:** 2026-07-30  
**Accepted artifact:** corrected epoch-4 fixed-seed control ledger  
**Decision:** `FREEZE_UBSE_P0A_FOR_A1_POCKET_PROPOSAL_ONLY`

## Outcome

UBSE-P0A passes all six preregistered gates after the frozen
weights-only shuffle-control correction. The result establishes a
source-closed, protein-only target-marginal pocket proposal. It does not
establish ligand-conditioned interaction information, binding affinity,
causality, or a deployable strict dual-cold affinity predictor.

The frozen P0A output may be used only as a ranking/proposal input for the
later A1 typed-event student. It does not change the UBSE-A0 remote
coordinate wait and does not unlock affinity, confirmation, or sealed data.

## Execution history

Attempt 1 was infrastructure-aborted without a scientific artifact or
decision. It was still at 100% GPU utilization after 57.3 minutes, then the
foreground execution session and process disappeared before the next check.
No checkpoint, ledger, result, traceback, temporary output, or Windows crash
event existed. The disposition is recorded separately as:

`ABORT_UBSE_P0A_ATTEMPT1_EXECUTION_SESSION_LOSS_NO_SCIENTIFIC_DECISION`

Attempt 2 used identical code, data, seeds, optimizer, controls, and gates in
a hidden detached process with persistent stdout/stderr. It completed in
23,698.993 seconds. Standard error was empty, and standard output contained
the complete raw gate result.

## Frozen substrate

| Quantity | Result |
|---|---:|
| source rows before closure | 66,660 |
| source targets before closure | 40,969 |
| held exact targets | 152 |
| held PubMed units | 152 |
| held scaffolds | 152 |
| homology targets excluded | 1,890 |
| retained rows after union closure | 62,849 |
| retained targets after union closure | 38,781 |
| parsed valid training targets | 38,781 |
| invalid training targets | 0 |
| training residues | 13,320,077 |
| validation panels | 64/64 |
| untouched audit panels | 88 |

Validation and audit have zero overlap on exact target, homology cluster,
PubMed, and scaffold. The corrected validation ledger contains exactly the
64 frozen validation targets and no audit target.

## Model execution

- Backbone: `facebook/esm2_t6_8M_UR50D`
- Revision: `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`
- Seeds: 1729, 1730, 1731
- Epochs: 4
- Evaluation epochs: 1, 2, 4
- Training windows per seed: 39,440
- Device: `cuda`
- GPU: NVIDIA GeForce RTX 4060 Laptop GPU
- PyTorch: `2.6.0+cu124`
- Precision: CUDA mixed float16
- CPU model offload: none
- Finite runs: 3/3

Mean training loss fell consistently for all three seeds:

| Seed | epoch 1 | epoch 2 | epoch 3 | epoch 4 |
|---:|---:|---:|---:|---:|
| 1729 | 0.4772 | 0.3798 | 0.3264 | 0.2889 |
| 1730 | 0.4801 | 0.3814 | 0.3295 | 0.2923 |
| 1731 | 0.4736 | 0.3792 | 0.3271 | 0.2897 |

## Shuffle-control correction

The raw process used `model_seed + epoch` for each diagnostic shuffled
sequence. The preregistration requires one deterministic shuffle per model
seed. This deviation was detected and frozen before any seed weight or result
existed.

The accepted correction:

- loaded the three final raw checkpoints;
- performed no retraining;
- reevaluated only the 192 epoch-4 shuffled-sequence rows with
  `control_seed == model_seed`;
- retained all epoch-1/2 trajectory rows as diagnostics;
- left all other 2,112 raw metric rows equal at the data-frame value level;
- recomputed every gate from the corrected 2,304-row ledger;
- preserved the raw result and ledger.

Independent ledger checks found:

| Check | Result |
|---|---:|
| corrected rows | 2,304 |
| epoch-4 rows | 768 |
| seed/epoch/control cells | 36 |
| rows per cell | 64 |
| duplicate primary keys | 0 |
| corrected primary shuffle rows | 192 |
| unaffected metric rows identical | 2,112/2,112 |
| nonfinite metric values | 0 |
| validation targets / panels | 64 / 64 |
| audit-target overlap | 0 |

## Accepted metrics

The accepted values are the median across the three model seeds.

| Control | AP | AUROC | oracle-size top-k recall | soft BCE |
|---|---:|---:|---:|---:|
| correct sequence | 0.3159 | 0.8398 | 0.2841 | 0.3695 |
| fixed-seed shuffled sequence | 0.0645 | 0.4975 | 0.0435 | 0.4494 |
| fit-only residue propensity | 0.1018 | 0.6578 | 0.0843 | 0.1365 |
| constant position | 0.0391 | 0.5000 | 0.0391 | 0.3000 |

The target-bootstrap AP contrasts are:

| Contrast | Estimate | 95% interval |
|---|---:|---:|
| correct - propensity | 0.2147 | [0.1604, 0.2692] |
| correct - fixed-seed shuffle | 0.2551 | [0.1971, 0.3120] |

Correct-sequence AP by seed was:

- seed 1729: 0.3149;
- seed 1730: 0.3159;
- seed 1731: 0.3188.

The AP range was 0.0040. Correct-minus-shuffle AP was positive for every
seed: 0.2504, 0.2512, and 0.2639.

## Gate decision

| Gate | Frozen requirement | Result | Pass |
|---|---|---|---|
| P0A-1 substrate | exact 62,849 rows/38,781 closed targets; >=35,000 valid targets; 64 validation panels | 62,849 / 38,781 / 38,781 / 64 | yes |
| P0A-2 absolute anchor | AP >=0.25; AUROC >=0.75; top-k recall >=0.25 | 0.3159 / 0.8398 / 0.2841 | yes |
| P0A-3 propensity increment | AP delta >=0.05; LCB95 >0 | 0.2147; LCB 0.1604 | yes |
| P0A-4 sequence destruction | AP delta >=0.05; LCB95 >0; all seeds positive | 0.2551; LCB 0.1971; 3/3 | yes |
| P0A-5 stability | AP range <=0.05; all finite | 0.0040; 3/3 finite | yes |
| P0A-6 execution/firewall | CUDA, hashes, no CPU model offload, no forbidden input/outcome | all satisfied | yes |

## Firewall

The run reports all of the following as false:

- ligand features loaded;
- affinity fields loaded;
- affinity values decoded;
- coordinates loaded;
- development/confirmation features loaded;
- development/confirmation labels loaded;
- sealed test consumed;
- audit contact labels loaded.

The source, panel, raw ledger, corrected ledger, and all three weight hashes
were verified. The corrected result binds the exact executed P0A code,
preregistration, correction protocol, and correction code hashes.

## GPU utilization

Across 11,279 telemetry samples:

| Measure | Mean | Median | p95 | Max |
|---|---:|---:|---:|---:|
| GPU utilization (%) | 98.40 | 100.00 | 100.00 | 100.00 |
| memory used (MiB) | 7,774 | 7,864 | 7,929 | 7,945 |
| power (W) | 32.66 | 26.96 | 71.15 | 97.11 |
| temperature (C) | 52.03 | 51.00 | 61.00 | 66.00 |

The device stayed highly utilized and thermally safe. Peak memory left only
243 MiB of reported framebuffer headroom, so later A1 workloads must retain
the explicit token/OT-cell budget and should not silently increase batch
size.

## Interpretation boundary and residual risk

P0A answers one narrow question: ordered protein-sequence context carries a
source-closed target-marginal pocket signal beyond residue propensity,
sequence composition, and fixed position. It does not answer which ligand
functional group interacts with which residue, nor whether a pair has a
particular affinity.

The model is a ranking proposal rather than a calibrated contact
probability. Its correct-sequence soft BCE (0.3695) is worse than the
fit-only propensity baseline (0.1365), even though ranking AP, AUROC, and
top-k recall pass. A1 must consume frozen ensemble logits/ranks, not interpret
them as calibrated probabilities or physical energies.

The oracle-size top-k recall passes by 0.0341, but this does not establish
the separate A1 teacher-event proposal-recall gate. A1 must still demonstrate
at least the preregistered event-residue recall on fit, validation, and audit.

P0A therefore advances the program from "no deployable target anchor" to
"frozen target-marginal proposal available." The central information deficit
after UBSE-G1 remains unresolved until a source-closed real
residue-functional-group-event teacher is transported to deployment inputs
and passes ligand/protein/structure destruction plus the exact-null rectangle
tests.

## Authoritative artifacts

- `reports/active/ubse_p0a_target_marginal_anchor_preregistration_2026-07-29.md`
- `reports/active/ubse_p0a_attempt2_detached_execution_amendment_2026-07-29.md`
- `reports/active/ubse_p0a_shuffle_control_execution_correction_2026-07-29.md`
- `reports/active/ubse_p0a_seed1729_1731.json`
- `reports/active/ubse_p0a_validation_ledger.parquet`
- `reports/active/ubse_p0a_seed1729_1731_corrected.json`
- `reports/active/ubse_p0a_validation_ledger_corrected.parquet`
- `reports/active/ubse_p0a_seed1729.pt`
- `reports/active/ubse_p0a_seed1730.pt`
- `reports/active/ubse_p0a_seed1731.pt`

Accepted corrected result SHA-256:

`d4be1fe9fa1a87faa5db8f390587f742a629d634ef35d9364e89cb92163a7a61`

Accepted corrected ledger SHA-256:

`81583fc6127e93e40207ecfdc97eab283e611b6647d3ebd273def1b84dd46587`
