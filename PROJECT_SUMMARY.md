# MetaSieve Current Research Summary

This is the canonical snapshot of the active project state as of 2026-08-11.
`history.md` remains the complete chronology; formal failed experiments and
superseded runs are retained under `archive/retired_research_20260811/`.

## Objective and status

MetaSieve targets Cold Target few-shot drug-target affinity prediction. Each
target is a task, source learning uses support/query episodes, and evaluation
adapts from `k=1/2/3/5` measured support ligands without reading query labels.

The system is trainable and its GPU path is verified. It is not yet an
admitted biological or production model. The best retained development model
improves Cold Target prediction, but absolute quality, support-specificity and
partner-identity Gates remain failed.

```text
NO_VALIDATED_END_TO_END_FEWSHOT_DTA_MODEL
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
```

## Frozen and retained design

- Target-as-task source episodes and strict support/query separation.
- The uncentered positive-ridge Identifiable Meta-Section. K1 showed centered
  sections were worse on BindingDB, Davis and KIBA.
- Task freedom `d <= 5`, restricted to the support-observable row space.
- Frozen ESM2 protein and graph-derived ligand inputs.
- CUDA-batched episodes, dual ridge solves, validation and prediction. The
  vectorized implementation owns new checkpoints because float32 training
  trajectories are not bitwise invariant.
- The downstream operator `A(F,z)=K(B(z)F(z))` remains unchanged. No candidate
  biological statistic enters `z` without independent admission.

AdaMBind is a core methodological reference for target-as-task training,
support/query episodes and explicit robustness interventions. MetaSieve does
not adopt its full MAML architecture, public split, scheduler implementation or
label-noise implementation without matched falsification.

## Result ladder

| Stage | Result | Decision |
| --- | --- | --- |
| Main v0 | Real support adaptation, but wrong-protein cluster specificity failed | Retain meta-learning mechanics; replace biological pair representation |
| K1 | Centering degraded all three datasets | `KEEP_UNCENTERED_POSITIVE_RIDGE` |
| Cold Target V1 | Three-seed CUDA run; RMSE at k1/2/3/5 `1.582/1.429/1.363/1.313` | `COLD_TARGET_FEWSHOT_V1_NOT_YET_GOOD` |
| GPU vectorization | Full run 2066 s to 388 s (`5.32x`); inference difference <= `5.25e-6` | Retain; fresh checkpoints required |
| Targeted V1 repair | Support-only section + ligand MLP64 + `d=4,ridge=1`; RMSE `1.495/1.357/1.290/1.230`, CI `0.549/0.553/0.558/0.567` | Best development candidate; biological Gates still fail |
| R0-B exact geometry | 2,845 complexes and 26,044,068 exact cells, but underpowered confirmation | `R0B_NOT_RUN_FAIL_CLOSED` |
| R0-C fresh confirmation | Full RPS `.04052`, additive N2 `.03930`; NLL guard failed | `MARGINAL_OR_SLOT_RECALIBRATION_ONLY`; block R1/V1 integration |
| AdaMBind scheduler Gate | Gate0 passed; utility and matched-null/biology Gates failed every k | `REJECT_TASK_SCHEDULER_GATE1_FAIL_CLOSED` |
| Fixed support noise | Did not beat clean matched development arm | Reject as current repair |

## Active implementation surface

- `model/`: verified primitives and `MetaSieveV1`; scheduler/noise paths remain
  only as controlled development baselines.
- `research/meta_fewshot/`: corpus, Meta-Section, V0/V1 runner, selected V1
  studies, GPU benchmark and vectorization verifier.
- `research/crossed_interaction/` and `research/e0_identifiability/`: active
  dependencies and legacy baselines, not admitted biological modules.
- `research/correspondence_router/` and structure scripts: mapping, geometry,
  closure and pre-fit infrastructure reusable for a new observable. The failed
  R0-C field/trainer is archived.
- `report/meta_fewshot/main_v1_support_only_mlp64_d4_final/`: best retained
  development run.

The root `main.py` is the supported orchestration entry point and exposes only
active workflows:

```powershell
conda run -n drug python main.py status
conda run -n drug python main.py archive status
conda run -n drug python main.py verify tests
conda run -n drug python main.py verify v1-vectorization
conda run -n drug python main.py v1 train-evaluate --output report/meta_fewshot/<new-run>
conda run -n drug python main.py data prepare --help
conda run -n drug python main.py data verify --help
```

Formal V1 execution is CUDA-only. Output paths are explicit and existing
destinations fail closed so prior evidence is not silently overwritten.

## Archive policy

Formal negative or superseded evidence is archived, never discarded. This
includes the AdaMBind task-reliability experiment, failed R0-C pair
field/trainer, R0/R0B/R0C preregistrations/results and superseded V1 runs. The
consolidation manifest records original path, archive path, size and SHA256.

Reproducible ephemeral artifacts are removed from the active tree: pytest
caches, smoke runs, empty smoke directories and interrupted outputs without a
terminal RESULT. The current environment disallowed recursive deletion, so
they are retained under `ephemeral_quarantine/` rather than destroyed.

## Binding constraints and next work

The next scientific change must target a new affinity-directed,
partner-specific observable that cannot be explained by independent protein
and ligand marginals. It must beat capacity-matched additive, ligand-only and
wrong-partner controls before connection to V1. Opened R0-C results may inform
diagnosis but may not be reused as fresh confirmation.

Do not reconnect the archived residual to V1, tune against opened R0-C, restore
centered ridge/full MAML/test-time scheduling, or claim biological admission
from engineering speed or development-only results. Fixed synthetic noise is
not a proxy for missing assay uncertainty.

See `task.md` for requirements, `history.md` for chronology, `experiment.md`
for the experimental contract, and the archive manifest for evidence policy.

Post-consolidation verification: **159 active tests passed** in the `drug`
environment; all **229** formally archived files passed SHA256 revalidation.
