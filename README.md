# MetaSieve-DTA

MetaSieve is a trainable bioinformatics model for **few-shot drug-target
affinity prediction on unseen targets**. Its central learning problem is not
generic pocket detection: it must learn transferable target-ligand knowledge
from large open source datasets and adapt that knowledge from `k` measured
support ligands for a new target.

## Current status

The episodic few-shot stage reached its first precondition and stopped:

```text
FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE
```

The governed BindingDB Ki corpus supports episodic *training* (442 source
targets, 220 usable at `k=5`) but not unseen-target *evaluation*: only 16
held-out targets can carry `k=5`, below the declared minimum of 30, and the
resulting `MDE_d = 0.622` exceeds the declared `0.600` ceiling. Leakage is zero
on target, ligand, scaffold, document and protein-homology axes. No model has
been trained and target-coefficient heterogeneity remains untested.

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

The design is informed by Wan et al., [A meta learning and task adaptive
approach for drug target affinity prediction](https://www.nature.com/articles/s41467-026-70554-5)
(Nature Communications, 2026): targets define tasks and source knowledge is
learned through support/query episodes. MetaSieve does not copy its full MAML,
task-LSTM and label-noise stack. It replaces free inner-loop adaptation with a
small identifiable section so biology and the repository's mathematics remain
equally load-bearing.

## Biology and mathematics

- Frozen ESM2, ligand graph states and P1B geometry provide biological inputs.
- The currently auditable candidate is the 288D radial chemistry T-BASIS. It
  passed structural reconstruction and partner controls, but has not passed an
  affinity-admission Gate.
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
OPEN_BINDINGDB_QUOTIENT_TRAINING_PIPELINE_EXECUTABLE
CQ_R1_DEVELOPMENT_INTERACTION_OBSERVED
CQ_TBASIS_LINEAR_AFFINITY_WITNESS_NOT_OBSERVED
TARGET_COEFFICIENT_META_LEARNING_NOT_YET_TESTED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_FEWSHOT_DTA_MODEL
```

BindingDB Articles 202608 produced 12,457 governed Ki cells in 320 panels.
The first real training run tested one population-shared linear direction on
the 288D basis and explained only `0.000709` of development quotient variance;
correct-pair performance did not beat zero, foreign-ligand or deranged-protein
controls with a positive confidence bound. This rejects the shared direction,
not target-conditioned meta-learning. The next experiment must change only the
coefficient-sharing assumption and test a `d<=5` source-learned task subspace.

## Repository boundaries

- `theory/FINAL_FROZEN_THEORY/`: authoritative probability-law mathematics.
- `model/`: verified operator, encoder and geometry primitives; no assembled
  validated few-shot DTA model yet.
- `scripts/`: governed data, sealing, structure and training utilities.
- `research/crossed_interaction/`: current open-data training programme.
- `report/`: current status and compact evidence; terminated detail is archived.
- `history.md`: chronological decisions and failure lessons.

Read [task.md](task.md), [current status](report/CURRENT_RESEARCH_STATUS.md) and
the [evidence ledger](report/EXPERIMENTAL_EVIDENCE_LEDGER.md) first.

## Verification

```powershell
conda run -n drug python -m pytest -q
```

Large third-party releases, embedding banks and caches are not redistributed;
see `DATA_AVAILABILITY.md`.
