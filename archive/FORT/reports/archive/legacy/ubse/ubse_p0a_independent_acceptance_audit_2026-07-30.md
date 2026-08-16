# UBSE-P0A Independent Acceptance Audit

**Date:** 2026-07-30  
**Mode:** read-only independent reproduction  
**Verdict:** accept

`FREEZE_UBSE_P0A_FOR_A1_POCKET_PROPOSAL_ONLY`

No defect was found that invalidates the corrected P0A decision. This
acceptance is limited to protein-side pocket ranking and does not establish
ligand-conditioned information, affinity prediction, or an independent
confirmation result.

## Frozen order

- preregistration: 2026-07-29 18:46;
- shuffle-correction protocol: 2026-07-29 19:22;
- detached attempt-2 amendment: 2026-07-29 20:03;
- attempt-2 launch: 2026-07-29 20:04;
- first seed weight: 2026-07-29 21:17;
- raw result: 2026-07-30 02:39;
- corrected result: 2026-07-30 02:41.

Code, protocol, source, panel, raw/corrected ledger, raw result, and all
three checkpoint hashes match their recorded identities.

## Independent reconstruction

- Closure reproduced `66,660 -> 62,849` rows and
  `40,969 -> 38,781` targets.
- Training has zero exact target, PubMed, or scaffold overlap with the 152
  held units.
- The full source has zero `target_key -> sequence` conflicts.
- All 64 validation panels contain exactly the ligand connectivity set
  frozen in the panel manifest.
- Validation and audit have zero target, homology, PubMed, or scaffold
  overlap.
- Training reconstructed to 38,781 valid targets, zero invalid targets,
  13,320,077 residues, and 39,440 windows per seed.
- Every batch respected the 64-window, 12,000-token, and 4,000,000-attention
  cell limits.
- Each checkpoint has 110 tensors and 7,737,770 finite parameters with the
  correct model revision and seed metadata.

## Correction reproduction

The raw and corrected ledgers contain 2,304 rows, 36 complete
seed/epoch/control cells, 64 targets per cell, and zero duplicate keys.
Only the 192 epoch-4 shuffled-sequence rows were reevaluated. The other 2,112
metric rows are equal.

All 192 accepted shuffle rows use `control_seed == model_seed`; all epoch-1/2
shuffle rows retain the explicitly diagnostic `model_seed + epoch` protocol.
All 576 shuffled-input hashes reconstruct from the frozen seeds.

The auditor loaded each checkpoint and reran:

- 192 epoch-4 correct-sequence examples;
- 192 epoch-4 fixed-model-seed shuffled examples.

Every metric matched the corrected ledger exactly; maximum numerical error
was zero. The correction did not make the control artificially easier: the
correct-minus-shuffle point estimate changed from raw 0.2560 to corrected
0.2551.

## Recomputed gates

| Quantity | Result |
|---|---:|
| AP | 0.315853 |
| AUROC | 0.839754 |
| top-k recall | 0.284139 |
| AP minus propensity | 0.214695 |
| propensity contrast LCB95 | 0.160384 |
| AP minus fixed-seed shuffle | 0.255132 |
| shuffle contrast LCB95 | 0.197060 |
| seed shuffle deltas | 0.250353 / 0.251160 / 0.263885 |
| seed AP | 0.314856 / 0.315853 / 0.318841 |
| seed AP range | 0.003986 |

The bootstrap reproduction matches the corrected JSON exactly.

GPU telemetry contains 11,279 samples: mean utilization 98.40%, median/p95/max
100%, peak memory 7,945/8,188 MiB, peak power 97.11 W, and peak temperature
66 C. Nominal telemetry coverage is approximately 95.2% of the 6.583-hour
wall time.

## Residual risks

These risks do not invalidate the frozen proposal role:

1. The 64-target validation was already involved in the G1 research route.
   It is development validation, not an independent confirmation set.
2. ESM2 pretraining membership was not closed against held targets.
   "Source-closed" applies to the BioLiP contact supervision, not the
   foundation-model corpus.
3. Correct-sequence soft BCE 0.3695 is worse than propensity 0.1365. P0A is a
   rank/proposal model; its logits are not calibrated probabilities or
   energies.
4. The top-k gate margin is only 0.0341. A1 must independently pass its
   teacher-event proposal-recall gate.
5. The executed code did not explicitly compare validation ligand
   connectivity sets, although this run's 64/64 sets were independently
   confirmed.
6. Some firewall/CUDA fields are declarative. Static column-path review and
   checkpoint reexecution support them, but a future version should add
   runtime file-access auditing.
7. The 6.58-hour runtime exceeded the original estimate. This is an
   execution-estimation error, not a scientific gate failure.

The accepted use remains frozen ensemble ranking/proposal only.
