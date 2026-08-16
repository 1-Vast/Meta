# Meta-learning innovation gate

Date: 2026-07-31

## Decision

`STOP_PROTEIN_CONDITIONED_EXPANSION_NO_INNOVATION_ATTRIBUTABLE_GAIN`

The requested literature review and a new, leakage-safe candidate were run in
the `drug` CUDA environment. No model is promoted to the active runtime. The
failure is diagnostic rather than a capacity-tuning result: the smallest
protein-conditioned residual operator is worse than its protein-free control,
and the available crossed measurements do not identify a clean protein-
conditioned estimand.

## Literature result

The following primary sources were checked against their papers and public
implementations before designing the probe:

- AdaMBind (Wan et al., Nature Communications 2026,
  [DOI 10.1038/s41467-026-70554-5](https://doi.org/10.1038/s41467-026-70554-5))
- MAML (Finn et al., ICML 2017,
  [arXiv:1703.03400](https://arxiv.org/abs/1703.03400))
- TADAM (Oreshkin et al., NeurIPS 2018,
  [arXiv:1805.10123](https://arxiv.org/abs/1805.10123))
- Deep Kernel Transfer (Patacchiola et al., NeurIPS 2020,
  [arXiv:2008.05414](https://arxiv.org/abs/2008.05414))
- MetaDTA public implementation (commit `d348c033ad153a1d9a4eeb28dfaaf25b547ec729`)

AdaMBind adds a query-loss/gradient-similarity task scheduler around a
conventional MAML DTA model. It does not add crossed labels or a new
protein-ligand interaction operator. Its random/CD-HIT split still permits
within-target support/query chemical and assay overlap, so its reported gain is
not evidence for FORT's strict target/homology/scaffold/document/assay-cold
estimand. TADAM changes task-conditioned metric scaling; MAML changes inner
updates; Deep Kernel Transfer changes covariance/posterior adaptation. All
three require a transferable task signal that has not passed FORT's
identifiability gate. A scheduler cannot create missing crossed observations.

## Data audit

The corrected audit is model-free and TRAIN-only. pKi and pKd remain separate;
rectangles are collapsed to target-pair/document units rather than counted as
independent rows.

| Run | pKi units | pKd units | pKi reversal | pKd reversal | noise source |
| --- | ---: | ---: | ---: | ---: | --- |
| raw, strict cross-homology | 12,661 | 26,153 | 0.360 | 0.400 | within-cell replicate |
| registry-closed + `[0,14]` clip | 12,661 | 22,640 | 0.360 | 0.413 | registry `replicate_sd` |

The registry-closed run contains pKi `248,775` and pKd `25,625` raw rows. It
explicitly excludes `4,635` pKd rows removed by the registered pKi-first
target/ligand collision policy. In this conservative run:

- exact document-local rectangles remain abundant, but `same_assay_rectangles = 0`;
- DD/noise(q90) is about `0.80` for pKi and `0.30` for pKd;
- 27.3% of raw assay-panel pair combinations in the raw sensitivity are in
  multi-panel units, so median unit aggregation can mix assay panels;
- the largest coarse document source family accounts for about 69% of units;
- strict runs exclude target pairs sharing a registered homology component.

The zero same-assay count is treated as a ChEMBL assay-ID limitation, not as a
biological null: assay IDs are target-specific in this extract. Protocol,
condition, and provenance comparability are not established. The raw audit
also has 220 pKi and 2 pKd out-of-range records; clipping is reported only as a
sensitivity.

## Candidate and strict result

`research/residual_bilinear_probe.py` implements IDG-RBP, an explicitly
bounded diagnostic candidate:

- fresh fit-only normalization and seeded Gaussian projections;
- a 16 x 16 bilinear matrix fit to residual ligand differences on fit
  components only;
- exact antisymmetric score `u_t^T B (x_q - x_s)`;
- strict gate closure on scaffold, ligand connectivity, document, and assay;
- B0, calibration/protein-free, correct protein, homology-matched wrong protein
  (when available), and random wrong protein controls;
- no old checkpoint, target-ID feature, structure, development label, or sealed
  label.

Full pKi result (`residual_bilinear_probe_pKi_seed1729.json`, 58 strict gate
episodes):

| Arm | RMSE | Spearman | Pairwise |
| --- | ---: | ---: | ---: |
| B0 | 1.3598 | 0.0865 | 0.5334 |
| protein-free calibration | 1.3550 | 0.0865 | 0.5334 |
| correct protein | 1.4266 | 0.0520 | 0.5182 |
| random wrong protein | 1.3965 | 0.0946 | 0.5359 |
| matched wrong protein | 2.1465 | -0.4382 | 0.3075 |

The correct-minus-protein-free component-bootstrap RMSE gain is
`-0.0709 [-0.1221, -0.0200]`; MAE gain is
`-0.0613 [-0.1088, -0.0099]`. Correct-minus-matched-wrong RMSE is positive,
but only 5 matched targets across 2 components are available and all ranking
gates fail. The smoke configuration independently produced correct RMSE
`1.5364` versus protein-free `1.4174`, with gain `-0.1190
[-0.2499,-0.0308]`.

This rules out using a deeper cross-attention block, graph ligand encoder,
TADAM metric, AdaMBind scheduler, or more epochs as a defensible rescue. The
new mechanism did not produce the required performance, and its protein effect
is not attributable under the strict controls.

## Checkpoint and closure correction

The earlier AnchorDelta numbers are retained only as non-decisional diagnostics.
The old checkpoint was trained on broad TRAIN episodes without a target/component
provenance manifest before a new holdout was carved, and the smoke selector did
not close support/query documents and assays. `scripts/anchor_delta.py` now
defaults to a fresh fit-component encoder and rejects uncertified checkpoints;
the new probe does not load a checkpoint at all.

## Reopening rule

Do not alter width, epochs, scheduler, loss weights, or representation to chase
this failed gate. Reopen only after one of the following is supplied:

1. a source/assay/protocol-comparable crossed panel with independent provenance
   blocks; or
2. a new authorized interaction label with the same target/homology/scaffold/
   document/assay closure.

Then require, before any high-capacity model:

- correct-protein IDG-RBP beats protein-free and homology-matched wrong protein;
- RMSE, Spearman, and pairwise component-bootstrap lower bounds are positive;
- the effect survives provenance-family and ligand-similarity bins;
- a matched additive/random-feature null fails;
- multi-seed results are consistent.

Only after those gates pass may a TADAM-style metric or AdaMBind-style scheduler
be evaluated as an equal-budget secondary ablation. Until then, claiming
excellent protein-conditioned performance would be attributing a result to an
innovation that the data and controls do not support.
