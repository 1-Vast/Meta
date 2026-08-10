# Phase 2B S4R evidence consolidation

## Terminal result

```text
S4R-A ligand representation audit ......... GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE
S4R single-axis graph-aware transfer ...... REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED
heldout-B .................................. NOT CREATED AND NOT READ
R6 amplitude/B5 integration ................ NOT RUN
affinity value reads ....................... 0
```

S4R changed exactly one axis of the failed S3R experiment: the ligand
statistic. Everything else — frozen ESM2 residue states, the gauge, the direct
`W` ordinal estimator, the unit-Frobenius constraint, the pairwise ordinal
loss, the hierarchical pair/construct/component weighting, the closure split,
the 210-update stream, the seeds, the control maps and the R1-R5 margins — is
byte-identical to S3R and is proved so below.

## Part 1 — the representation audit

The mean-pooled 41-D atom-marginal basis is measurably collapsed. On 39,435
label-blind heldout-A ligand-graph pairs:

| statistic | mean-pooled 41-D |
|---|---:|
| effective rank of the embedding | 5.336 of numerical rank 33 |
| effective rank of the pair-difference matrix | 6.183 |
| distinct ligand graphs sharing a bit-identical vector | 687, in 307 groups |
| variance of the difference norm explained by heavy-atom-count difference | 0.8517 |

The collapse has an exact chemical witness. `m`-xylene and `p`-xylene have
identical atom-local feature multisets and therefore **bit-identical** 41-D
means; only their connectivity differs. That case is asserted in
`tests/test_s7_l2b_phase2b_s4r.py`.

Six frozen graph-aware candidates were measured. All six cleared every
preregistered A-gate:

| radius | d | eff. rank of differences | INC | RET | coverage | `W` parameters |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 128 | 20.93 | 0.3552 | 0.0780 | 0.99997 | 163,840 |
| 1 | 256 | 28.44 | 0.3910 | 0.0514 | 1.00000 | 327,680 |
| 1 | 512 | 35.08 | 0.4143 | 0.0262 | 1.00000 | 655,360 |
| 2 | 128 | 21.78 | 0.3799 | 0.0802 | 0.99997 | 163,840 |
| 2 | 256 | 30.85 | 0.4236 | 0.0510 | 1.00000 | 327,680 |
| 2 | 512 | 41.09 | 0.4588 | 0.0283 | 1.00000 | 655,360 |

`INC` is the fraction of candidate pair-difference energy that no linear
function of the baseline pair difference can express; `RET` is the converse
loss. The registered capacity-parsimony rule selected the smallest admissible
representation: radius 1, `d = 128`, per-heavy-atom Morgan environment counts
over a train-only vocabulary, `163,840` parameters against the baseline's
`52,480`.

The audit verdict was therefore `GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_
INFORMATIVE`. The representation is not the excuse.

## Part 2 — the primary panel

46,818 unordered ligand pairs across 112 closure components. Every arm used the
identical pair, construct and component map, common-mask SHA-256
`98ddcb111d3f76c73d47f521069fab5412a5c5c75f0dc3035920443da2cec614` — the same
value S3R recorded.

| arm | component-macro AP_bidir |
|---|---:|
| candidate, graph-aware `d=128` | 0.046856 |
| baseline41, mean-pooled `d=41` | 0.035880 |
| trained permuted-label learner | 0.036293 |
| frozen B5 differential | 0.031582 |
| foreign ligand pair | 0.046212 |
| residue-context corruption | 0.027357 |
| ligand-only | 0.025472 |
| zero-`W` chance | 0.025472 |
| within-construct chemistry shuffle (non-gating) | 0.051322 |

| Gate | observed delta | one-sided LCB95 | required | result |
|---|---:|---:|---:|:---:|
| R1 candidate - chance | +0.021384 | +0.016064 | +0.05 | FAIL |
| R2 candidate - B5 | +0.015273 | +0.008173 | +0.03 | FAIL |
| R3 candidate - foreign pair | +0.000644 | -0.009226 | +0.03 | FAIL |
| R3b candidate - ligand-only | +0.021384 | +0.016064 | +0.03 | FAIL |
| R4 candidate - context corruption | +0.019498 | +0.013599 | +0.03 | FAIL |
| R5 candidate - permuted learner | +0.010563 | +0.003880 | +0.05 | FAIL |

| non-gating contrast | delta | LCB95 |
|---|---:|---:|
| C1 candidate - baseline41 | +0.010976 | +0.004939 |
| C2 baseline41 - chance | +0.010408 | +0.006920 |

The earliest failed boundary is R1, so the terminal verdict is
`REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED`.

## Part 3 — what actually changed, and what did not

**The representation change was real and it helped.** The above-chance gain
roughly doubled, from `+0.010408` to `+0.021384`, and the direct contrast
C1 `+0.010976 [LCB +0.004939]` is above zero. S3R lost to its own
capacity-matched permuted-label learner by `-0.001245`; S4R beats it by
`+0.010563 [LCB +0.003880]`. So the graph-aware statistic carries information
that is genuinely linked to the real labels and not reachable by capacity
alone. That is a new fact, and it was not true before this stage.

**It is nevertheless not ligand-conditioned residue selection.** R3 is the
decisive control. Replacing each pair's own two ligands with a frozen foreign
ligand pair, holding the protein, the estimator and `W` fixed, costs only
`+0.000644 [LCB -0.009226]`. The within-construct chemistry shuffle scores
`0.051322`, *above* the candidate. The learned residue direction is therefore
almost independent of which ligand difference is fed to it. Consistently, the
singular spectrum of `W` is nearly flat over its leading directions
(`0.2091, 0.1887, 0.1865, 0.1771, 0.1702, ...`), so `W dg` does not concentrate
on a ligand-specific subspace.

The right reading is that the richer ligand statistic sharpened a
**construct-level** residue-change prior — which residues in this protein are
liable to change at all — rather than a ligand-specific one. That prior is
worth more than the atom marginal could express, which is why C1 and R5 moved,
and it is also why R3 did not.

Two further boundaries were confirmed by construction rather than by
measurement. Residue-context corruption drops the arm to `0.027357`, close to
chance, so the residue-resolved protein context is doing the work. The
ligand-only arm is exactly the chance arm: a residue-constant field lies in
`span{1} subset span{Q_P}` and the gauge annihilates it, so no ligand-only
shortcut can survive this estimand at all. R3b is consequently numerically
identical to R1 and is a structural proof, not an independent empirical
contrast.

## Part 4 — this is not an implementation failure

Module participation passed in full:

```text
gradient non-zero .................. min |grad W| = 3.2617
relative W movement ................ 1.4139   (>= 0.05)
unit Frobenius norm ................ 1.0000000207
held-out raw score variance ........ 1.1635e-04
zero-W arm equals analytic chance .. 0.0254715980 == 0.0254715980
ligand-only equals chance .......... 0.0254715980, pre-projection field 0.03571
context and foreign both degrade ... True
one shared stream across all arms .. True
repeat W / prediction / metric ..... 0.0 / 0.0 / 0.0 max absolute difference
repeat prediction SHA-256 .......... identical to the candidate
```

Three independent cross-stage anchors show that nothing but the ligand
representation moved:

1. the training stream's semantic SHA-256 is
   `4bc68d54884437ded999fbb5f8fc8997b47b456f2bdde4e0fbafd9df3dcdc3ef`, the
   value S3R registered, and the stream file SHA-256 also matches;
2. the common-mask SHA-256 equals S3R's;
3. the `baseline41` arm reproduces the S3R candidate **exactly** —
   `|0.03588006089257408 - 0.03588006089257408| = 0.0` — and C2 reproduces the
   S3R R1 interval to every printed digit.

A failed run whose baseline replicates the prior stage bit-for-bit is strong
evidence that the comparison is clean.

## Part 5 — governance

- Heldout-B was not created, not written and not read. Only `train` and
  `heldoutA` residue views exist in this stage.
- Heldout-A had already been consumed by S3R, so every number here is
  development evidence and none of it is confirmation.
- Label views opened: `train_residue_masks.json.gz` by each of the four
  training processes, `heldoutA_residue_masks.json.gz` by the scoring process.
  No process opened both.
- Real structural residue-edge label reads: `103,116` across the two created
  views (`train` and `heldoutA`). Affinity value reads: `0`. DAVIS, KIBA,
  recipient, ChEMBL and BindingDB reads: `0`.
- No threshold, seed, margin, budget, capacity or representation was changed
  after any S4R metric was read. One run, as registered.

## Part 6 — remaining boundary

Under the registered stopping rules, `REAL_RESIDUE_DIRECTION_STILL_NOT_
IDENTIFIED` closes the pose-free representation repair route. It does not
authorize attention, a larger PLM, a second protein encoder, a parallel branch,
3D pose or geometry, typed interaction channels, affinity supervision,
knowledge graphs, PU learning or few-shot adaptation, and it does not authorize
re-running this stage at `d = 256` or `d = 512`.

What the two stages jointly establish is narrower and more useful than either
alone. Phase 2A showed the labels carry same-construct ligand conditionality.
S2R showed the estimator is trainable. S4R now shows that the ligand
representation was a real, measurable bottleneck and that removing it doubles
the above-chance signal — but that the residual signal is a construct-level
residue-change prior, not ligand-conditioned residue selection, because it
survives foreign ligands intact.

The unresolved question is therefore no longer "is the ligand representation
too poor"; it is whether any pose-free, sequence-plus-2D estimand can bind a
specific ligand substructure to a specific residue context without the
correspondence information that only geometry supplies. Answering that would
require a separately governed information stage, not another repair of this
one.

No biological statistic is admitted to `z`. Affinity, selectivity, few-shot
sectioning, heldout-B, R6, DAVIS/KIBA/recipient labels, CSMO, Band and the
frozen operator `A(F,z)=K(B(z)F(z))` remain untouched.
