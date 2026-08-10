# Research boundary

## Current registered stage

`crossed_interaction/` contains the completed X1A-R path. The amended X1A ICC
PASS is historical. The repaired direct-DD audit used exact-assay-aligned
rectangles and returned `X1A_R_DEPENDENCE_PRECONDITION_FAILED` for both Ki and
Kd. X1B was not run. X2, GPU training, support adaptation and z remain
unauthorized; no active trainable research stage exists.

A label-blind cycle-space feasibility audit is also complete. ChEMBL has large
raw panel quotient dimension but zero exact-assay crossed dimension, and nearly
half of each endpoint's quotient coordinates remain in one dependency
component. It does not reopen training. The only registered research action is
the audit-only BindingDB curated-article source census in
`PREREG_CQ_R0_BINDINGDB_CENSUS.md`.

## Completed structural path

`s7_l2b_r0r/` contains the completed MONN provenance/mapping workflow, frozen
ESM2 B5 exact-residue localizer, Phase 2A attribution audit, synthetic S2R
direct-W witness, real structural S3R transfer, the S4R-A label-blind ligand
representation audit and the S4R single-axis graph-aware repair.

```text
Phase 2A  LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
S2R       BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED
S3R       REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
S4R-A     GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE
S4R       REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED
S5D       LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED
```

S4R is the structural stopping point. It changed exactly one axis of S3R — the
ligand statistic — replacing the mean-pooled 41-D atom marginal with a frozen
radius-1 Morgan per-heavy-atom statistic over a train-only 128-entry
vocabulary. The protein branch, gauge, estimator, loss, sampler, split, seeds,
control maps and 210-update stream were byte-identical, proved by a stream
SHA-256 match, a common-mask SHA-256 match and a bit-exact reproduction of the
S3R candidate by the `baseline41` arm.

The change was real: the above-chance gain doubled to `+0.021384` and the
candidate beat its capacity-matched permuted-label learner. It was not
sufficient: R1 needs `+0.05`, and a foreign ligand pair costs only `+0.000644`,
so the recovered signal is a construct-level residue-change prior rather than
ligand-conditioned residue selection.

S5D then trained nothing and reused those checkpoints to ask why R3 failed. It
registered the collapse mechanism and falsified it: candidate residue fields
have top principal energy fraction `0.4793` against a data-side bound of
`0.4550`, and true-versus-foreign field cosine is `0.4487`. The estimator does
steer on the ligand. Its symmetric-difference conditional estimand, which
cancels pocket membership exactly, then found nothing —
`E1 = +0.011285 [LCB -0.007749]`, `E2 = -0.000440`.

## No eligible structural repair

The registered S4R stopping rule closes the pose-free ligand representation
route. Re-running the stage at `d = 256` or `d = 512`, at radius 2, at another
seed or with another budget is explicitly excluded, as are attention, a larger
PLM, a second protein encoder, a parallel branch, geometry, pose, typed
channels, affinity supervision, knowledge graphs, PU learning and few-shot
adaptation. The S5D stopping rule closes the conditional estimand route and
forbids a fourth estimand variant on heldout-A, consumed three times already.

Ligand information is neither lost upstream nor diluted by the metric; it
arrives, rotates the residue field substantially, and points somewhere
biologically wrong. The remaining hypothesis is that the missing ingredient is
correspondence — which ligand substructure sits against which residue — and
that a pose-free sequence-plus-2D estimand has no channel to supply it. That is
a separately governed information stage about geometry, not a repair of these
stages, and it requires its own preregistration.

## Promotion rule

Research code may enter `model/` or `scripts/` only after its own frozen Gate,
independent structural confirmation, source-affinity increment over ligand-only
and wrong-protein controls, and sealed transfer. No current research statistic
is admitted to production `z`.

Historical failed and superseded stages are evidence in `history.md` and Git,
not current execution instructions.
