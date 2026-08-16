# DCST-R14 label-blind transport-support decision

Date: 2026-07-29  
Decision: `STOP_ABSOLUTE_PLINDER_TRANSFER_TEST_SISMT_AND_DTIOD_GATES`

## Result

R14 reproduced all 2,106 exact-target, firewalled R6 source-train rows and
loaded no ChEMBL affinity column. The frozen privileged teacher remained
pair-responsive on ChEMBL train: its median anchor-centered moment RMS was
`0.141195`, versus `0.149125` in PLINDER, a ratio of `0.946823`.

Responsiveness did not imply transport support:

- only `5.0089%` of ChEMBL-train targets had maximum source sequence 4-mer
  containment at least `0.40`;
- only `9.8500%` of 20,000 audited ChEMBL-train ligands had maximum source
  Morgan Tanimoto at least `0.40`;
- source-versus-train domain AUC was `0.908641` for target features,
  `0.882452` for ligand features, and `0.941130` for the privileged
  anchor-centered moment.

The NoPriv centered-moment domain AUC was lower (`0.778072`) but its
ChEMBL-train centered RMS (`0.106305`) was not negligible. Thus pair
variation alone does not identify privileged transport.

The formal R14 route was `STOP_PLINDER_SOURCE_EXPAND_STAGE1`: both entity
overlap gates failed, while pair responsiveness passed. Runtime was
`123.454 s` and peak allocated CUDA memory was `981.8 MiB`. Confirmation
features, confirmation labels, and sealed test were not consumed.

## Scope of the stop

This stops transferring the complete or raw R6 absolute representation from
PLINDER to ChEMBL. It also rules out global MMD/OT forcing, another nonlinear
fusion head, or importance weighting over the full source as a justified
rescue.

It does not yet rule out the two user-requested information-object changes:

1. SISMT may retain a small but powered intersection of target-supported and
   privileged-certified spectral directions. It must first report retained
   dimension, stability, transported coverage, and effective source sample
   size without affinity.
2. DTIOD may replace absolute states with local mixed finite-difference
   responses. It must first show privileged tangent semantics, held-mask
   student recovery, and improved target support without affinity.

Either route stops before Stage 2 if its independent label-blind gate fails.

