# FACTOR-C F0-C1 decision

Date: 2026-07-26  
Preregistration SHA-256:
`10DBA441B5193891BDE898E05A45E46C3090C8095E14A402FB23111E14E6B9F1`

## Frozen verdict

`FACTOR_F0C1_REPRESENTATION_UNIDENTIFIED_STOP`

F1-C remains locked. This is not a real-chemical-support failure and is not a verdict that FACTOR as
a whole is false.

## What passed

- All losses and gradients were finite in all three folds.
- Masked structural prediction was strong. Directed bond accuracy was 1.000 in every fold;
  attachment-bin margins were +0.347 to +0.390; masked-element margins were +0.150 to +0.170.
  The Reinecke-held element margin was +0.14969, only 0.00031 below its frozen +0.15 gate.
- Frozen atom embeddings encoded pharmacophore role strongly: 5-NN role macro-F1 was 0.753 to
  0.804, with margins of +0.724 to +0.776 over the frequency baseline.
- Every chemistry-corruption and role-derangement carrier changed. Environment/motif mismatch
  changed 99.45% to 99.66% of carriers, above the 95% gate.
- True coverage exceeded mean broken-chemistry decoy coverage by +0.5091, with source-stratified
  LCB95 +0.5004.
- Corpus firewalls were zero on connectivity and scaffold in every fold; each source had exactly
  one-third total external weight. The current run did not open the quarantined ChEMBL confirmation
  partition and did not consume the sealed test.

## What failed

The learned representation did not meet the preregistered identifiability conditions:

1. Atom-embedding participation rank was only 8.15 to 10.05, below 16; rank ratios were 0.064 to
   0.078, below 0.10. The encoder reconstructs local chemistry well but concentrates it in too few
   directions for the frozen distance atlas.
2. Atom-level broken-chemistry false coverage was 0.0660 to 0.0748, above 0.05. Pair and motif
   calibration passed in all folds, so the remaining calibration defect is specifically atom-level.
3. Inner scaffold-OOD medians were 0.7735, 0.8093 and 0.8203, all below 0.85. Inner q10 passed for
   KIRHub-held (0.6131) and Papyrus-held (0.6211) but failed for Reinecke-held (0.5732).
4. Source-balanced external coverage was median 0.6710 and q10 0.4789, below 0.90/0.70. Source-level
   medians were 0.6946 KIRHub, 0.6393 Reinecke and 0.6879 Papyrus-Christmann.

Because rank and atom-decoy calibration failed before the external gate can be interpreted, the low
external coverage cannot be classified as genuine absence of transferable chemical primitives.
Increasing epochs, width, seeds or changing thresholds is not authorized.

## Resource interpretation

The 96-pair cap is a label-blind audit weighting bound, not a claimed full-training limit. The local
RTX 4060 and host memory completed the fixed run. Since the scientific gates did not pass, an
uncapped high-memory replication is not yet warranted; if a future carrier mechanism passes, that
replication must occur on the larger machine before full DTA training.

## Integrity note

The first formal run was stopped before any fold result because SciPy exact 5-NN degraded severely
in 128 dimensions. It was restarted from scratch using chunked CUDA `cdist`, preserving exact
Euclidean distance, k, data, model, seed and gates. No partial checkpoint or result was reused.

The existing ChEMBL confirmation partition remains permanently quarantined:
`project_historical_confirmation_labels_read=true`. For this run,
`current_run_confirmation_labels_read=false`; `sealed_test_consumed=false`.
