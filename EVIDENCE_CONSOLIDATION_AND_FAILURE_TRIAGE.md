# Evidence consolidation and failure triage

## Earliest unresolved boundary

```text
sequence + ligand graph
    -> correct-protein contact/distance geometry       PASS
    -> partner-sensitive structural compatibility     PASS
    -> aggregate mechanism readout                     LIGAND-DOMINATED
    -> actual pair-local P1B mechanism observability   NOT TESTED
    -> correct-protein affinity increment              NOT IDENTIFIED
    -> k<=5 identifiable support section               NOT IDENTIFIED
    -> biological z admission                          NOT AUTHORIZED
```

## What is not broken

- frozen mathematical operator contracts;
- canonical data and label firewall;
- open holo acquisition/governance;
- P1B correct-protein contact/distance geometry;
- ChEMBL37 release provenance and closure-safe folds.

## What failed historically

- index-wise and compressed MIF readouts lacked stable affinity increment;
- nonlinear global MIF probes did not rescue the signal;
- task-local radial headroom was primarily ligand/series SAR;
- consumed kinase panels did not identify a `k<=5` double-held-out section;
- external replication failed;
- broad public panels had assay noise exceeding the estimated interaction;
- pose-free aggregate typed features and S4 aggregate pseudo-labels were not
  protein-specific.

The executable implementations were removed after their findings were recorded
in `history.md`.  Their exact trees remain recoverable from commits `3281780`,
`12a2765`, and `608decf`.

## Correct interpretation of S4

S4 tested mean-pooled ESM + ECFP with Ridge.  It did not test P1B's atom-local,
residue-local and pair-local states and therefore cannot close the sequence+2D
model class.  Its valid conclusion is limited to ligand-dominated aggregate
features.

## S5 decision tree

```text
mapping/chain contract fails
  -> repair data contract

oracle teacher or slot ceiling fails
  -> repair teacher/mapping/slot representation

synthetic head fails
  -> repair objective/optimization

oracle passes, actual frozen P1B fails
  -> pose-free P1B inputs insufficient for named channel
  -> separately register pose-aware stage

actual frozen P1B passes structural controls
  -> freeze channels
  -> separately register real source-affinity Gate
```

## Production admission

No structural statistic enters `model/` or `z` merely because it reconstructs a
3D pseudo-label.  It must also show source closure-OOF affinity increment over
ligand-only and wrong protein, then sealed transfer.  Few-shot adaptation must
be restricted to the support design row space, with rank, conditioning, query
coverage and abstention reported.  The frozen theory does not automatically
certify pairwise ranking claims.
