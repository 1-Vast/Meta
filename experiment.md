# Active experiment protocol

## Status

```text
S7/L2B Phase 1 B5 ................. COMPLETE, DEVELOPMENT PASS 6/6
Phase 2A teacher audit ............ COMPLETE
S2R synthetic direct-W witness .... COMPLETE, PASS
S3R real structural direct-W ...... COMPLETE, FAIL AT R1
heldout-B / R6 / affinity ......... NOT OPENED
active training stage ............. NONE
```

## S2R

The gauge-free direct matrix estimator passed three calibration seeds and one
sealed synthetic seed. Sealed held-out component-macro `AP_bidir = 0.6620`.
This established trainability only; it did not establish biology.

## S3R

Frozen inputs: ESM2 residue states and mean-pooled 41-D ligand atom features.
Trainable object: one `1280 x 41` direct matrix, unit Frobenius norm. Training:
210 fixed Adam updates, no weight decay, hierarchical pair/construct/component
aggregation, no hyperparameter or seed selection.

The primary 46,818-pair, 112-component panel returned:

```text
candidate AP                  0.035880
chance AP                     0.025472
R1 candidate - chance        +0.010408  [LCB +0.006920]  FAIL (< +0.05)
R2 candidate - B5            +0.004298  [LCB -0.001630]  FAIL (< +0.03)
R3 candidate - foreign       +0.000145  [LCB -0.002651]  FAIL (< +0.03)
R4 candidate - context       +0.003544  [LCB -0.003156]  FAIL (< +0.03)
R5 candidate - permuted      -0.001245  [LCB -0.006595]  FAIL (< +0.05)
```

Terminal verdict: `REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED`.

Module participation and deterministic replay passed. Heldout-B and R6 were not
opened. Affinity reads were zero. No production code or frozen mathematics was
modified.
