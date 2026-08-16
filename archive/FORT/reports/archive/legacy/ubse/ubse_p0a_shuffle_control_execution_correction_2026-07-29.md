# UBSE-P0A shuffled-sequence control execution correction

Date: 2026-07-29  
Frozen at: 2026-07-29T19:21:39+08:00  
Status: frozen before any seed weight, final result, validation ledger, or
performance metric was available

## Trigger

An independent code audit found that the running implementation calls
`evaluate_controls(..., seed + epoch)`. The preregistration specifies one
deterministic within-sequence permutation for each training seed. The running
code therefore changes the shuffled-sequence realization across reported
epochs 1, 2, and 4.

The source file loaded by the running process had SHA-256
`c1732e6d4cfacd8e9ee42ab03bde20d14d98a5ea9d9dff0c6f247cf998211ad7`.
At correction freeze, `ubse_p0a_seed1729_1731.json`, the validation ledger,
and all `ubse_p0a_seed*.pt` files were absent.

This deviation does not change training batches, gradients, learned weights,
correct-sequence predictions, propensity predictions, or constant-position
predictions. It changes only the shuffled-sequence control realization.

## Binding correction

1. Let the running three-seed training finish without changing its optimizer,
   epochs, batches, thresholds, or weights.
2. Preserve hashes of the raw result and raw validation ledger before any
   correction.
3. Load each final saved seed weight and re-evaluate the same 64 validation
   targets with exactly one within-sequence permutation keyed by the training
   seed: 1729, 1730, or 1731.
4. Construct a corrected accepted ledger by replacing only the epoch-4
   `shuffled_sequence` rows. Add an explicit `control_seed` field. All other
   epoch-4 rows must be byte-equivalent in their metric values to the raw
   ledger.
5. Recompute the complete P0A gate and summary from the corrected accepted
   ledger. Do not hand-edit the decision or select a favorable permutation.
6. Retain epoch-1 and epoch-2 raw shuffled rows only as diagnostic trajectory
   records with their actual seeds (`training_seed + epoch`). They are not
   used in any gate or cross-epoch shuffled-control claim.
7. Verify exactly 64 targets, four controls, and three seeds at epoch 4;
   require no duplicate `(seed, epoch, target_key, control)` rows and no audit
   target or contact label.

The accepted P0A-4 decision must use the corrected fixed-seed epoch-4
shuffled control. A disagreement in any unaffected metric, missing weight,
or inability to reconstruct the frozen validation substrate invalidates the
correction and stops P0A rather than authorizing a retraining rescue.

## Claim boundary

This correction cannot change any threshold, model, dataset, split, or
scientific claim. A P0A pass still freezes only a target-marginal pocket
proposal for a later interaction-source experiment; it cannot establish
pair-specific semantics or unlock affinity, Stage-2, confirmation, or sealed
access.
