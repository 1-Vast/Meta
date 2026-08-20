# S1 ADDENDUM ADD-4 — device switch to CUDA + vectorized training (frozen 2026-08-20)

Frozen BEFORE the GPU run, per PREREGISTRATION B.1 ("a dated addendum records
otherwise BEFORE the run"). Triggered by user review of host utilization: the
whole stage had been running on CPU (my own frozen B.1 choice for seed
determinism), and the training loops were Python-overhead-bound (per-batch
list comprehensions + torch.stack), leaving the GPU idle.

## Frozen changes

- Device: cuda if torch.cuda.is_available() else cpu (all tensors: model,
  ligand features, precomputed arm feature matrices, Form-1 per-pair protein
  tensors). RTX 4060 Laptop GPU verified available 2026-08-20.
- Performance fix (semantics unchanged): per-arm feature matrices and target
  vectors are precomputed ONCE per (arm, estimand); epoch loops become pure
  index gathers over tensors. Form-1 precomputes per-pair protein/ligand
  tensors once. No hyperparameter, architecture, rng-stream, estimand, or
  metric changes: LR 1e-3, WD 1e-4, 200 epochs, batch 512, clip 10, keyed
  streams S1.order/S1.boot/... unchanged.
- Checkpoint selection stays best-val-MSE per estimand; the val criterion is
  computed exactly as before (per-pair finite-mask MSE, pair-summed).

## Determinism and disposition of prior runs

- Seeds remain keyed (stable_rng + torch.manual_seed). cuBLAS matmul is
  deterministic on this machine for identical shapes/weights; run-to-run
  reproducibility on this machine is expected. Cross-device (CPU vs GPU)
  bit-exact equivalence is NOT guaranteed and is not required by any frozen
  rule.
- All prior CPU runs (SMOKE, SEED1_RESULT_BROKEN_T2.json,
  SEED1_RESULT_DIAG_T3BUG.json, SEED1_RESULT_ADD2_T3SIGNBUG.json) remain
  archived diagnostic artifacts; the final adjudicated cells come from the
  GPU runs of the final code state. T1/T0m/T2 cells from the ADD-3 CPU run
  are expected to match the GPU run within float32 tolerance (same code
  paths, same hyperparameters).

## Verification gate before the full run

- tests/test_s1_structure.py 17/17 PASS on the vectorized code.
- CPU smoke (2 pairs, 5 epochs) PASS on GPU.
- F8f floor (0/9 nonconstant), C-perm ~0, and permutation-drop behavior must
  reproduce in the final seed-1 run before adjudication (spot-checked against
  the ADD-3 CPU diagnostic run).
