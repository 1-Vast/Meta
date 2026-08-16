# K-LBP v2 R3 decision

**Verdict:** `R3_ESTIMATOR_INSENSITIVE_NO_DECISION`

## Gates

- `G_R3_0_execution_valid`: `False`
- `G_R3_1_recovery`: `False`
- `G_R3_2_null_control`: `True`
- `G_R3_3_N_biased_upward`: `True`
- `G_R3_4_graceful_degradation`: `True`

## Claim boundary

This stage certifies estimator behavior on synthetic coefficients under empirical TRAIN-only noise.
It is not affinity evidence and does not authorize R5-R7. A non-200 run is a smoke run only and can
never certify the estimator.
