# K-LBP v2 R3 decision

**Verdict:** `R3_ESTIMATOR_INSENSITIVE_NO_DECISION`

## Gates

- `G_R3_0_execution_valid`: `False`
- `G_R3_1_recovery`: `False`
- `G_R3_2_null_control`: `False`
- `G_R3_3_N_biased_upward`: `False`
- `G_R3_4_graceful_degradation`: `False`

## Execution Failure

The registered run stopped after 1 replicate in `S1_gamma_0.0` because
`alternating_gls_did_not_converge_within_frozen_200_iterations`. The five Variant-E fold fits
converged 3/5; N and N2 converged 1/1 and 1/1.

Not run after the execution gate failed: `S1_gamma_0.5`, `S1_gamma_1.0`,
`S2_heteroscedastic`, `S3_degenerate_signal`, `S4_rank2_misspecification`,
`S5_null_coordinate`.

## Claim boundary

The execution failure means estimator sensitivity was not certified; G-R3-1 through G-R3-4 were not
interpreted. A non-200 run is a smoke run only and can never certify the estimator. This is not
affinity evidence and authorizes neither R4 nor R5-R7.
