# R-MAON G0 decision

**Overall verdict:** `RMAON_G0_TOPOLOGY_OR_POWER_STOP`
**G0-A topology:** `RMAON_G0_TOPOLOGY_OR_POWER_STOP`
**G0-B estimator:** `RMAON_G0_NULL_SCORE_AND_RECOVERY_PASS__MODULE_ONLY`

## Frozen results

| regime | replicates | rejection / power | median recovered / truth |
| --- | ---: | ---: | ---: |
| `S1_null` | 200 | 0.055 | -- |
| `S1_active` | 200 | 1.000 | 1.0027 |
| `S2_heteroscedastic` | 200 | 1.000 | 0.9993 |
| `S3_degenerate_signal` | 200 | 0.995 | 0.9610 |
| `S5_null_coordinate` | 200 | 0.045 | -- |

## G0-B gates

- `G0B_NULL`: `True`
- `G0B_WRONG_COORDINATE`: `True`
- `G0B_POWER_S1`: `True`
- `G0B_POWER_S2`: `True`
- `G0B_GRACEFUL_S3`: `True`
- `G0B_RECOVERY`: `True`

## Boundary and numerical audit

- The affinity-bearing Parquet scan was reader-filtered to TRAIN:
  `True`; rows: `12574`;
  components: `101`.
- Empirical covariance ranks were 9--64
  (median 64.0); 30 target covariances
  were rank deficient. The PSD factor added jitter: `False`.
- Component-balanced whitening retained 11 of
  14 coordinate dimensions; weighted Gram maximum error was
  1.688e-14.
- The prospective manifest existed: `False`.

## Result boundary

G0-B tests the regular-null direct-operator estimator on synthetic outcomes under empirical Metz TRAIN
noise. It is not a real-affinity performance result. G0-A is independently binding: without a frozen
powered prospective multi-family manifest, A0, M1, strict dual-cold training, and prediction remain
blocked. The assay-monotone multi-fidelity innovation was not tested.

No development, Davis, confirmation, or sealed label was read. `sealed_test_consumed=false`.
