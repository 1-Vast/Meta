# Stage Q2d-1b report - STOP at feature-space oracle precheck (2026-08-19)

Prereg SHA: 872bc4402f228d940776e7efe2fee6b91e8310badb4e8830f653ca5e5d2e998e.
Forensic tests: 5/5 passed (Q2D1_CORRECTIONS.md).

## Frozen rule 3.7 fired

The feature-space oracle (rank-4 ALS, train-only fit, closed form, diagnostic
only) fails to reach dz >= 0.70 on any cold surface for mechanism M1:

| surface | seed0 | seed1 | seed2 |
|---|---|---|---|
| protein-cold | 0.466 | 0.577 | 0.505 |
| ligand-cold | 0.472 | 0.319 | 0.472 |
| double-cold | 0.563 | 0.343 | 0.444 |

M2 and M3 are equally unidentified (Q2D1B_ORACLE_PRECHECK.json). Per the
frozen rule, the stage STOPS BEFORE any training; no representation
comparison (Q2d-2) is authorized. This is a strict FAIL of the frozen
identifiability gate, not a model result.

## Attribution diagnostic (seed 0, 30-iteration ALS, _diag_centring.py)

| centring variant | pc | lc | dc |
|---|---|---|---|
| train-only ID centring (frozen) | 0.541 / 0.121 | 0.415 / -0.205 | 0.392 / -0.211 |
| global all-cell centring | 0.542 / 0.137 | 0.499 / -0.052 | 0.642 / 0.206 |
| grand-mean only | 0.515 / 0.130 | 0.681 / 0.354 | 0.680 / 0.453 |

(dz / sp per surface.) Two defects are attributed:

1. ID centring injects structure the feature bilinear cannot represent.
   Train-only per-row/per-column means are ID-specific offsets; subtracting
   them pollutes the train targets and biases the factor estimates.
   Grand-mean-only centring recovers dc to 0.68.
2. The protein feature map is unidentified on the row complement. P_t is
   510-dim but the train row span has rank <= 80, so A (510x4) is identified
   only on that span; cold families' rows live partly outside it and the
   protein-cold surface stays at chance (sp ~0.13) under every variant.
   A third defect is pending: unresolved-SMILES ligands carry hash-fallback
   ECFP bits and dominate the cold scaffold clusters.

## Consequences and next stage

The mechanism claim is neither confirmed nor falsified; the harness was
unidentifiable as frozen. Corrected successor Q2d-1c (new frozen prereg):
- protein truth features = global-PCA-compressed pocket physicochemical
  vector (510 -> 32 dims, frozen projection computed before truth
  generation; keeps dims <= train-row rank so the factor map is identifiable
  and cold rows project into the same space);
- ligand pool restricted to SMILES-resolved ligands (no hash-fallback bits);
- feature-smoothed double centring (per-row and per-column means removed
  through feature-linear projections fit on train, so every subtracted
  offset is representable by the learner) with sd over train cells;
- everything else (ladder A-E, 8 arms, 3 seeds, budget, checkpoint rule,
  value-level reproduction, 0.30/0.70/0.05 gate, oracle precheck before
  training) carried over unchanged.
