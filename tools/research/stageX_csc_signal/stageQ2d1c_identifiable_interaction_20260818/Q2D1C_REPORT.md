# Stage Q2d-1c report - STOP at feature-space oracle precheck (2026-08-19)

Prereg SHA: 25b8b9129120d0a770ba353cd56a8a388dd847778e4e5cb2b488ee0cbfee7106.

## What was fixed vs Q2d-1b, and what was learned

The Q2d-1b oracle failure was traced to a broken ALS implementation (it could
not even fit the training set: in-fit dz 0.50) plus the train-only ID
centring. Q2d-1c switched the oracle to a closed-form truncated-SVD fit on
the complete train grid (exact in-sample, dz/sp 1.0/1.0 in-fit), restricted
the ligand pool to SMILES-resolved ligands, compressed protein features to
PCA-32, and used feature-smoothed double centring (row/column offsets
regressed on features; empirically the fitted offsets are ~0).

## Frozen rule 3.7 fired again - and this time the limit is structural

SVD oracle on M1 (train-only, closed form, diagnostic only):

| surface | seed0 dz/sp | seed1 dz/sp | seed2 dz/sp |
|---|---|---|---|
| protein-cold | 0.715 / 0.398 | 0.621 / 0.299 | 0.798 / 0.582 |
| ligand-cold | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| double-cold | 0.767 / 0.554 | 0.587 / 0.251 | 0.779 / 0.533 |

Seed 1 misses the frozen dz >= 0.70 on protein-cold and double-cold, so the
stage STOPS BEFORE any training; no representation comparison is authorized.
This is a strict FAIL of the frozen gate, not a model result.

## Attribution (tested, deterministic)

- TRUE generative weights reach dz 0.968-0.994 on all three surfaces for all
  three seeds (truth is internally recoverable; the mechanism exists).
- The train-row feature submatrix has rank 28 < 32 (PCA-32 still over-bases
  the train rows), and 8.8% of the true A's energy (seed 1) lies in the
  unidentifiable null space of the train rows. No train-only estimator,
  closed form or gradient-based, can recover that component; it is
  information-theoretically absent from the training cells. The oracle's
  recovered map therefore systematically underestimates cold-row predictions
  (correlation with the true map 0.52, scale 0.30) on the protein-cold and
  double-cold surfaces.
- The ligand side is fully identified (rank 48/48; ligand-cold dz 1.0).

## Consequences

Q2d-1d (new frozen prereg) restricts the truth protein factor map to the
train-row span of the protein features: A_t = V_train^T @ C with V_train the
rank-28 right-singular basis of the train-row feature submatrix and C
(28x4) Gaussian, QR-orthonormal, scales [1.0, 0.8, 0.6, 0.4]. The mechanism
remains feature-conditioned (protein factor = linear function of the pocket
physicochemical features, observable subspace); what changes is that the
truth is now drawn entirely within the information that the training design
can identify. Everything else (ligand features, splits, arms, ladder A-E,
budget, gate 0.30/0.70/0.05, negative-arm rules, value-level reproduction,
oracle precheck before training) is carried over unchanged.
