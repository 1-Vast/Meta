# Preregistration - P1R2B-PHASE2B-S3R

## Gauge-free real structural ligand-conditioned residue residual

Stage identifier: `P1R2B-PHASE2B-S3R_REAL_STRUCTURAL_DIRECT_W`

Written on 2026-08-10 before any S3R implementation or S3R biological metric.
This is a single development run on the already governed MONN structural
labels. It reads no affinity value.

## 1. Authorization and question

Phase 2A established
`LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING`. The sealed S2R
synthetic control then established `BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED`
for a bounded direct matrix estimator. S3R asks one question:

> Do frozen ESM2 residue states and the frozen 41-D ligand chemistry summary
> support a ligand-conditioned residue ranking beyond the generic B5 pocket
> prior on unseen protein closure components and unseen ligand graphs?

The historical factorized Phase 2B run remains
`PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED`; it is not
retroactively rescued. S3R is a new estimator contract.

## 2. Frozen biological inputs

S3R reuses the Phase 2B R1 data and evaluation contract byte-for-byte:

- frozen ESM2-650M residue states `H_P in R^(L x 1280)`;
- frozen permutation-invariant ligand summary `g(L) in R^41`, the arithmetic
  mean of the existing deterministic atom features;
- frozen protein-only B5 prior `b^P(P)` and nuisance basis
  `Q_P = orth([1, b^P(P)])`;
- exact sequence/construct identity, graph identity, non-empty distinct Murcko
  scaffolds, protein closure components and held-out A/B splits;
- MONN binary residue masks, used only in this real structural stage;
- the already materialized foreign-ligand and within-construct derangement
  controls;
- the sealed B5 pair predictions and all R1-R6 metric definitions.

No PLM, atom encoder, geometry branch, attention block, typed interaction,
affinity head, PU loss, knowledge graph or parallel module is added or trained.

## 3. The only trainable object

One direct matrix, no bias:

```text
W in R^(1280 x 41), 52,480 trainable parameters
delta_raw(P,L) = (I - Q_P Q_P^T) H_P W g(L)
d_ab            = delta_raw(P,La) - delta_raw(P,Lb)
d_dir           = d_ab / sqrt(mean_r(d_ab^2) + 1e-12)
```

The projection is accumulated in float64. Positive rescaling leaves the AP and
`d_dir` unchanged. After every optimizer update:

```text
W <- W / ||W||F
```

The candidate evaluated by R1-R6 is the full direct `W`. A rank-8 SVD
truncation is reported only as a non-gating compression diagnostic and cannot
change the verdict or select a model.

This is a trainable deep bioinformatics model: ESM2 supplies frozen deep
residue representations, while `W` is the learned ligand-conditioned biological
measurement. It is not yet a production model or a biological `z`.

## 4. Eligible real pairs and labels

An eligible unordered pair `{La,Lb}` has one exact construct, distinct graph
keys, distinct non-empty scaffolds and a non-empty residue-mask symmetric
difference. For masks `Ra,Rb`:

```text
G_ab = Ra \ Rb
L_ab = Rb \ Ra
```

The reverse order is not a second observation. Training and inference preserve
the existing hierarchy:

```text
residue -> unordered ligand pair -> construct -> closure component -> macro
```

Residues and ligand pairs are never inference units.

## 5. Objective and optimization

For each pair, use the all-residue bidirectional ordinal loss on `d_dir`:

```text
L_pair = 0.5 * [
  mean_{g in G_ab, j notin G_ab} softplus(-(d_g-d_j))
  +
  mean_{l in L_ab, j notin L_ab} softplus(-((-d_l)-(-d_j)))
]
```

The objective is balanced pair -> construct -> component -> batch. The frozen
sampler visits every training component per epoch, samples at most two
constructs per component and eight pairs per construct, and batches 16
components. The sampled identifiers are materialized and hashed.

```text
optimizer: Adam
learning rate: 1e-3
weight decay: 0
gradient clipping: 5.0
epochs: 6
expected updates: 210 (must be verified exactly)
parameter seed: 20260901
sampler seed: 20260902 + epoch
model selection: none
early stopping: forbidden
device: CUDA in the `drug` environment when available
```

The identical estimator, sampler and budget train the R5 permuted-label arm.
One same-seed repeat must reproduce the checkpoint tensor and held-out
predictions within `1e-7` max absolute error; exact byte identity is reported
separately and is not required across CUDA serialization.

## 6. Frozen controls

The following controls are inherited without re-selection:

- `R2`: frozen B5 residue differential;
- `R3`: both ligands replaced by the materialized scaffold-distinct foreign
  training pair;
- `R4`: contextual ESM2 states shuffled within amino-acid type in the same
  protein, with `b^P`, `Q_P` and the candidate recomputed consistently;
- `R5`: identical direct-W learner trained on the frozen within-construct
  ligand-label derangement;
- chemistry-shuffle and wrong-protein arms remain secondary diagnostics;
- `R6`: `G_2B(P,L) = G_B5(P,L) + delta_raw(P,L)` with the same atom-residue
  mask and tie-aware estimator used by Phase 2A.

Wrong ligands are nuisance controls, not biological non-binders.

## 7. Metrics and Gates

For every eligible pair over all aligned residues:

```text
AP_gain:   score d_ab,  positives Ra \ Rb
AP_loss:   score -d_ab, positives Rb \ Ra
AP_bidir = mean(AP_gain, AP_loss)
```

Aggregation and inference use closure components. One-sided 95% lower bounds
use 10,000 paired component bootstrap resamples with seed 20260903. The original
R1-R6 thresholds are unchanged:

| Gate | Contrast | Margin | Interval requirement |
|---|---|---:|---|
| R1 | candidate - exact per-pair chance | >= 0.05 | LCB95 > 0 |
| R2 | candidate - frozen B5 differential | >= 0.03 | LCB95 > 0 |
| R3 | candidate - foreign-ligand pair | >= 0.03 | LCB95 > 0 |
| R4 | candidate - residue-context corruption | >= 0.03 | LCB95 > 0 |
| R5 | candidate - trained permuted-label learner | >= 0.05 | LCB95 > 0 |
| R6 | `G_2B` pair AP - sealed B5 pair AP | >= -0.005 | LCB95 >= -0.005 |

Held-out B is scaffold-strict and secondary. It must not show a sign reversal,
but cannot rescue a failed primary Gate.

## 8. Fail-closed preconditions and module participation

Before metrics are opened, S3R must verify:

1. preregistration, S2R verdict, inputs, predictions, controls and code hashes;
2. train/held-out component overlap and train/held-out ligand graph overlap are
   zero;
3. the real pair census exactly matches the frozen R1 census;
4. float64 projection error <= 1e-12, projection orthogonality <= 1e-8,
   antisymmetry <= 1e-10 and positive-scale normalization <= 1e-8;
5. `||W||F = 1` within 1e-5 after every update;
6. S2R calibration and sealed artifacts all contain PASS and their hashes match;
7. no affinity-marked source is opened.

Before R1-R6 determine a verdict, module participation requires:

- nonzero W gradient in every epoch;
- relative W movement from initialization >= 0.05;
- non-degenerate candidate score variance >= 1e-8;
- zero W produces chance AP;
- context corruption and ligand substitution both reduce candidate AP;
- the same-seed repeat reproduces predictions within `1e-7`;
- R5 uses the identical trainable estimator and update stream;
- the final checkpoint, sampled stream and complete held-out prediction table
  are materialized and hashed.

Only output-level `delta_raw` and residue rankings may be interpreted. Entries
of W and rank-8 singular vectors are not physical energy channels.

## 9. Earliest-failure verdict

Exactly one verdict is written:

```text
S3R_CONTRACT_OR_ARTIFACT_FAIL_CLOSED
S3R_NUMERICAL_OR_REPLAY_PRECONDITION_FAILED
PHASE2B_MINIMAL_RESIDUAL_NOT_IDENTIFIED
PHASE2B_SHORTCUT_DEPENDENCE
PHASE2B_RESIDUE_DIFFERENTIAL_IDENTIFIED_BUT_B5_INTEGRATION_FAILED
STRUCTURAL_LIGAND_CONDITIONED_RESIDUE_STATISTIC_IDENTIFIED_IN_DEVELOPMENT
```

Rules:

1. Contract, hash, firewall or census failure stops before training.
2. Numerical, norm, replay or module-artifact failure stops before biological
   interpretation.
3. R1 or R2 failure gives `PHASE2B_MINIMAL_RESIDUAL_NOT_IDENTIFIED`.
4. R3, R4, R5 or module participation failure gives
   `PHASE2B_SHORTCUT_DEPENDENCE`.
5. R1-R5 PASS and R6 FAIL gives the integration verdict.
6. All conditions PASS gives development evidence only and authorizes one
   separately preregistered sealed structural confirmation.

## 10. Frozen downstream boundary

S3R does not read affinity and cannot identify binding energy, selectivity,
off-target activity, k-shot adaptation or a production biological statistic.
It does not modify:

```text
A(F,z) = K(B(z)F(z))
```

Even a complete PASS admits no raw residue map to `z`. Independent structural
confirmation, a source-affinity correct-vs-ligand and correct-vs-wrong-protein
Gate, support-rank/coverage analysis and a separate finite-dimensional `z`
bridge remain mandatory in that order.
