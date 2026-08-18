# Stage Q2c report (2026-08-18)

Frozen prereg: PREREGISTRATION.md, SHA-256 1027ccde8c8946aa8314ebd7642af89a6abbc3366afd965e8ab43f0da5a26a5c.

## Status

| Item | Status | Artifact |
|---|---|---|
| Q2c-0 harness self-audit | DONE | Q2C0_PROJECTION_AUDIT.json, Q2C0_REPORT.md |
| Oracle alignment | DONE | ORACLE_ALIGNMENT_TABLE.md |
| Q2c-1 representation x learner matrix | DONE | Q2C1_MATRIX.json |
| Q2c-1b oracle tau*=2.0 control | DONE | Q2C1B_ORACLE_TAU2.json |
| Q2c-2 frozen gate rerun (pair-centered local ESM) | NOT STARTED (precondition unmet) | - |
| Q3b pairability audit | DONE (partial) | Q3B_PAIRABILITY_AUDIT.json |

## Q2c-0 headline findings

1. The X0c anova_projection Pearson of 0.51-0.54 on dead heads is explained,
   not by projection-manufactured correlation, but by (a) projecting the full
   fitted output (mu+pm+lm+head) instead of the interaction head, and (b) the
   no_interaction_head arm's inter_scale drifting during training (raw dz of
   the nominally dead head after training: 0.457; it was only zeroed at
   predict time). Interior eval cells = 0/545: eval parents AND eval ligands
   are never seen in train, so the projection operator degenerates to
   grand-mean subtraction on eval. Dead-head + true main effects reproduces
   0.093, close to the random-graph null (-0.064, p=0.169).
2. Projection negative controls: tau*=0 PASS (0.000, p=1.00);
   no_interaction_head PASS (0.021, p=0.637); random graph PASS; ligand_only
   FAIL (0.114, p=0.015) - the projection is therefore removed from all
   diagnostics and gates (frozen rule).
3. Endpoint distortion is material: latent interaction vs observable-implied
   interaction Pearson 0.59 / Spearman 0.555 on determinate eval cells. The
   sigmoid + quantization + censoring pipeline destroys ~40% of the
   interaction's rank information at the observable level.
4. Minimal linear no-censoring Q2: dead-zone sign accuracy 0.52 - linear
   learners cannot read the interaction from one-hot pockets even without
   endpoint distortion.

## Q2c-1 matrix (Q2C1_MATRIX.json)

Median over seeds 0-2, gate point (tau*=1.0, rank 4, dense), dead-zone sign
accuracy (dz) / Spearman (sp) / gap vs ligand_only:

| representation | learner | dz | sp | gap |
|---|---|---|---|---|
| one_hot_pocket | linear | 0.558 | 0.122 | +0.055 |
| one_hot_pocket | mlp (artifact) | 0.504 | 0.033 | -0.022 |
| one_hot_pocket | mlp z-scale | 0.504 | 0.033 | -0.022 |
| pocket_esm (pair-centered local ESM) | linear | 0.478 | 0.034 | -0.035 |
| pocket_esm | mlp | 0.467 | -0.017 | -0.064 |
| pocket_esm | mlp z-scale | 0.508 | -0.017 | -0.029 |
| oracle_PU | linear | 0.532 | 0.041 | -0.011 |
| oracle_PU | mlp (artifact) | 0.664 | 0.375 | +0.128 |
| oracle_PU | mlp z-scale | 0.641 | 0.370 | +0.095 |
| random | linear / mlp | 0.500 / 0.478 | ~0 | ~0 |
| shuffled | linear / mlp | 0.478 / 0.518 | ~0 | ~0 |

Frozen interpretation rules applied:
- oracle-linear pass + E2E fail -> optimization/routing: NOT THE CASE
  (oracle linear 0.532 also fails).
- oracle pass + pocket_esm fail -> representation gap: NOT THE CASE, because
  the oracle itself does not reach dz 0.70 (0.664; sp 0.375 does clear 0.30).
- all probes fail -> harness definition / graph power / truth generation:
  CONSISTENT - every non-oracle cell is at chance (0.47-0.56) and the oracle
  tops out at 0.664.
- Q2c-2 (pair-centered local ESM gate rerun) is therefore NOT STARTED: its
  frozen precondition (oracle passes, ESM fails) is not met, and the
  Q1-passing ESM representation is at chance on the interaction task
  (dz 0.467-0.508).

Material negative result: the Q1-certified representation does not transfer
to the Q2 interaction task; the identity-link sensitivity shows the endpoint
transform is NOT the blocker (z-scale dz within 0.02 of logit-scale).

## Q2c-1b control (Q2C1B_ORACLE_TAU2.json)

oracle_PU at tau*=2.0, rank 4, 3 seeds, official protocol: dz 0.749 / 0.576
/ 0.733, median 0.733 (>= 0.70); Spearman median 0.485. The 0.70 threshold
is therefore REACHABLE by the frozen protocol when the signal-to-noise ratio
is doubled; the failure at tau*=1.0 is a graph-power / SNR limitation, not a
protocol, threshold, endpoint or representation-only artifact.

## Q3b

- S13 carries one column per variant; the paper's duplicate measurement at
  1 uM is not resolved as a second matrix in the committed supplement, so
  duplicate agreement is not computable from S13 alone (data limitation).
- Bemis-Murcko scaffold census limited by name resolution: 14/92 inhibitors
  matched to PubChem SMILES (13 unique scaffolds).
- Effective sample size: WT median activity 95.0 (25.4% exactly 100);
  median 89 variant cells per variant; median 340 cells per inhibitor.
- WT vs variant construct-background equality and per-pair ATP protocol
  remain unresolved dimensions (Q3 census).

## Verdict (final for stage Q2c)

Pipeline qualification: FAILED at Q2 (X0c), failure mode now decomposed by
Q2c-0 / Q2c-1 / Q2c-1b: ANOVA projection was never the source of the
artifact's high dead-arm correlations (full-yhat metric + drifting
negative-control arm); the endpoint transform is not the blocker
(identity-link dz within 0.02); the Q1-passing local-ESM representation does
not transfer to the interaction task (all ESM cells at chance); the oracle
factor input reaches dz 0.664 at tau*=1.0 and dz 0.733 at tau*=2.0, so the
frozen 0.70 dead-zone threshold is UNREACHABLE at the frozen gate point's
SNR on this graph even with oracle factors.

Biological conclusion: UNRESOLVED.
B1/B2/C/D: NOT AUTHORIZED.
Single highest-information next step: a new preregistration (Q2d) that
either moves the gate point along the empirical power curve (e.g. tau*=2.0,
where the oracle clears 0.70 and negative controls can be re-verified) or
enlarges the synthetic graph to the sample size the power curve requires for
tau*=1.0; the existing frozen gate may not be moved retroactively.
