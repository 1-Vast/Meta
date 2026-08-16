# KIRHub WTPAIR candidate-1 decision

Verdict: **`WTPAIR_C1_FAIL_STOP_FAMILY`**.

The frozen 25 strict homology-component x chemical-component folds retained 308 independent
homology components. Every fold had at least 59 evaluable target profiles and 51 components; each
WTPAIR arm fit 20,000 mixed-difference examples from at least 5,473 eligible same-group,
cross-homology target pairs.

Component-macro Spearman was 0.0323 [0.0112, 0.0538] for WTPAIR, versus 0.0429
[0.0236, 0.0622] ligand-only, 0.0911 [0.0703, 0.1119] KLIFS-group centroid, and 0.0419
[0.0227, 0.0611] for the matched 256-parameter cellwise bilinear ridge. The destruction arms were
0.0099 [-0.0102, 0.0303] for within-group protein shuffle and 0.0302 [0.0093, 0.0512] for matched
random protein.

The paired gains are decisive:

- WTPAIR minus ligand-only: -0.0107 [-0.0326, 0.0107];
- WTPAIR minus group centroid: -0.0588 [-0.0838, -0.0341];
- WTPAIR minus matched cellwise bilinear: -0.0097 [-0.0231, 0.0038];
- WTPAIR minus within-group protein shuffle: +0.0224 [-0.0025, 0.0477];
- WTPAIR minus random protein: +0.0021 [-0.0247, 0.0300].

The required effect was +0.030 (component MDE80 +0.016 at paired SD 0.10). Only fold support and
the permissive RMSE guard passed; every performance and mechanism-specificity gate failed.

Interpretation: directly supervising same-group protein-pair x ligand-pair rank rearrangements does
not make the pooled full-sequence ESM geometry load-bearing under strict dual-cold evaluation.
Failure against the same-capacity cellwise model excludes lack of bilinear parameter count as the
explanation. Failure against random protein excludes a usable continuous ESM direction. The strong
group centroid shows that coarse taxonomy remains predictive, while within-group extrapolation is
not identified by this representation and supervision.

No alpha change, larger rank, Transformer, pocket augmentation, second seed, or WTPAIR mechanism
revision is authorized. Candidate 1 consumes 1/3 slots in the reopened round.

`sealed_test_consumed=false`; `confirmation_labels_read=false`.
