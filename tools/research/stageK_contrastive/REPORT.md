# Stage K / K2 report — contrastive coembedding (K-REG configuration)

Development evidence; meta_val read once per seed after freezing; meta_test
sealed and never opened. Authorities: KREG_vs_T2.contrast.json,
K2_multiseed_contrast.json, per-seed row summaries, PREREGISTRATION.md,
PREREGISTRATION_K2.md.

## Screen (single seed, 20260815)

Two preregistered training objectives on the coembedding branch: K
(episodic InfoNCE) FAILS (k=0 MSE 2.6242 vs T2 2.5961; mixed). K-REG
(positive/negative regression alignment) PASSES the screen gates with
resolved single-seed gains: k=2 MSE -0.0744, k=3 -0.0739, k=5 -0.0465,
k=0 centered -0.0601, k=0 Pearson +0.1025; no resolved degradation;
controls clean. K-REG was therefore selected for the three-seed
confirmation.

## Confirmation (3 seeds: 20260815/16/17, pooled component bootstrap)

K-REG minus T2, pooled across seeds:

| k | MSE | centered | Spearman | CI |
|---|---|---|---|---|
| 0 | -0.1118 [-0.1851, -0.0490] R | -0.0154 [-0.0304, +0.0001] | +0.0052 | -0.0038 |
| 1 | -0.0480 [-0.0841, -0.0147] R | -0.0154 | +0.0052 | -0.0038 |
| 2 | -0.0273 [-0.0465, -0.0095] R | -0.0114 | +0.0101 | +0.0002 |
| 3 | -0.0218 [-0.0363, -0.0075] R | -0.0100 | +0.0122 | +0.0017 |
| 5 | -0.0122 [-0.0222, -0.0017] R | -0.0068 | -0.0004 | -0.0045 |

Per-seed k=0 MSE: 2.4592 / 2.8012 / 2.4395 (T2: 2.5961 / 2.9811 / 2.4581);
k=5 MSE: 0.9394 / 0.9924 / 0.9704 (T2: 0.9859 / 1.0072 / 0.9458). No
support/protein control inversion in any seed.

## Gates and verdict

K2-1 PASS (k=2/3/5 MSE pooled resolved); K2-2 FAIL (the k=0 centered gain
did not survive pooling: -0.0154 [-0.0304, +0.0001]); K2-3 PASS (no pooled
resolved degradation in any metric); K2-4 PASS (no inversion in any seed).

**Per the preregistered stop rule, the configuration is NOT CONFIRMED and
nothing is promoted: no meta_test is opened, no code moves to model/ or
scripts/.** The record is nevertheless the strongest in the project: the
first mechanism whose ALL-k MSE improvements resolve across three seeds
with ranking preserved and honest controls.

## Mechanism reading (for the ledger)

The coembedding branch with regression alignment reduces the within-target
representation collapse (pairwise cosine 0.99859 -> 0.9908), and the gains
are concentrated in the level/calibration direction (k=0 level^2 pooled
unresolved but MSE resolved) rather than in shape: the centered term barely
moved across seeds. The mechanism is therefore a mild training-time
representation regularization that improves calibration consistency, not a
new information source - consistent with every measurement in D0/D0b/D0c
that the level is assay-history-limited. The k=0 MSE stays >= 2.44 in
every seed: the <=1.00 k=0 target remains out of reach, and the bounded
conclusion in report/BOUNDARY_20260817_NIGHT.md stands.
