# Transformer + Bayesian + meta-learning architecture decision (2026-07-25)

## Required architecture

The active model must retain all three components:

1. a Transformer interaction basis over frozen ESM-2 segment tokens and ligand features;
2. an exact Bayesian posterior over a target-specific residual function;
3. episodic meta-learning in which support defines the posterior and query loss trains shared
   parameters across targets.

`model/bayesian_meta.py::TransformerBayesianMetaLearner` implements this contract. Support labels
enter only the centered Gaussian likelihood. They never enter the Transformer or ligand encoder.
The prior is a learned full covariance, target adaptation uses FP32 Cholesky solves, and `k <= 1`
falls back exactly to the ligand base because it contains no ranking contrast.

## First short run

`transformer_bayes_meta_short.json` used five frozen homology-component folds, 600 base steps, 300
inner cross-fit steps and 800 meta steps per fold. It reused the already-spent panel development rows
and is not a confirmatory gate.

| arm | Spearman | RMSE |
|---|---:|---:|
| B0 | 0.2737 | 0.7984 |
| Transformer-Bayes-Meta | 0.2983 | 0.8076 |
| wrong support | 0.2675 | 0.8350 |
| label permuted | 0.2452 | 0.8510 |
| protein shuffled | 0.2929 | 0.8099 |
| protein random | 0.2819 | 0.8122 |

TBM minus B0 was +0.0253 with component-bootstrap interval [-0.0078, 0.0597]. Label specificity
passed (+0.0551, LCB95 +0.0075), while wrong-support and protein specificity did not. Exact `k=0`
fallback, support-order invariance, support-offset invariance and positive finite variance all passed.

## Result-driven adjustment

The two failed specificity controls motivated one localized training change: a weight-0.1 soft
contrast requiring the true task to outperform wrong support and shuffled protein, with meta steps
increased from 800 to 1200. The architecture and evaluation rows were unchanged. This is a post-result
development run, not independent evidence.

`transformer_bayes_meta_adjusted.json` showed that the contrast was harmful: TBM fell to 0.2779,
the gain over B0 fell to +0.0024, RMSE safety failed, and label specificity disappeared. The contrast
did not create protein specificity. It is disabled in the retained configuration.

## Decision

```text
RETAIN_TRANSFORMER_BAYES_META_ARCHITECTURE
RETAIN_PURE_QUERY_META_OBJECTIVE
TRANSFORMER_BAYES_META_SHORT_FAIL_REVIEW
```

The first configuration is retained because it is the only tested model satisfying the required
three-part architecture while improving both ranking and label specificity. It is not yet a passed
design: the effect misses the 0.03 minimum, its LCB95 crosses zero, wrong-support specificity is not
resolved, and protein destruction does not remove the gain significantly.

The next valid improvement must strengthen protein-conditioned interaction information structurally,
not add epochs or a loss that merely makes destruction controls worse. No confirmation, Davis or
sealed data was read. Panel development rows were read and were already spent before these runs.
