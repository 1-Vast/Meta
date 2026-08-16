# LPCT one-seed architecture decision (2026-07-26)

Preregistration: `reports/active/lpct_preregistration.md`. Result:
`reports/active/lpct_short.json`. Five frozen homology-component folds, seed 1729, CUDA; 600 base,
300 inner cross-fit and 800 episodic meta steps per fold. The already-spent Metz development panel
was used. No confirmation, Davis or sealed label was read.

## Result

| arm | target-macro Spearman | RMSE |
|---|---:|---:|
| B0 | 0.2737 | 0.7984 |
| retained Transformer-Bayes-Meta | 0.2983 | 0.8076 |
| **LPCT** | **0.2892** | **0.8347** |
| LPCT with wrong support | 0.2618 | 0.8357 |
| LPCT with permuted support labels | 0.2801 | 0.8318 |
| LPCT with shuffled protein | 0.2850 | 0.8329 |
| LPCT with random protein | 0.2847 | 0.8340 |

LPCT minus B0 was `+0.0161` with grouped component-bootstrap 95% interval
`[-0.0144, +0.0477]` over 96 paired components. Correct support beat wrong-target support by
`+0.0274 [LCB95 +0.0080]`, but did not beat label-permuted support
(`+0.0103 [LCB95 -0.0168]`). Protein specificity also failed. RMSE exceeded the frozen
`1.02 * B0` safety boundary (`0.8347 > 0.8144`), and negative transfer was 0.538 versus 0.356 for
the retained TBM.

Exact `k<=1` fallback, support-order invariance, support-label-offset invariance and positive finite
variance all passed. The full 146-test repository suite passed before this run.

## Verdict and localization

```text
LPCT_SHORT_FAIL_STOP
FREE_VALUE_NETWORK_BYPASSES_SUPPORT_LABEL_CONTRAST
```

LPCT is a genuine end-to-end learned predictive module, but it is not a successful one on this
substrate. Its learned value path can predict from ligand/protein pair tokens without making the
support affinity contrast load-bearing. The small correct-versus-wrong-support advantage shows that
the support chemical set matters, while failure against label permutation shows that much of the
network's result is support-chemistry conditioning rather than transport of target-specific SAR.

No post-result width, learning-rate, epoch, loss-weight or scheduler rescue is allowed on these spent
rows. A future learned operator must make its prediction algebraically linear or otherwise
non-bypassable in centered support labels, and must be distinguished from existing attentive neural
processes and FS-CAP-style context encoders. Before another architecture is scored, the newly acquired
Reinecke 2024 fine-resolution pKd panel should complete its homology/firewall/power audit so that
further development is not driven by the repeatedly used, rounded Metz panel.
