# KIRHub CAPIT/ASPIRE candidate-2 P0 decision

Verdict under the frozen gate: **`ASPIRE_P0_FAIL_STOP_ALIGNED_POCKET`**.

This is a near-threshold failure with strong positive mechanism evidence, not a null result.
Across 302 strict full-sequence homology components, component-macro Spearman is 0.4539
[0.4252, 0.4821] for the aligned 85-residue pocket oracle, versus 0.4253
[0.4001, 0.4493] for the KLIFS-group centroid, 0.3983 [0.3712, 0.4238] for pooled ESM,
0.3562 [0.3307, 0.3810] for within-group pocket shuffle, and 0.3670
[0.3435, 0.3903] for random same-group targets.

Paired gains:

- aligned pocket minus group centroid: +0.0286 [+0.0116, +0.0459];
- aligned pocket minus pooled ESM: +0.0556 [+0.0403, +0.0714];
- aligned pocket minus within-group pocket shuffle: +0.0977 [+0.0761, +0.1188];
- aligned pocket minus random same-group targets: +0.0869 [+0.0665, +0.1068].

All information-destruction and fold-support gates passed. The sole failed gate was the frozen
minimum substantive gain over group centroid: +0.0286 is 0.0014 below +0.030. The shared-panel
oracle is deliberately easier than dual-cold prediction because it reads held-ligand measurements
on training targets. Therefore it cannot be upgraded to a predictive claim, and the threshold is
not relaxed after observing the result.

The correct conclusion is that aligned active-site sequence carries real within-group
target-specific ligand-ranking information and is substantially better aligned to this mechanism
than pooled ESM. Evidence is nevertheless insufficient to authorize the CAPIT strict operator
under its preregistration. No BLOSUM tuning, learned residue weighting, fragment tensor, pocket
graph, Transformer, or extra seed is run.

Candidate ledger: 2/3. `sealed_test_consumed=false`; `confirmation_labels_read=false`.
