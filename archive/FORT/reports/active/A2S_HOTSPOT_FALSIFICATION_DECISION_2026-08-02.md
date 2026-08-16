# A2S-HOTSPOT heterogeneous-sparsity falsification

Date: 2026-08-02  
Branch: `research/a2s-hotspot-sparse-20260802`  
Artifact: `reports/active/a2s_hotspot_falsification_2026-08-02.json`  
Records: `reports/active/a2s_hotspot_falsification_records_2026-08-02.parquet`  
Roles opened: source `fit` and source `probe` only. `locked` and recipient labels were not requested.

**Decision: `STOP_HOTSPOT_SPARSITY_NOT_REPRODUCED`.**

## Binding result

The original 26-dimensional descriptor-plus-Morgan-PCA basis does not contain a stable sparse target
state. Across three deterministic scaffold-disjoint splits of 52 probe targets:

- eight named coordinates retain 63.7% of the full head gain;
- eight coordinates after arbitrary orthogonal rotation retain 76.8% on average;
- the top-three coordinate set has mean Jaccard 0.303 across split seeds (median 0.20);
- effective `s60`, the first coordinate budget retaining 60% of the full gain, is eight;
- matched source-PCA directions retain 24.4%, confirming non-low-rankness but not sparsity.

The rotation is an exact full-model reparameterisation: it preserves every full linear prediction.
Its superior truncated performance means the earlier top-coordinate result is a generic consequence
of truncating a noisy dense ridge head, not evidence for privileged interaction determinants.

## Alternative bases

| fit-only label-free basis | dimension | full head gain | effective `s60` | top-3 split Jaccard | top-8 retained |
|---|---:|---:|---:|---:|---:|
| descriptors | 10 | +0.0410 | 5 | 0.465 | 0.984 |
| Morgan PCA | 26 | +0.0629 | 8 | 0.341 | 0.702 |
| Morgan random projections | 26 | +0.0545 to +0.0674 | 5-8 | 0.297-0.329 | 0.702-0.812 |
| pharmacophore counts + descriptors | 18 | +0.0575 | 8 | 0.443 | 0.740 |

The descriptor basis received one bounded rescue because it was both semantic and more stable. Five
descriptor coordinates retain 72.6%. Five arbitrary rotations of those same descriptors retain
82.4% on average, with no loss of full prediction and comparable or better support stability. The
rescue therefore fails the identical rotation control.

## Consequences

1. The proposed binding-hotspot interpretation is withdrawn. The data support neither biological
   hotspots nor a 2-3 determinant target state.
2. The flat source-head spectrum remains real. It can no longer be reinterpreted as heterogeneous
   coordinate sparsity in any tested basis.
3. Gate G4's k approximately 10 knee cannot be explained by `s log(d/s)` from a real sparse state;
   that numerical agreement was coincidental.
4. Per the user-supplied stop rule, no sparse representation, LISTA, hard-concrete selector, LASSO,
   OMP, or sparse meta-operator will be trained on this branch.

## Successor condition

The sequential evidence instead exposes a measurement-coherence hypothesis. Affinity is a
conditional thermodynamic observable: target construct, assay format and conditions define the
measured state. Pooling assay contexts inflates the residual variance that governs k-shot
identifiability. A metadata-only census found 311 k=5 exact-assay probe groups across 74 targets and
68 components, sufficient for a source-only test.

The successor may proceed only if exact-assay correct support carries assignment-dependent ranking
information in support-query Tanimoto strata below 0.35. This keeps the original passive A2S-DTA
task and does not turn it into active learning, retrieval, selective prediction or calibration.
