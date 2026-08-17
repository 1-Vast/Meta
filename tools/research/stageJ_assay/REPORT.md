# Stage J report — assay-aware level head + paired level alignment: REJECTED

Development evidence, single seed, meta_val read once after freezing;
meta_test sealed. Authorities: J_vs_T2.contrast.json, per-arm row summaries,
PREREGISTRATION.md, D0c_JOURNAL_IDENTIFIABILITY.json, per-arm RESULT.json.

## Verdict

**Rejected by gates G2/G3; nothing promoted.** G3 fails with RESOLVED
ranking degradation: J minus T2 k=2 Spearman -0.0624 [-0.1195, -0.0044],
k=2 CI -0.0282 [-0.0578, -0.0022], k=3 Spearman -0.0598 [-0.1153, -0.0048].
G2 fails: k=2/3/5 MSE all worse (unresolved). Stop rule S2 fires.

The level head itself works as measured: k=0 level^2 1.7314 -> 1.3207
(-0.4107, unresolved; MSE -0.3941 [-0.8918, +0.0312]) - the largest k=0
level improvement ever recorded in this project - but the candidate trades
ordering for calibration, exactly the trade the preregistration forbids.

## Attribution

- Journal covariates are real: D0c probe 1.619 vs 2.155 constant (shuffle
  2.522); the trained head with journal features reaches level^2 1.297-1.321
  vs 1.502 without them (J-NOJRNL), at k=0.
- The paired alignment term (I2j) adds nothing measurable at this budget
  (J vs J-NOPAIR within 0.004 at k=0; both 1.30-1.32 level^2).
- The k>=1 substitution cost persists despite the residual-role supervision:
  at k>=1, TRUE support labels calibrate level near-optimally, and any
  learned zero-shot level estimate the head supplies replaces part of that
  calibration while also shrinking the Tanimoto transport's shape residuals
  - the second time this interaction has degraded a candidate (Stage E).

## Numbers (frozen meta_val banks, component-weighted, restored pK^2)

| arm | k=0 MSE | k=0 level^2 | k=1 MSE | k=2 MSE | k=3 MSE | k=5 MSE |
|---|---|---|---|---|---|---|
| T2 (frozen) | 2.5961 | 1.7314 | 1.7712 | 1.3245 | 1.2197 | 0.9859 |
| J-NOJRNL | 2.3686 | 1.5016 | 1.6753 | 1.2933 | 1.2232 | 1.0802 |
| J-NOPAIR | 2.2061 | 1.2968 | 1.7404 | 1.4169 | 1.3396 | 1.2072 |
| J | 2.2021 | 1.3207 | 1.7046 | 1.3877 | 1.3113 | 1.1857 |

Controls: matched-wrong and permuted support above correct at every k in
every arm; wrong-protein above correct; no inversions.

## What this closes

The legal covariate space for the k=0 level is now essentially exhausted:
sequence (frozen and LoRA-tuned), panel composition, journal/publisher
provenance, pocket structure, assay counts - the best trained head reaches
level^2 ~1.30-1.32 (J-NOPAIR) against the 0.1239 budget, and every
mechanism that couples such a head to the k>=1 transport degrades ordering
with resolved intervals. The remaining untested framework family is
contrastive coembedding (ConPLex-style); it operates on ordering structure
and cannot by itself create cross-component level information beyond the
measured shares. MSA remains blocked on a governed UniRef snapshot.
