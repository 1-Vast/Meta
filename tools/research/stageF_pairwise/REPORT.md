# Stage F report — pairwise learned interaction transport: REJECTED

Development evidence, single seed, meta_val read once after freezing;
meta_test sealed (logical exclusion after parsing; 768 cells withheld).
Authorities: F_vs_T2.contrast.json, FABS_vs_T2.contrast.json,
F_meta_val.rows.summary.json, FABS_meta_val.rows.summary.json, per-arm
RESULT.json, PREREGISTRATION.md.

## Verdict

**Rejected by gates G2/G3; nothing promoted.** Against the frozen T2 baseline,
every F minus T2 interval crosses zero at every k. The only favourable means
are k=5 MSE -0.0076 and k=5 centered -0.0062 — far below any resolved
threshold — while Spearman/CI are lower at k=0,1,3,5 (unresolved). The
framework-only arm (F-ABS) is worse at every k (k=0 +0.1618, k=5 +0.0225,
unresolved). Stop rule S2 applies (G3 fails on the means; G2 fails outright).

## Numbers (frozen meta_val banks, component-weighted, restored pK^2)

| arm | k=0 MSE | k=2 MSE | k=3 MSE | k=5 MSE | k=5 Spearman | k=5 CI |
|---|---|---|---|---|---|---|
| T2 (frozen baseline) | 2.5961 | 1.3245 | 1.2197 | 0.9859 | 0.3141 | 0.6188 |
| F (both innovations) | 2.6192 | 1.3555 | 1.2270 | 0.9783 | 0.2951 | 0.6063 |
| F-ABS (framework only) | 2.7580 | 1.3628 | 1.2505 | 1.0084 | 0.2903 | 0.6121 |

Paired F minus T2 (component bootstrap): k=0 MSE +0.0230 [-0.17, +0.22];
k=5 MSE -0.0076 [-0.04, +0.03]; k=5 Spearman -0.0190; k=5 CI -0.0124.
F-ABS minus T2: k=0 +0.1618 [-0.16, +0.54]; k=5 +0.0225 [-0.01, +0.06].

## Attribution

- The learned pairwise edge logits add nothing over the fixed Tanimoto
  kernel: the supervision term (F minus F-ABS) recovers only the small
  regression the framework itself introduced, and against the frozen
  baseline nothing resolves. This is the fifth learned-kernel family in this
  project's record (R3/R4, R6/R7, A2's moment form, Stage D's level head,
  Stage F's edge operator) that fails to beat the fixed Morgan/Tanimoto
  weighting — now including the Stage L pairwise direction as input.
- Consistent with the repeated observation that the trunk's ligand-varying
  subspace carries r ~ 0.22 at most: no downstream operator, however
  well-matched to the measured signal, can extract more than the
  representation contains.

## Ledger update

Falsified framework families to date: analytic/legacy operators; BPSF/CIPF
pair trunk; contact-grammar trunk with moment-form adaptation (A2); inner/
outer-loop meta-learning (Stage A/B); centered-objective protein
conditioning (Stage P); panel-set level head + orthogonal routing
(Stage D/E); pairwise learned transport (Stage F). External representations
tested: ESM-150M pooled, ESM-650M pooled, panel composition, assay
covariates. Blocked lanes: MSA (no governed UniRef snapshot), structure
pocket (15/499 exact holo coverage locally).
