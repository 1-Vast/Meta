# Stage E report — panel-set level head + orthogonal level/shape routing: REJECTED

Status: development evidence, single seed, meta_val read once after freezing.
meta_test sealed (logical exclusion after parsing; 768 cells withheld).
Authorities: RESULT.json per arm under
report/meta_fewshot/stageD_level_panel_20260817/, T2_meta_val.rows.summary.json,
LSP_meta_val.rows.summary.json, LSP_vs_T2.contrast.json,
LSP_PANEL_SHUFFLE.json, PREREGISTRATION.md.

## Verdict

**LSP is rejected by the preregistered gates; nothing is promoted.** G1 fails:
the k=0 MSE gain is -0.2026 [-0.5195, +0.0869] (unresolved) and k=5 degrades
with a RESOLVED interval +0.2269 [+0.0685, +0.4098]. Stop rule S1 fires.
G2 also fails on the means: Spearman/Pearson/CI are lower at every k
(k=0 Spearman -0.0905, CI -0.0447), although the intervals are unresolved.
The trained panel level head is inert: its own level error (1.438 pK^2) is
only 0.10 pK^2 better than a panel-shuffled control (1.539), against the D0
frozen probe's 1.887 vs 5.075 — the head learned almost none of the panel
signal the no-training probe found. The k>=1 degradation is level
interference: level^2 +0.1895 resolved at k=5 while the centered term
+0.0374 is unresolved — panel_level and the Tanimoto transport fit the same
target level twice (the Stage A substitution confound, now resolved in the
harmful direction).

## Measured numbers (frozen meta_val banks, component-weighted, restored pK^2)

| arm | k | MSE | level^2 | centered | Spearman | Pearson | CI | cliff sign |
|---|---|---|---|---|---|---|---|---|
| T2 | 0 | 2.5961 | 1.7314 | 0.8648 | 0.0790 | 0.0871 | 0.5372 | 0.5299 |
| LSP | 0 | 2.3935 | 1.5153 | 0.8782 | -0.0115 | 0.0392 | 0.4925 | 0.5741 |
| T2 | 1 | 1.7712 | 0.9065 | 0.8648 | 0.0790 | 0.0871 | 0.5372 | 0.5299 |
| LSP | 1 | 1.8497 | 0.9715 | 0.8782 | -0.0115 | 0.0392 | 0.4925 | 0.5741 |
| T2 | 2 | 1.3245 | 0.5257 | 0.7988 | 0.1956 | 0.1964 | 0.5780 | 0.5710 |
| LSP | 2 | 1.4484 | 0.6306 | 0.8177 | 0.1562 | 0.1854 | 0.5562 | 0.5809 |
| T2 | 3 | 1.2197 | 0.4283 | 0.7914 | 0.2355 | 0.2285 | 0.5924 | 0.5677 |
| LSP | 3 | 1.3662 | 0.5538 | 0.8125 | 0.1877 | 0.2193 | 0.5658 | 0.5635 |
| T2 | 5 | 0.9859 | 0.2637 | 0.7222 | 0.3141 | 0.3127 | 0.6188 | 0.6086 |
| LSP | 5 | 1.2128 | 0.4532 | 0.7596 | 0.2446 | 0.2808 | 0.5817 | 0.6997 |

Paired LSP-minus-T2 contrasts (component bootstrap, per target across draws):

| k | MSE | level^2 | centered | Spearman | CI |
|---|---|---|---|---|---|
| 0 | -0.2026 [-0.52, +0.09] | -0.2161 [-0.55, +0.09] | +0.0135 | -0.0905 | -0.0447 |
| 1 | +0.0785 | +0.0651 | +0.0135 | -0.0905 | -0.0447 |
| 2 | +0.1239 | +0.1049 | +0.0190 | -0.0394 | -0.0217 |
| 3 | +0.1465 | +0.1255 | +0.0210 | -0.0477 | -0.0266 |
| 5 | +0.2269 [+0.07, +0.41] | +0.1895 [+0.05, +0.36] | +0.0374 | -0.0695 | -0.0370 |

Controls (both arms, label-bound and protein-sensitive, no inversions):
matched-wrong and permuted support MSE >> correct at every k; wrong-protein
MSE > correct at every k. T2's k=5 already sits at 0.9859 with correct
support dependence.

## Why it failed (measured, not speculated)

1. The D0 panel signal did not survive end-to-end training: the trained head's
   panel association is ~0.1 pK^2 versus the frozen probe's ~3.2 pK^2.
2. The level term and the support transport both consume the target mean; at
   k>=1 their interference costs more than the k=0 gain (resolved at k=5).
3. The orthogonal routing itself did not improve shape (centered +0.014 to
   +0.037, unresolved) — level contamination of the contact dictionary was
   NOT the binding constraint on within-target ordering; the representation's
   own ligand-varying information is.

## Attribution ablations (preregistered, single seed, frozen meta_val banks)

| arm | k=0 MSE | k=0 level^2 | k=0 centered | k=5 MSE | k=5 Spearman | k=5 CI |
|---|---|---|---|---|---|---|
| T2 (baseline) | 2.5961 | 1.7314 | 0.8648 | 0.9859 | 0.3141 | 0.6188 |
| T2-LEVEL (loss-only) | 3.4431 | 2.6361 | 0.8070 | 1.0240 | 0.3140 | 0.6154 |
| LSP-NOROUTE (framework-only) | 2.3259 | 1.4113 | 0.9146 | 1.2409 | 0.2511 | 0.5927 |
| LSP (both, routed) | 2.3935 | 1.5153 | 0.8782 | 1.2128 | 0.2446 | 0.5817 |

Attribution:

- The panel level head carries the k=0 level gain (T2 1.7314 -> LSP-NOROUTE
  1.4113, -0.32) whether or not routing is used, but its substitution with the
  Tanimoto transport costs more at k>=1 (k=5 1.24 vs 0.99).
- Routing (I2) helps centered only marginally (0.9146 -> 0.8782) and never
  beats the baseline's 0.8648: level contamination of the contact dictionary
  was not the binding constraint on shape.
- The level term without the panel head (T2-LEVEL) damages level (1.7314 ->
  2.6361) while slightly helping centered (0.8648 -> 0.8070): the level
  objective routed into the shape path is actively harmful, the opposite of
  the incumbent's implicit level training through smooth_l1.
- Ranking is worse than T2 in every LSP variant at every k (means); none of
  the differences in Spearman/CI are resolved except the k=5 MSE degradation.

The stage's two innovations are therefore attributed: the framework gives a
real but small k=0 level gain with a k>=1 interference cost; the training
innovation gives a small centered gain that does not recover the baseline's
shape term. G1/G2 fail; S1 fires; nothing is promoted.

## What remains true after this stage

- T2 (leak-free retrain) reproduces the incumbent band: k=0 2.596, k=5 0.986
  (below 1.00 already, with honest controls).
- Every k>=1 few-shot gain remains level calibration through support labels;
  the shape term has still never moved by a resolved amount in any arm.
- Best legitimate zero-shot level predictors remain ESM-650M linear (1.6875)
  and the D0 panel probe (1.887); the trained panel head reaches 1.438 on
  level output but its substitution with transport erases the k>=1 gains.

## Next direction (per the loop: a genuinely different hypothesis, not a fix)

The two measured frontiers: (a) zero-shot level is assay-history-dominated
(70% in-fold document variance, 6.8% out-of-document transfer, 23.9% panel
transfer, 11.9% sequence transfer); (b) within-target ordering information in
the trunk representations is capped near r ~ 0.22. The next stage must change
the INPUTS or the INTERACTION FAMILY, not the loss: candidates are the
ESM-650M residue-input trunk (external representation), the structure/pocket
lane (external data; needs a legal source), or a pairwise (query, support)
learned operator that uses the Stage L directional SAR signal directly.
All three are Stage F preregistration material; none runs before being
preregistered.
