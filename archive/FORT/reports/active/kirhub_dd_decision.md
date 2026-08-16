# KirHub mutation-anchored double-difference audit

Date: 2026-07-26  
Decision: **signal present, confirmatory gate blocked; do not implement a new loss yet**

## Scope correction

The public workbook does not match the proposed “61 inhibitors × 369 kinases” dose-response
design. Its relevant contents are:

- Table S3: 14 inhibitors × 10 concentrations × 369 wild-type kinases.
- Table S4: 92 inhibitors × 409 wild-type kinase preparations at 1 µM.
- Table S13: 92 inhibitors × 349 mutant/fusion preparations at 1 µM.
- Tables S4 and S13 contain published aggregate residual-activity values, not the underlying
  duplicate observations.

The primary audit therefore uses Tables S4/S13 only. It conservatively keeps point-mutant
constructs with an explicit wild-type counterpart, excludes fusion/ambiguous constructs, and
does not treat construct–ligand pairs as independent observations.

## Predefined estimand and units

\[
D_{g,m,i,j} =
(r_{g,m,i}-r_{g,m,j})-(r_{g,WT,i}-r_{g,WT,j})
\]

where \(r\) is percent residual activity at 1 µM. Primary non-saturation bounds are 5–95%;
10–90% is the sensitivity analysis. A ligand pair is rank-informative only when both the
wild-type and mutant within-pair separations are at least 10 percentage points.

Inference is macro-aggregated at HUGO gene level. KLIFS families and kinase groups are used
only as stricter holdout-unit counts. Pair and construct counts are descriptive.

## Result

| Quantity | 5–95% | 10–90% |
|---|---:|---:|
| Eligible point-mutant constructs (≥20 ligands) | 222 | 134 |
| Eligible genes | 34 | 26 |
| Eligible kinase families | 22 | 18 |
| Eligible kinase groups | 7 | 5 |
| Informative ligand pairs (descriptive) | 75,596 | 28,778 |
| Gene-macro rank-reversal rate | 0.1133 | 0.1396 |
| Rank-reversal 95% LCB | 0.0865 | 0.0974 |
| True-WT minus wrong-WT Spearman advantage | +0.5010 | +0.5216 |
| Pairing advantage 95% LCB | +0.4481 | +0.4630 |
| Median gene absolute double difference | 15.33 pp | 15.93 pp |

All four signal-presence checks pass. The effect also survives the narrower non-saturation
window, so it is not explained solely by values at 0/100% assay ceilings.

## Why this does not authorize a model run

The target success criterion is a paired Spearman gain of at least +0.03. With 34 independent
genes, the normal-approximation 80%-power detection envelope is:

- +0.024 if the future per-gene paired standard deviation is 0.05;
- +0.048 if it is 0.10;
- +0.096 if it is 0.20.

At family level (22 units), those limits are +0.030, +0.060, and +0.120. Thus a +0.03 claim is
only identifiable under an unusually low paired variance. The public workbook does not expose
raw duplicates, so that variance assumption cannot be checked from the released values.

Implementing a double-difference loss now would produce an exploratory result whose main gate is
underpowered. The correct next step is to obtain raw duplicate measurements or an independent
mutation panel, then pre-register a gene-held-out comparison. If exploratory work is explicitly
accepted, the first and only implementation should be a minimal double-difference loss ablation
against ligand-only, shuffled-protein, wrong-WT, and permuted-support controls.

## Reproduction

```powershell
D:\anaconda\envs\drug\python.exe research\kirhub_dd_audit.py
D:\anaconda\envs\drug\python.exe -m pytest tests\test_kirhub_dd_audit.py -q
```

Machine-readable result: `reports/active/kirhub_dd_audit.json`.

Sources:

- Saifudeen et al., *Nature Biotechnology* (2026):
  https://www.nature.com/articles/s41587-026-03090-8
- Supplementary workbook:
  https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-026-03090-8/MediaObjects/41587_2026_3090_MOESM4_ESM.xlsx
