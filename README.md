# MetaSieve-DTA

MetaSieve studies trainable meta-learning for zero- and few-shot affinity
prediction on protein targets absent from training.

## Current state

No model has yet met the governed excellence threshold. Fixed Morgan/Tanimoto
residual weighting is the strongest reproducible k>=2 query-specific comparator,
but it is ligand-only. Stage S rejected a global protein-conditioned SAR field:
a shuffled protein reproduced the complete measured protein gain.

Work now follows a five-phase gated programme. Stage T is the active Phase-1
test of whether true matched molecular-pair transformations expose transferable,
affinity-relevant protein conditioning. Relative SAR, level/shape training,
few-shot adaptation and final confirmation remain blocked in that order.
Ridge, closed-form adaptation, test-time query labels and ungoverned data mixing
are excluded.

The optional Cartesian module is tested algebraically but not a performance path:
none of the 17,717 current DTA cells has a legal common-frame protein-ligand pose.

## Read first

1. [Current task](task.md)
2. [Current evidence](report/CURRENT_MODEL_EVIDENCE.md)
3. [Measured boundary](report/BOUNDARY_20260817_NIGHT.md)
4. [Evidence ledger](report/EVIDENCE_LEDGER.md)
5. [Repository organization](docs/PROJECT_FILE_ORGANIZATION.md)

Historical details removed from the working tree remain recoverable through
[the archive index](archive/README.md) and Git history.

## Verification

```powershell
conda run -n drug python -m pytest tools/tests -q
conda run -n drug python -m scripts.audit_research_record --skip-loading
```

Local data availability and dataset governance are documented in
[dataset/README.md](dataset/README.md).
