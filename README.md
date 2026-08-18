# MetaSieve-DTA

MetaSieve studies trainable meta-learning for zero- and few-shot affinity
prediction on protein targets absent from training.

## Current state

No model has yet met the governed excellence threshold. The active work is
Stage X, a new Core-Task-1 qualification cycle using governed public
matched-variant and single-platform kinase panels. X0-D data acquisition is
complete enough to continue, but X0 has not passed: the first local
representation-capability implementation and planted-signal draft require
correction before any biological or model claim. See
`report/STAGE_X_ROUND1_REVIEW_20260818.md`.

The earlier BindingDB exact-MMP programme remains read-only evidence: it was
not estimable on the governed double-cold surface and did not establish
biological absence. Relative SAR, level/shape training, few-shot adaptation and
final confirmation remain gated behind a valid transferable protein signal.
Ridge, closed-form adaptation, test-time query labels and ungoverned data mixing
are excluded.

The optional Cartesian module is tested algebraically but not a performance path:
none of the 17,717 current DTA cells has a legal common-frame protein-ligand pose.

## Read first

1. [Current task](task.md)
2. [Stage X round-1 review](report/STAGE_X_ROUND1_REVIEW_20260818.md)
3. [Current evidence](report/CURRENT_MODEL_EVIDENCE.md)
4. [Measured BindingDB boundary](report/BOUNDARY_20260817_NIGHT.md)
5. [Evidence ledger](report/EVIDENCE_LEDGER.md)
6. [Repository organization](docs/PROJECT_FILE_ORGANIZATION.md)

Historical details removed from the working tree remain recoverable through
[the archive index](archive/README.md) and Git history.

## Verification

```powershell
conda run -n drug python -m pytest tools/tests -q
conda run -n drug python -m scripts.audit_research_record --skip-loading
```

Local data availability and dataset governance are documented in
[dataset/README.md](dataset/README.md).
