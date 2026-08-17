# MetaSieve-DTA

MetaSieve studies trainable meta-learning for zero- and few-shot affinity
prediction on protein targets absent from training.

## Current state

No model has yet met the governed excellence threshold. The retained development
Pareto set is A0/B3/C2; fixed Morgan/Tanimoto residual weighting is the strongest
k>=2 query-specific comparator. R14 closed the ranking-loss axis and localized
the next question to representation identifiability.

The next authorized family is a protein-conditioned, low-rank SAR moment update,
followed conditionally by correlation-preserving counterfactual meta-training.
It is preregistered but not implemented. Ridge, closed-form adaptation, inner
loops, test-time gradients and query-label inputs are excluded.

The optional Cartesian module is tested algebraically but not a performance path:
none of the 17,717 current DTA cells has a legal common-frame protein-ligand pose.

## Read first

1. [Current task](task.md)
2. [Current evidence](report/CURRENT_MODEL_EVIDENCE.md)
3. [Measured boundary](report/BOUNDARY_20260816.md)
4. [Next research plan](report/NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md)
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
