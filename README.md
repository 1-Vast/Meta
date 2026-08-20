# MetaSieve-DTA

MetaSieve studies trainable meta-learning for zero- and few-shot affinity
prediction on protein targets absent from training.

## Current state

No model has yet met the governed excellence threshold. The programme now has
two parallel tracks: Main Line P measures practical cold-target zero/few-shot
performance, while Main Line M / Core Task 1 decides whether a gain is caused by
transferable protein-conditioned interaction. A mechanism failure does not
invalidate a performance result; it prevents the corresponding mechanism claim
and production integration.

Core Task 1 remains **UNRESOLVED overall**. The latest CIIP-1A oracle-coordinate
local-ESM potential is **ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED** on its governed
49-pair Duong-Ly functional-assay subset. Correct and random local windows were
both nonconstant on 9/9 test pairs, and correct-minus-random pair-mean R2 was
-0.1217 with parent bootstrap [-0.4569, +0.0327]. This is a scoped model failure,
not biological falsification and not a Ki/Kd/pK cold-target DTA result.

A completed read-only audit found ESM contextual propagation outside the
mutation site, but did not fit a predictor; context-only ligand-conditioned
value is therefore **NOT EVALUATED**. CIIP-1B, the BindingDB interaction bridge,
a deployable protein router and production integration are not authorized. Any
mechanism successor requires a new preregistration and explicit authorization.

The earlier BindingDB exact-MMP programme remains read-only evidence: it was not
estimable on the governed double-cold surface and did not establish biological
absence. Practical few-shot baselines may continue independently, but cannot be
described as protein-conditioned interaction without passing the mechanism
controls. Ridge, closed-form adaptation, test-time query labels and ungoverned
data mixing are excluded.

The optional Cartesian module is tested algebraically but not a performance path:
none of the 17,717 current DTA cells has a legal common-frame protein-ligand pose.

## Read first

1. [Current task](task.md)
2. [Evidence ledger](report/EVIDENCE_LEDGER.md)
3. [Current model evidence](report/CURRENT_MODEL_EVIDENCE.md)
4. [Current CIIP-1A verdict](tools/research/stageCIIP_potential_bridge/CONTROL_REPORT.md)
5. [Context-propagation audit](tools/research/stageCIIP_context_propagation_20260820/CONTEXT_PROPAGATION_REPORT.md)
6. [Measured BindingDB boundary](report/BOUNDARY_20260817_NIGHT.md)
7. [Repository organization](docs/PROJECT_FILE_ORGANIZATION.md)

Historical details removed from the working tree remain recoverable through
[the archive index](archive/README.md) and Git history.

## Verification

```powershell
conda run -n drug python -m pytest tools/tests -q
conda run -n drug python -m scripts.audit_research_record --skip-loading
```

Local data availability and dataset governance are documented in
[dataset/README.md](dataset/README.md).
