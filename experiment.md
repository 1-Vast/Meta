# Experimental contract

The active state is summarized in `PROJECT_SUMMARY.md`; the full chronology is
in `history.md`. This file is only a compact contract reminder.

## Estimand

Predict quantitative affinity for query ligands of an unseen target after
observing `k=1/2/3/5` support affinities for that target. Each target is a
meta-learning task. Query labels must never be used for model selection.

## Source and episode rules

- Keep measurement modalities separate.
- Train on target episodes or complete panels, not IID activity rows.
- Use a low-dimensional task basis with `d <= 5`.
- Estimate the target section from support rows only with strictly positive
  ridge regularization.
- Report support rank, conditioning, query coverage and all controls.

## Required controls

```text
support-free population
zero section
correct support
foreign-target support
permuted support labels
ligand-only
wrong-protein / wrong-partner
endpoint-separated summaries
source-separated summaries
```

No raw pair map or arbitrary neural latent enters `z`. Only an independently
confirmed, partner-specific biological statistic may be proposed for a later
bridge to `A(F,z)=K(B(z)F(z))`.
