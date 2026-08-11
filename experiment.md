# Experimental contract

Current execution status and the retained/failed result ladder are in
`PROJECT_SUMMARY.md`. The Phase 0 record below is historical context, not the
current model result.

## Final estimand

Predict quantitative affinity for query ligands of a previously unseen target
after observing `k=1/2/3/5` support affinities for that target. A target is one
meta-learning task; support and query ligands must be scaffold-disjoint where
the data permit.

## Source learning

1. Keep measurement modalities separate.
2. Train on complete panels or target episodes, not IID activity rows.
3. Learn at most one low-dimensional task basis `U`, `d<=5`, on the frozen
   biological candidate `phi(P,L)`.
4. Balance panel, dataset and target contributions; record source provenance.
5. Never use confirmation targets or query labels for model selection.

Quantitative Ki/Kd/Kdapp data constrain within-modality values. Profiling data
with inhibition or displacement endpoints contributes only ordinal/ranking
information. These losses may train the same mechanism coordinates but may not
be numerically pooled as one label.

## Meta-learning episode

For each source target, split ligands into support and query. Learn `U,w0` from
query loss after computing a strictly positive-ridge support section. The
section is constrained to the support feature row space. Report its rank,
condition number and query coverage; return zero correction or abstain when the
query is off coverage.

The first experiment changes only coefficient sharing. It must not add a new
encoder, attention stack, learned support network, geometry branch or affinity
head.

## Required controls

- support-free population prediction;
- zero section;
- correct support;
- foreign-target support;
- within-task permuted support labels;
- ligand-only and wrong-protein biological controls;
- endpoint-separated and source-separated summaries.

Primary inference units are held-out targets or predeclared dependency units,
never ligand rows. Population admission requires positive lower confidence
bounds for correct over ligand-only and wrong-protein predictions, followed by
independent publication/time confirmation.

## Operator boundary

The episodic learner is upstream of `A(F,z)=K(B(z)F(z))`. No raw 288D vector,
pair map or arbitrary neural latent enters `z`. Only a bounded, independently
confirmed mechanism statistic plus support rank/coverage/certificate may be
proposed for a later bridge contract.

## Phase 0 execution record

`research/meta_fewshot/phase0_episode_census.py`, commit `f25e57b`, thresholds
declared in the script before the run.

```text
source      targets 442   usable at k=1/2/3/5   336 / 281 / 257 / 220
evaluation  targets  68   usable at k=1/2/3/5    24 /  19 /  18 /  16
                          scaffold-disjoint      24 /  18 /   9 /   8
evaluation cells 535   components 12   leakage 0 on all five axes
MDE_d  targets k=5 0.622   components 0.718   declared ceiling 0.600
```

Terminal verdict `FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE`. Phase 1 and Phase 2
were not entered; no model preregistration exists and nothing was trained.
Affinity label reads 0 (`pK` stripped on read). Endpoint Ki only.
