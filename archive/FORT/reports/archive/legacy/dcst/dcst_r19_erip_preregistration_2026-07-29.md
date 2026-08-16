# DCST-R19 Evidence-Rectangular Interaction Pretraining preregistration

Date: 2026-07-29  
Status: frozen before TRAIN-only feasibility audit

## Motivation

R14-R18 show that no available external structural source covers enough of the
strict ChEMBL dual-cold target and ligand domain. R19 therefore tests a
same-domain two-stage alternative. Stage 1 would train directly on
high-confidence ChEMBL-TRAIN interaction contrasts, not transfer an absolute
PLINDER/KLIFS/BindingDB representation.

For two targets and two ligands, define the rectangular contrast

`tau = y(t1,l1) - y(t1,l2) - y(t2,l1) + y(t2,l2)`.

Target-only and ligand-only main effects cancel exactly. A model that predicts
`tau` must carry pair-specific interaction information. This is a distinct
information objective, not another attention, router, alignment, or prompt
module.

## Audit data and label firewall

Read only rows whose frozen `dual_cold_split` is `train`. The audit may load:

- target, ligand connectivity, scaffold, endpoint, and homology component;
- replicate count and replicate standard deviation;
- assay and document identifiers.

The numeric `affinity` column must not be loaded. Development, confirmation,
buffer, and sealed-test rows must not be loaded. The audit may use
`replicate_sd` only as an already-frozen reliability summary; it may not use
the pair's central affinity value or its sign/magnitude.

## Frozen high-confidence definition

A TRAIN pair is Stage-1 eligible only when:

1. `n_records >= 2`;
2. `replicate_sd` is finite and at most `0.30` pK units;
3. it has at least two distinct assay IDs or at least two distinct document
   IDs.

Before rectangle construction, exclude a ligand if it is connected to more
than 50 eligible targets. This prevents a small set of promiscuous reference
compounds from determining the interaction objective.

Then iteratively take the bipartite 2-core: every retained target and ligand
must have degree at least two.

## Frozen feasibility gates

R19 may request a separate CUDA model preregistration only if all gates pass:

1. high-confidence source: at least 5,000 pairs, 150 targets, 3,000 ligands,
   and 100 homology components;
2. post-promiscuity bipartite 2-core: at least 2,500 edges, 100 targets, 1,000
   ligands, and 80 homology components;
3. exact rectangular topology: at least 10,000 rectangles, 500 target pairs
   with at least one rectangle, and 80 participating homology components;
4. provenance: a deterministic sample of at least 2,000 rectangles has at
   least 80% whose four cells collectively span at least two assays and two
   documents; the extrapolated qualified count is at least 8,000;
5. concentration: the largest target-pair block contributes at most 5% of
   rectangles and the largest ligand contributes at most 2% of eligible
   2-core edges;
6. firewall: zero non-TRAIN rows and zero numeric affinity values are loaded.

Failure returns `STOP_ERIP_RECTANGULAR_EVIDENCE_INADEQUATE`.

Pass returns `REQUEST_ERIP_STAGE1_MODEL_PREREGISTRATION`. A pass does not
authorize development/confirmation scoring. The subsequent model
preregistration must freeze:

- a matched direct-affinity curriculum control;
- shuffled-cell, additive-main-effect, and label-permuted rectangle controls;
- provenance-clustered sampling;
- a CUDA Stage-1 contrast loss followed by full-TRAIN Stage-2 affinity
  fine-tuning;
- a development-only gate against B0 before confirmation remains sealed.
