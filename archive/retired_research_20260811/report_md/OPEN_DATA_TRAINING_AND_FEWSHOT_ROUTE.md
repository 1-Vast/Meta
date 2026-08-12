# Open-data training and few-shot route

Updated: 2026-08-10.

## Outcome

MetaSieve can now execute governed development training on open affinity data.
This does not mean that an affinity statistic has been identified.

BindingDB Articles 202608 produced 320 usable Ki panels with 12,457 averaged
cells. Strict document, protein-40%-identity and Murcko-scaffold union closure
left 31 components. The largest component contains 85.86% of cells, so the
corpus is suitable for development optimization but not population inference.
The frozen development split contains 12 components and quotient rank 220.

The first minimal witness used the structurally validated 288D T-BASIS and one
shared linear response. It failed:

```text
correct quotient RMSE              0.580314
zero-interaction RMSE              0.580520
explained fraction                 0.000709
correct - zero loss reduction      0.000239 [-0.000981, 0.001496]
correct - foreign loss reduction   0.000870 [-0.001045, 0.002847]
correct - deranged loss reduction -0.000817 [-0.003419, 0.001498]
```

Terminal verdict:

```text
CQ_TBASIS_LINEAR_AFFINITY_WITNESS_NOT_OBSERVED
```

This rejects one population-shared linear affinity direction on the fixed
radial basis. It does not reject target-dependent coefficients on that basis,
nor does it establish that a more complex representation would work.

## Few-shot objective

The intended model remains few-shot drug-target affinity prediction for an
unseen target. Open datasets must therefore be assigned different roles rather
than pooled into one absolute-label table.

1. BindingDB curated Articles: quantitative Ki/Kd quotient and calibration.
2. Kinobeads, PKIS and PKIS2: dense within-panel ordinal/selectivity pretraining;
   their measurement modalities must not be called Ki/Kd.
3. Klaeger full-dose Kdapp: dense quantitative kinase calibration, closed
   against overlapping compounds and publications.
4. PDSP Ki: non-kinase development stratum.
5. Davis and any recipient panel: frozen external evaluation only unless a new
   authorization explicitly changes their role.

## Highest-potential minimal next model

Do not add another PLM, GNN, attention stack, pose branch or knowledge graph.
Test coefficient heterogeneity on the same frozen biological coordinates:

```text
q_t(P,L) = w0^T phi(P,L) + a_t^T U^T phi(P,L)
U in R^(288 x d), d <= 5
```

Dense profiling targets provide source-task coefficients `a_t`. For an unseen
target, only support-observable directions may be estimated. With support
matrix `X_S`, adaptation is restricted to `row(X_S)`; report rank,
conditioning and query coverage, and return zero adaptation/abstention outside
that coverage. No adapted intercept is added.

Training order is fixed:

1. register and build cross-source role, licence and overlap manifests;
2. ordinal pretrain `U,w0` on dense profiling panels, balancing panel then
   dataset then closure component;
3. add quantitative BindingDB/Klaeger constraints without mixing modalities;
4. evaluate target-held-out and scaffold-held-out `k=1/2/3/5` episodes against
   population-only, foreign-support, permuted-support and ligand-only controls;
5. only after correct support improves both ligand-only and wrong-protein arms
   may a bounded low-dimensional statistic be proposed for biological `z`.

The frozen operator `A(F,z)=K(B(z)F(z))` remains unchanged. This route makes
the biological interaction coordinates and mathematical support-identifiability
constraint equally load-bearing, while keeping the trainable model small.

## Engineering status

- BindingDB metadata projection, exact-value extractor, corpus builder, frozen
  ESM2/P1B/T-BASIS feature generator and complete-panel quotient trainer run in
  the `drug` environment.
- Feature generation processed three 12,457-cell arms in 117.14 seconds
  (106.35 cells/s) on CUDA; the feature step itself read zero affinity values.
- The linear solve is intentionally CPU-based because it has only 288
  coefficients. GPU utilization is not a success criterion for this solve.
- Production `model/`, biological `z`, CSMO, Band and the law operator were not
  modified.
