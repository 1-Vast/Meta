# Gauge-fixed ensemble SSL test-plan review

Updated: 2026-08-09.

## Status of the supplied plan

`METASIEVE_SSL_TEST_PLAN_2026-08-08.md` is explicitly a prospective protocol,
not an experimental artifact. It improves the proposed decision tree, but it
does not supply the missing S5-S9 source commits, manifests, checkpoints,
predictions, or metric files. Formal status remains:

```text
S5_REGISTERED_NOT_RUN
EXTERNAL_S5_S9_CLAIMS_NOT_REPRODUCED
GAUGE_FIXED_ENSEMBLE_SSL_NOT_AUTHORIZED_FOR_TRAINING
```

The plan is accepted as a conditional research proposal after the corrections
below. It does not replace the active S5 preregistration.

## Strong parts retained

1. A0/A1/A2 separate single-pose structure, ensemble information, and
   quantitative label calibration.
2. U1/U2/U3 distinguish matched apo structures, relaxed-unbound proxies, and
   bound-only trajectories.
3. Label-late join, heterogeneous closure graphs, whole-series components, and
   component-level inference are mandatory.
4. Kd and Ki remain separate; IC50 is excluded from the first test.
5. ECFP retrieval, additive encoders, wrong protein, chemistry perturbation,
   mean-structure, unbound swap, and support controls are retained.
6. Failure at an earlier Gate stops later training and production admission.

## Required corrections

### 1. S5-S9 are not a verified starting point

G1 cannot be called an S8 reproduction until the exact S8 implementation,
thresholds, split, control maps, and original artifacts are recovered. Without
them, A0 is a newly registered baseline and must obtain its own development and
sealed confirmation blocks.

### 2. `d_z <= 5` has two different meanings

The frozen theory limits support-derived member identity to at most `k`
continuous dimensions. It does not prove that every population statistic must
have dimension at most five. A five-coordinate population state is an
engineering/capacity restriction for the exact tensor sieve. These claims must
remain separate.

### 3. Not every proposed quantity is a statistic coordinate

- `wrong-protein margin` depends on a nuisance-control distribution and a pose
  search budget. It is an admission metric, not a deployable biological state.
- `coverage` is a validity predicate/partiality flag. It must travel beside the
  statistic, not be treated as a continuous biological coordinate.
- frame-energy entropy is gauge-invariant only to additive shifts after a
  temperature and normalization convention is frozen; its physical entropy
  interpretation is not automatic.
- `-tau*logmeanexp(-E/tau)` is an aggregation of learned scores until the score
  scale, reference measure, and gauge are identified. It must not be labelled a
  free energy before then.

The initial candidate state is therefore at most a set of unnamed, executable
observables plus separate validity flags. Biological names are earned by
interventions, not assigned by architecture.

### 4. G5 is not executable with the current context registry

The current `context_index` uses abstract coordinates 13, 17, and 27 of the
28-dimensional placeholder interface. A new `d_z <= 5` exact sieve can be
constructed, but it cannot be assembled through the existing `BandOperator`
context map without a separately admitted context registry. No index remapping
is authorized before a biological statistic passes G4.

### 5. A2 is empirical calibration, not a physical gauge proof

Within-series `Delta log K` removes a target-level offset but does not remove
ligand-dependent assay error. A small absolute source loss aligns the model to
an experimental label scale; it does not reconstruct a universal standard-state
free energy. Endpoint calibration must be separate from the mechanistic state,
and assay/source IDs may not enter that state.

The admissible wording is `SSL-pretrained, source-difference-calibrated`, not
`thermodynamic gauge identified`, until independent thermodynamic controls and
cross-source replication pass.

### 6. The edge-direction control is ambiguous

Reversing both an ordered edge and its target difference changes nothing.
Instead, freeze a node-label permutation or an inconsistent sign perturbation
inside each series component while preserving graph topology. Compare with a
node-potential model and report the induced cycle inconsistency. The exact
perturbation must be frozen before evaluation.

### 7. G2 controls must respect U1/U2/U3

An unbound-ensemble swap is meaningful only for U1/U2. U3 contains bound
trajectories only and cannot be required to pass an unbound-reference Gate.
Results and coverage must remain stratified; U2 cannot confirm a U1 claim.

### 8. Few-shot theory requires an outer object

Constraining a support vector to at most `k` dimensions is necessary but not
sufficient for a theory-certified adapter. G6 must also emit conservative outer
envelopes/radii, accept noise tolerance, expose off-coverage partiality, and
pass monotonicity and gauge-invariance tests. Otherwise it is a regularized
few-shot engineering model, not an implementation of the frozen general
identifiability theorem.

## Data reality

- MISATO supplies explicit-water bound-complex MD and quantum descriptors for
  almost 20,000 complexes. It supports U3 and ensemble observability. It does
  not automatically supply matched U1 or absolute binding free energies.
- PLINDER v2 supplies pocket/protein/ligand similarities, matched molecular
  pair/series metadata, and links to apo/predicted structures. Linked apo
  structures are not MD ensembles and require construct/pocket matching.
- The PLINDER repository currently warns that BindingDB affinity annotations
  are disabled because of a parsing bug. Affinity must come from an independently
  governed release.
- The Open Force Field protein-ligand benchmark supplies curated congeneric
  targets and relative-free-energy graph metadata. It is suitable for a small
  physics replication/control, not a source-wide population claim.
- PDBbind-derived rows remain prohibited until redistribution and label
  provenance are independently licensed and sealed.

## Correct execution order

```text
R-VERIFY
  recover S5-S9 evidence or retain external-only status
    |
    v
S5
  actual frozen-P1B pair-local observability
    |
    +-- FAIL with oracle PASS -> separately register pose-aware baseline
    +-- PASS -> freeze structural channel, then source calibration
    |
    v
G0-CENSUS (label-blind metadata only)
  MISATO U3 / PLINDER U1 links / quantitative series graph / closure / licence
    |
    +-- insufficient -> NOT_RUN_DATA_UNDERDETERMINED
    |
    v
G1 new pose-aware baseline confirmation
    v
G2 ensemble observability, stratified U1/U2/U3
    v
G3 source-difference calibration
    v
G4 fresh absolute-risk admission
    v
G5 separate biological and context registry admission
    v
G6 outer-section few-shot evaluation
```

Only `R-VERIFY`, the existing S5, and a metadata-only G0 census are currently
permitted. A0/A1/A2 GPU training, affinity label joins, production `z`, DAVIS,
and P2-P4 remain frozen.

## Minimal G0 outputs before model code

1. release/license/checksum manifests for every candidate source;
2. counts of usable U1, U2, and U3 systems kept separate;
3. exact holo/apo construct and pocket mapping coverage;
4. matched-series nodes, edges, connected components, endpoints, assays, and
   independent target components without reading numeric values;
5. overlap against all development-exposed PDB, protein, pocket, ligand, and
   scaffold registries;
6. storage, frame count, I/O throughput, and GPU memory smoke estimates;
7. a power plan based on a declared development source or an external variance
   bound, never on the sealed confirmation labels;
8. explicit `NOT_RUN` reasons for undercoverage, licensing, or mapping failure.

No ensemble model should be implemented until this census shows that the
required comparison is estimable.

## Primary sources checked

- MISATO: https://www.nature.com/articles/s43588-024-00627-2
- MISATO archive: https://zenodo.org/records/7711953
- PLINDER data layout and linked apo structures: https://plinder-org.github.io/plinder/dataset.html
- PLINDER repository and known issues: https://github.com/plinder-org/plinder
- OpenFF protein-ligand benchmark: https://github.com/openforcefield/protein-ligand-benchmark
