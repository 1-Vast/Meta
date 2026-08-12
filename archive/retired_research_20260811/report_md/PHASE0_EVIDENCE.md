# Phase 0 — few-shot episode feasibility, evidence report

## Terminal verdict

```text
FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE
```

No model was preregistered and none was trained. The stage stopped at the
earliest failed precondition, as required.

Corpus status marker: `DEVELOPMENT_ONLY_CLOSED_COMPONENT_CORPUS`.

## 1. What was run

A label-blind census of the governed BindingDB Articles 202608 Ki corpus. The
`pK` field was stripped on read and never entered any statistic, so
**affinity label reads for this stage are zero**. Ki is the only endpoint; no
Kd, Kdapp, IC50, inhibition or displacement value was touched.

Feasibility thresholds were declared inside
`research/meta_fewshot/phase0_episode_census.py` and committed (`f25e57b`)
before the census ran, so the outcome cannot have been produced by choosing a
threshold afterwards.

## 2. The design is clean

| axis | train | development | shared |
|---|---:|---:|---:|
| targets | 442 | 68 | **0** |
| ligands | 4,144 | 248 | **0** |
| scaffolds | 1,772 | 107 | **0** |
| documents | 295 | 20 | **0** |
| protein homology-40 groups | 154 | 25 | **0** |

Leakage is exactly zero on every audited axis. The unseen-target requirement,
scaffold disjointness, document disjointness and homology disjointness all hold
by construction. This is a positive finding about the existing closure
governance and is not the reason the stage failed.

## 3. The source side is sufficient

```text
source targets                                442
source targets usable at k=1/2/3/5   336 / 281 / 257 / 220
      of which scaffold-disjoint     328 / 266 / 227 / 182
median ligands per source target                7
```

Episodic training of `U` and `w0` over target-wise support/query episodes is
comfortably supported by the source split. Nothing about the meta-learning
*training* design is blocked.

## 4. The evaluation side is the binding constraint

```text
development targets                            68
development targets usable at k=1/2/3/5  24 / 19 / 18 / 16
      of which scaffold-disjoint          24 / 18 /  9 /  8
median ligands per development target            3
development cells                              535   (4.3% of the corpus)
development dependency components               12
```

Only **16** held-out targets can carry a `k=5` episode, and only **8** of those
with a scaffold-disjoint support set. The declared requirement was 30.

## 5. Power

```text
MDE_d = (z_0.95 + z_0.80) / sqrt(N)

primary unit    held-out target      N = 16   MDE_d = 0.622
secondary unit  dependency component N = 12   MDE_d = 0.718
declared maximum acceptable                   MDE_d = 0.600
```

At `k=5` the panel can only resolve a standardized effect of `0.62` on the
decisive `correct support > zero adaptation` contrast. That is a large effect
to demand of a first few-shot witness, and the requirement was set before the
count was known. Both failing checks — 16 targets against 30, and `0.622`
against `0.600` — are close to their thresholds, which is precisely why they
must not be moved after the fact.

| check | observed | required | result |
|---|---:|---:|:---:|
| evaluation targets at k=5 | 16 | >= 30 | **FAIL** |
| source targets at k=5 | 220 | >= 100 | PASS |
| evaluation dependency components | 12 | >= 5 | PASS |
| target leakage | 0 | 0 | PASS |
| powered at primary unit | 0.622 | <= 0.600 | **FAIL** |

## 6. Failure localization

The cause is **insufficient independent target episodes in the evaluation
split**, and specifically the depth of held-out targets at larger `k`.

Explicitly *not* the cause, with evidence:

- **endpoint or assay heterogeneity** — excluded. One endpoint (Ki), one scale,
  one contract; no modality mixing occurred.
- **target, scaffold, document or homology leakage** — excluded. All five
  leakage counts are zero.
- **source-side episode supply** — excluded. 220 source targets support `k=5`.
- **pipeline or contract validity** — not implicated. The corpus manifest,
  file hashes and feature manifest all resolve, and the census reproduced the
  manifest cell count exactly.

Explicitly **untested**, and therefore not claimable in either direction:

- whether the frozen 288D T-BASIS contains the required signal;
- support rank sufficiency, conditioning, and query row-space coverage;
- optimizer or implementation adequacy;
- a ligand-only shortcut;
- whether target-specific coefficient heterogeneity exists at all.

The last point matters most: **this run did not test the scientific hypothesis.**
`TARGET_COEFFICIENT_HETEROGENEITY` remains exactly as unresolved as it was
before, and nothing here counts as evidence against it.

## 7. Structural note on any re-split

The obvious remedy — component-wise cross-validation over all 31 dependency
components instead of the fixed 19/12 split — would raise the number of
evaluable targets. It is **not** applied here, because it changes the frozen
split contract and requires its own preregistration.

It is also not a guaranteed fix. The largest dependency component holds
**85.86%** of all cells. Any component-respecting partition must either place
that component in training, leaving a small evaluation fold, or hold it out,
leaving little to train on. The corpus has roughly 31 nominal independent units
but is dominated by one, so re-splitting redistributes the problem rather than
removing it.

That is why the registered outcome for this verdict is to acquire or govern a
better open target-panel corpus, not to re-cut this one.

## 8. Governance

- Terminal verdict: `FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE`.
- Affinity label reads: **0**. `pK` stripped on read.
- Endpoint: Ki only.
- Training: none. No preregistration of a model was written; Phase 1 and
  Phase 2 were not entered.
- No threshold was lowered, no requirement relaxed, no split re-cut.
- Corpus manifest sha256 `89072e0a557d3d6f4b9c8189f0a38915d3db689cdb87bff0687416000125a7a9`.
- Frozen 288D T-BASIS features sha256
  `0bc94a70c40780bfba0046b166ededa6eb7855361d55c11ed222a8328c48a03c`, 288
  dimensions, arms `correct`, `foreign_ligand`, `deranged_protein`.
- All code remains under `research/meta_fewshot/`. `model/`, production
  `scripts/`, production `z`, CSMO, Band, the mesh and `A(F,z)=K(B(z)F(z))` are
  unmodified.
- DAVIS, KIBA, recipient and external confirmation labels remain closed.
