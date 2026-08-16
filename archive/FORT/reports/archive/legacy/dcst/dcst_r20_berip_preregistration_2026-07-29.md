# DCST-R20 Balanced Evidence-Rectangular Interaction Pretraining preregistration

Date: 2026-07-29  
Status: frozen before balanced-manifest construction

## Motivation

R19 established a large high-confidence TRAIN-only rectangular topology but
failed because one target-pair block contributed `38.7527%` of all
rectangles. R20 tests whether a scientifically valid, endpoint-consistent and
provenance-separated balanced manifest remains large after removing that
pseudo-replication.

R20 is a new audit. R19 remains stopped; its concentration threshold is not
changed.

## Frozen source

Reuse the R19 high-confidence definition and TRAIN-only firewall:

- `n_records >= 2`;
- finite `replicate_sd <= 0.30` pK;
- at least two assay IDs or at least two document IDs;
- remove ligands connected to more than 50 eligible targets;
- take the bipartite 2-core.

The numeric `affinity` column remains forbidden.

## Rectangle validity

A four-cell rectangle is eligible only if:

1. all four cells have the same endpoint (`pKi` or `pKd`);
2. the union contains at least two assay IDs and two document IDs;
3. no single document occurs in all four cells;
4. its two targets and two ligands are distinct.

## Frozen hierarchical balancing

Assign each valid rectangle a SHA-256 ordering key from seed `1729`, endpoint,
the sorted target pair, and the sorted ligand pair.

1. retain at most 256 lowest-hash rectangles per target pair;
2. from that result, retain at most 4,096 lowest-hash rectangles per unordered
   homology-component pair.

No label value may affect ordering, inclusion, or a cap.

## Frozen gates

All must pass:

1. before caps: at least 50,000 valid rectangles, 500 target pairs, 150
   targets, and 130 homology components;
2. after both caps: at least 25,000 rectangles, 500 target pairs, 150 targets,
   and 130 homology components;
3. at least 40% of selected rectangles and 200 target-pair blocks cross
   homology components;
4. largest target-pair share at most 1% and largest homology-pair share at
   most 10%;
5. at least 5,000 selected `pKd` rectangles or at least 20% of the manifest
   is `pKd`; otherwise the model preregistration must be pKi-only and the gate
   fails here rather than silently mixing endpoints;
6. zero non-TRAIN rows and zero numeric affinity values loaded.

Failure returns `STOP_BERIP_BALANCED_MANIFEST_INADEQUATE`.

Pass returns `REQUEST_BERIP_STAGE1_MODEL_PREREGISTRATION`. It authorizes only
writing a label-free rectangle-key manifest and drafting a separate CUDA
model preregistration. Development, confirmation, and sealed scoring remain
forbidden.
