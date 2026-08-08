# Processed Data

This directory contains data derived from `dataset/raw/`, including validated
protein embedding caches and split-sealed experiment inputs. It is not a source
of record for raw datasets.

Current mechanism artifacts:

- `open_structures/pilot20k_holo_governed_v2`: 14,906 P1A-passing complexes;
- `open_structures/pilot20k_homology_split_v2`: 11,926/1,490/1,490 split;
- `open_structures/pilot20k_structure_supervision_v2`: corrected 128-slot labels;
- `open_structures/pilot20k_esm2_t30_slots128_v1`: frozen P1B ESM-2 bank; and
- `davis_mechanism_p1c_v1`: source/metaval-only ESM banks used by failed P1C.

P1R0/P1R1/P1R2A/P1R2B0/P1R2B1 reused these immutable banks and wrote reports/cache
only; they created no replacement dataset or checkpoint. The P1R2A feature
cache and B0/B1 metrics are report artifacts, not new training datasets. F0
used existing manifests only and emitted no processed affinity corpus.
F0R created the incomplete cache
`source_affinity/chembl37_f0_rehydrated_v1/`: five files passed frozen hash
verification before live-API drift stopped recovery. It is not a training
dataset and must not be treated as one.

The `v2` structure sidecar supersedes earlier compact-order mapping artifacts.
See `report/mechanism_refactor/P1_EXECUTION_REPORT.md` before reusing them.
