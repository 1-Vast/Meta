# Stage R1: a governed double-cold development protocol

Artifacts: `dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1/`
(`manifest.json`, `assignment.json`, `ligand_similarity.json`), built by
`scripts/build_double_cold_split.py`. The existing CD-HIT40 protocol is
untouched and remains the default everywhere.

Assignment hash `9d8c7289c1b6162f0e39c0c7ff2222bb45305fe3193bfee7b0b8214c0baf5684`,
frozen before any training. `QPSMPData(..., split_directory=...)` recomputes and
verifies that hash on load and refuses to run if it moves.

## Why

Stage R0 established that the old protocol cannot separate chemical
generalisation from ligand recall: 48.9% of its `meta_val` k=0 query cells
contain a ligand that appears verbatim in `meta_train`, and that single fact
accounted for the entire Stage 10 result. Every architectural decision in this
project so far was made on that population.

## Construction

Two axes, both label-blind. Nothing in the builder reads `pK`.

* **protein**: `protein_group_40` (CD-HIT40) components, unchanged;
* **ligand**: Bemis-Murcko scaffold clusters. A cell is evaluation only if its
  component **and** its scaffold cluster are both on the evaluation side;
  off-diagonal cells are discarded.

Then, in order: **document/assay closure** — an evaluation cell is dropped if
any of its `panel_ids` documents also occurs in the training block; and
**eligibility** — an evaluation target is kept only if it retains at least 9
unique ligands, so the nested k=5 protocol with a real query panel is defined.

**Rejected alternative.** Single-linkage Morgan Tanimoto clustering at 0.4 was
tried first and rejected: it chains, and the largest cluster absorbs **74.6%**
of the 9,880 corpus ligands. No two-axis assignment can be built on it. The
`< 0.4` guarantee is therefore obtained by *measurement* per evaluation ligand
rather than by construction, which is the honest version of the same control.

## The protocol

| block | cells | targets | components | ligands | scaffold clusters | median ligands/target |
|---|---:|---:|---:|---:|---:|---:|
| `meta_train` | 5,643 | 346 | 258 | 3,825 | 1,628 | 9 |
| `meta_val` (development) | 1,411 | **41** | **19** | 1,014 | 410 | 21 |
| `meta_test` (confirmation) | 768 | **22** | **10** | 518 | 223 | 22 |

Both evaluation blocks exceed the protocol requirement of >= 6 independent
components and >= 30 eligible targets on the development side, and every
evaluation target supports k=5 with a query panel.

## Closure verification

Recomputed from the assignment, not asserted:

| block against `meta_train` | exact ligand | scaffold cluster | protein component | document |
|---|---:|---:|---:|---:|
| `meta_val` | **0** | **0** | **0** | **0** |
| `meta_test` | **0** | **0** | **0** | **0** |

`meta_val` and `meta_test` also share **0** protein components.

## Measured chemical similarity to the training block

Per evaluation ligand, maximum Morgan(r=2, 1024) Tanimoto to any training-block
ligand:

| tier | `meta_val` ligands | `meta_test` ligands |
|---|---:|---:|
| **< 0.40 (low-similarity tier)** | **827 (81.6%)** | **456 (88.0%)** |
| 0.40 - 0.60 | 172 | 57 |
| 0.60 - 0.80 | 14 | 5 |
| >= 0.80 | 1 | 0 |

For comparison, the old protocol's k=0 bank was 48.9% exact-identity overlap and
the largest single stratum was near-duplicate. The one `>= 0.80` ligand here is
a genuine analogue with a different Murcko scaffold; it is retained and
reported rather than silently removed.

## What it costs

44.1% of the corpus is retained; 9,109 cells are lost off-diagonal, 567 to
document closure and 219 to eligibility. The training block falls from 12,633
cells to 5,643 and its median target has 9 unique ligands rather than 32.
Episodes are therefore thinner, and `draw_episode` gained an optional
`min_query_size` so that a caller using within-target centered or ranking
objectives can require a panel wide enough for them to be defined.

This is a real reduction in training signal and it will show up as a level
shift in absolute numbers. It does not affect any comparison made **within**
this protocol, which is what every Stage R2-R5 decision is.

## Status

* `meta_val` here is the **development** population for Stages R2-R4.
* `meta_test` here is **untouched** and is the confirmation population for
  Stage R5. It has never been read.
* Checkpoints trained under the old CD-HIT40 protocol may be used for
  diagnostics only. Any double-cold claim requires training under this split.
