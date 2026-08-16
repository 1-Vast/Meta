# DCST-R19 ERIP feasibility decision

Date: 2026-07-29  
Decision: `STOP_ERIP_RECTANGULAR_EVIDENCE_INADEQUATE`

## Result

Five of six frozen gates passed. The candidate is large and
provenance-corroborated, but its raw rectangle population is dominated by one
target-pair block, so R19 does not authorize affinity loading or model
training.

| Quantity | Observed |
| --- | ---: |
| High-confidence pairs | 19,712 |
| Targets / ligands / homology components | 441 / 15,082 / 408 |
| Promiscuous ligands excluded | 3 (192 rows) |
| Bipartite 2-core edges | 7,500 |
| 2-core targets / ligands / homology components | 237 / 3,101 / 210 |
| Exact rectangles / target-pair blocks | 527,654 / 994 |
| Rectangle-participating targets / homology components | 224 / 199 |
| Provenance-qualified sample | 99,897 / 100,000 = 99.897% |
| Largest ligand edge fraction | 0.44% |
| Largest target-pair rectangle fraction | 38.7527% |

The largest target-pair block exceeds the frozen 5% limit by a factor of
`7.75`. Treating every rectangle as an independent Stage-1 example would
therefore create pseudo-replication and allow one assay panel/family block to
determine the objective.

## Firewall

Only ChEMBL TRAIN rows were loaded. The `affinity` column was not requested or
loaded. Development, confirmation, and sealed rows were not loaded or scored.
The audit was CPU-only metadata/graph computation and completed in `1.886 s`.

## Consequence

R19 is stopped without model fitting. The positive topology result supports a
separately preregistered balanced successor, but does not permit changing the
R19 gate post hoc. A successor must additionally require endpoint-homogeneous
four-cell contrasts, exclude rectangles that can be wholly attributed to one
document, and cap both target-pair and homology-pair contributions before any
affinity value is loaded.

Authoritative machine result:
`reports/active/dcst_r19_erip_seed1729.json`.
