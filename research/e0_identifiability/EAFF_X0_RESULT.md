# E-AFF-X0 Crossed Target-Ligand Census Result

## Verdict

```text
STOP_SOURCE_INTERACTION_UNDERDETERMINED
```

The governed ChEMBL37 source contains many nominal crossed rectangles but too
few independent panel/closure components to identify a real affinity
interaction at the preregistered design sensitivity.

## Label Firewall

X0 used `152,737` governed activity IDs from the label-blind E0 input and
queried only document, endpoint, assay, context and protocol-parameter metadata
from the pinned ChEMBL37 SQLite. The SQL selected zero affinity value fields.
No pAffinity, published value, pChEMBL, DAVIS or recipient label was read.

## Results

| Endpoint | Panels | Nominal rectangles | Effective components | Replicate-supported components | Required |
|---|---:|---:|---:|---:|---:|
| Ki | 597 | 1,059,169 | 36 | 18 | 245 |
| Kd | 34 | 232,875 | 12 | 4 | 245 |

Ki contains 53,673 eligible cells, 235 targets and 22,258 ligands. Kd contains
4,995 cells, 75 targets and 1,537 ligands. These raw counts are not independent
sample sizes.

After joining rectangles that share a structured panel or D1
homology-document closure component, the largest single dependency component
contains 56.49% of all Ki rectangles and 76.76% of all Kd rectangles. The
million-scale nominal count is therefore strong pseudoreplication rather than
million-scale crossed evidence.

The replicate-supported stratum is even smaller. Consequently neither a
noise-aware X1 nor an additive-null-only X1 reaches the frozen lower-bound
requirement of 245 independent components.

## Meaning

X0 stops the current ChEMBL interaction-identification route before any label
look. It does not show that protein-by-ligand affinity interaction is absent.
It shows that this governed source cannot distinguish interaction variance from
panel, assay, document and homology dependence at the registered sensitivity.

Therefore H0C cannot be reinterpreted as proof that the radial basis lacks
biology. X1, X2, T-ANISO, RFSA and production integration remain unauthorized.
If the project continues, the next admissible action is a separately registered
label-blind census of a genuinely crossed source/selectivity corpus, not a
weaker ChEMBL split or a lower threshold.

## Audit

The independent audit reproduced every panel, rectangle, dependency and
endpoint count exactly; verified all existing artifact hashes; confirmed that
the SQL and output artifacts contain no affinity value fields; and confirmed
zero DAVIS/recipient reads.
