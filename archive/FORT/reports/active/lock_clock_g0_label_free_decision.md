# LOCK/CLOCK G0-L decision

Accepted artifact:
`reports/active/lock_clock_g0_label_free.json`, SHA-256
`3548f0a165168fed576b008a71c8473c772f6f5a34722a62bd70f7b7483e8df8`.

The fixed LOCK kernel passed the four registered execution gates on 372 genes, 319 strict homology
components, 93 KLIFS families, and 9 groups:

| quantity | result | gate |
| --- | ---: | --- |
| normalized minimum eigenvalue | `2.07e-16` | `>= -1e-8` |
| family+composition residual energy | `0.3851` | `>= 0.05` |
| centered alignment with pooled ESM-2 | `0.1633` | `<= 0.95` |
| non-constant within-family pair fraction | `0.9992` | `>= 0.80` |

The separate low-dimensional claim failed. The centered effective rank was `289.36`, and the top
16 dimensions retained only `0.2340` of centered Frobenius energy versus the frozen `0.80` bar.
No dimension increase is authorized.

Verdict:
`LOCK_G0_LABEL_FREE_PASS_AUTHORIZE_REORDERING_AUDIT`.

This is a geometry result only. It does not show affinity prediction, ligand reordering, cross-family
transfer, or a structure-conditioned CLOCK effect. No affinity label was read in G0-L.

