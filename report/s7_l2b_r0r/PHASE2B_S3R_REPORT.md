# Phase 2B S3R real structural ordinal result

Terminal verdict: `REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED`

| arm | component-macro AP_bidir |
|---|---:|
| candidate | 0.035880 |
| b5diff | 0.031582 |
| foreign | 0.035735 |
| context | 0.032336 |
| permuted | 0.037125 |
| chemistry_shuffle | 0.043652 |
| zero_W | 0.025472 |

| Gate | delta | LCB95 | PASS |
|---|---:|---:|:---:|
| R1_vs_chance | 0.010408 | 0.006920 | False |
| R2_vs_frozen_B5_differential | 0.004298 | -0.001630 | False |
| R3_vs_two_ligand_foreign_pair | 0.000145 | -0.002651 | False |
| R4_vs_residue_context_corruption | 0.003544 | -0.003156 | False |
| R5_vs_trained_permuted_label_learner | -0.001245 | -0.006595 | False |

R6 was not opened because pairwise ordinal training does not identify
absolute amplitude, ligand-feature origin or the difference-design nullspace.
No affinity value was read and the frozen law operator was not modified.
