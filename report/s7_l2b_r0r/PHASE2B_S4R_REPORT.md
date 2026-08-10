# Phase 2B S4R graph-aware ligand representation result

Terminal verdict: `REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED`

| arm | component-macro AP_bidir |
|---|---:|
| candidate | 0.046856 |
| baseline41 | 0.035880 |
| b5diff | 0.031582 |
| foreign | 0.046212 |
| context | 0.027357 |
| ligand_only | 0.025472 |
| permuted | 0.036293 |
| chem_shuffle | 0.051322 |
| zero_W | 0.025472 |

| Gate | delta | LCB95 | margin | PASS |
|---|---:|---:|---:|:---:|
| R1_vs_chance | 0.021384 | 0.016064 | 0.05 | False |
| R2_vs_frozen_B5_differential | 0.015273 | 0.008173 | 0.03 | False |
| R3_vs_two_ligand_foreign_pair | 0.000644 | -0.009226 | 0.03 | False |
| R3b_vs_ligand_only | 0.021384 | 0.016064 | 0.03 | False |
| R4_vs_residue_context_corruption | 0.019498 | 0.013599 | 0.03 | False |
| R5_vs_trained_permuted_label_learner | 0.010563 | 0.003880 | 0.05 | False |

| non-gating contrast | delta | LCB95 |
|---|---:|---:|
| C1_candidate_minus_baseline41 | 0.010976 | 0.004939 |
| C2_baseline41_minus_chance | 0.010408 | 0.006920 |

Heldout-B was neither created nor read. R6 was not opened. No affinity
value was read and the frozen law operator was not modified.
