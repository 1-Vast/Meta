# Phase 2B S5D estimand and collapse diagnostics

Terminal verdict: `LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED`

Trains nothing; reuses the frozen S4R checkpoints. The S4R verdict
`REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED` is unchanged by this stage.

## D1 ligand-steering collapse

131 heldout-A constructs have at least three eligible pairs.

| quantity | median |
|---|---:|
| `rho_dg`, unit ligand differences (data-side upper bound) | 0.4550 |
| `rho_graph`, candidate residue fields | 0.4793 |
| `rho_base`, baseline41 residue fields | 0.5758 |
| `rho_graph - rho_dg` | 0.0138 |
| true-vs-foreign field cosine | 0.4487 |

Rule: median rho_graph >= 0.8 and median rho_graph >= median rho_dg + 0.1. Collapse confirmed: `False`.

## D2 symmetric-difference conditional estimand

40157 eligible pairs across 107 closure components, from 46818 primary pairs. Median changed residues per pair 7.0, median gain fraction 0.5000.

| arm | component-macro AP_cond |
|---|---:|
| candidate | 0.655030 |
| baseline41 | 0.638830 |
| foreign | 0.655470 |
| permuted | 0.628586 |
| chance | 0.643744 |

| Gate | delta | LCB95 | margin | PASS |
|---|---:|---:|---:|:---:|
| E1_vs_conditional_chance | 0.011285 | -0.007749 | 0.05 | False |
| E2_vs_foreign_ligand_pair | -0.000440 | -0.021814 | 0.03 | False |
| E3_vs_trained_permuted_learner | 0.026444 | -0.002977 | 0.03 | False |

| non-gating contrast | delta | LCB95 |
|---|---:|---:|
| E4_candidate_minus_baseline41 | 0.016199 | -0.005947 |
| E5_baseline41_minus_chance | -0.004914 | -0.023844 |

Heldout-A was already consumed by S3R and S4R, so every number here
is development evidence and none of it confirms anything. Heldout-B
was neither created nor read, R6 was not opened, no affinity value
was read and the frozen law operator was not modified.
