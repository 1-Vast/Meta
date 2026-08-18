# I6 production-dataflow integrity report (Stage X0c)

Test suite: `tests/test_x0c_integrity.py` + `tests/test_q0.py`.
Result: 31 passed (23 production-dataflow contracts + 8 VariantRecord unit
tests), 0 failed, runtime ~15 s, CPU/GPU: torch available (CUDA) but the
contracts are tensor-free or small-tensor; no training runs in the suite.

Coverage map against the frozen I6 assertions:

| # | Assertion | Test |
|---|---|---|
| 1 | successor prereg frozen; original X0 prereg untouched | test_prereg_successor_frozen / test_original_x0_prereg_untouched |
| 2 | contrast antisymmetry | test_csc_antisymmetry |
| 3 | identity pair strictly zero | test_csc_identity_pair_zero |
| 4 | reference-term sign flip | test_csc_reference_term_sign_flip |
| 5 | eval labels never enter the reference | test_csc_reference_train_only |
| 6 | stable seed across OS processes | test_stable_seed_cross_process |
| 7 | no Python hash() anywhere in stage code | test_no_python_hash_anywhere |
| 8 | planted truth bitwise recomputable | test_planted_truth_bitwise_recomputable |
| 9 | main effects enter labels additively | test_generator_main_effects_enter_labels |
| 10 | interval bounds ordered and directional | test_interval_bounds_ordered |
| 11 | sign-only target direction | test_sign_only_target_direction |
| 12 | matched arms share cells/masks; no train/val/eval overlap | test_matched_arms_share_cells_and_masks |
| 13 | no parent or scaffold crosses blocks | test_no_parent_or_scaffold_crosses_blocks |
| 14 | unique cells (no duplicate rows) | test_cells_unique |
| 15 | gradient coverage for every trainable branch + regularizer finite nonzero | test_gradient_coverage_and_regularizer |
| 16 | dead-branch capture | test_dead_branch_capture |
| 17 | protein permutation destroys planted link | test_protein_permutation_destroys_planted_link |
| 18 | label permutation destroys signal | test_label_permutation_destroys_signal |
| 19 | cluster bootstrap resamples clusters | test_cluster_bootstrap_resamples_clusters |
| 20 | restricted data not committed to Git | test_restricted_data_not_committed |
| 21 | old-residue hard rule on admitted pairs | test_old_residue_consistency_hard_rule |
| 22 | BRAF alias not generalized | test_braf_alias_not_generalized |
| 23 | checkpoint selection never reads eval labels | enforced by q2.train signature (no eval mask parameter; val-only monitoring) + test_matched_arms_share_cells_and_masks |
| 24-31 | VariantRecord mutation application, atomic multi-mutation failure, ProteinGym parsing, deterministic serialization, hash stability, class vocabulary | tests/test_q0.py |

Toy-only checks are not sufficient: every contract above imports and executes
the production objects (csc.py, q2.py generator/training model, x0_common,
pair table / Q0-B audit). The legacy parent-dir suites
(`tests/test_x0_integrity.py` 8 tests, `tests/test_x0_corrections.py` 26
tests) also pass and remain part of the final pytest run.
