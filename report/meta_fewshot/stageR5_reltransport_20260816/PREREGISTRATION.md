# Stage 0/1 (R5) preregistration: contract repairs and structural gates

Frozen before any real-data result of the new model existed. This stage
delivers (a) the 2026-08-16 experimental-contract repairs and (b) the
Stage 1 structural and synthetic falsification gates for Core Innovation A
(the protein-conditioned relative interaction potential shared by the
zero-shot shape and the few-shot transport, `model/reltransport.py`).

## Part 1: contract repairs (fixed before any new training)

1. **Wrong-protein donor pool.** At evaluation the donor candidates come from
   the *same* evaluation split (`meta_val` for meta_val targets), so the
   contrast varies protein identity alone and never confounds it with
   seen-versus-unseen. Training-time counterfactual donors stay in
   `meta_train` (both proteins seen). Implemented as two separate pools in
   `matched_donors(donor_pool=..., whitening_pool=...)`.
2. **Whitening fit pool.** The whitening mean and covariance are fitted on
   `meta_train` only, **always** — the evaluation population never contributes
   its own whitening statistics.
3. **Donor descriptions.** Every artifact records its actual donor pools
   (`RESULT.json` `donors` block); the wrong description in the R3/R4
   COMPARE payload (which claimed meta_train donors while the code drew
   meta_val donors with meta_val-fitted whitening) is corrected and the
   corrected comparison is rerun in R6.
4. **Gradient-cosine reporting.** No single-episode/single-step gradient
   snapshot may be reported as a mechanism claim. Every diagnostic step
   records the mean over all of its episodes plus the per-episode records;
   `gradient_summary` aggregates means and conflict frequencies across all
   recorded steps of a seed; seed aggregation happens at the stage level.
5. **meta_test seal.** `meta_test` is dropped physically by default
   (`QPSMPData(include_meta_test=False)`); `train_qpsmp.py` and
   `evaluate_qpsmp.py` require an explicit `--include-meta-test`/`--eval-meta-test`
   opt-in; `evaluate_qpsmp.py` gains `--split` (default `meta_val`) and
   `--split-directory`. The A0 files that were auto-generated on the old
   meta_test remain sealed and are never read.
6. **Artifacts for every new run:** config, split assignment sha256, seed,
   checkpoint sha256, per-target predictions (jsonl), component bootstrap
   (compare script), peak GPU memory, wall time, parameter count, gradient
   coverage census, activation statistics, donors block, meta_test block.

## Part 2: Stage 1 structural gates (`tests/test_reltransport_synthetic.py`)

Each gate is a falsification of a specific way the mechanism could be
cosmetic. All run before real-data training; none reads meta_val labels for
selection.

1. delta antisymmetry (exact); 2. anchor set has exactly zero mean shape;
3. endpoint == ligand_prior + target_level + shape (exact); 4. level branch
constant within target; 5. ligand prior protein-blind; 6. k=0 returns the
endpoint exactly; 7. support permutation invariance; 8. query permutation
equivariance; 9. query independence (no transductive statistic); 10. support
labels enter only as residuals (shift-by-c contract); 11. query labels never
an input; 12. geometry refused; 13. k=1 correction varies with the query;
14. k=1: functional label effect + trainable transport parameters + nonzero
autograd gradient through the support ligand; 15. support-ligand replacement
changes the prediction; 16. no dead trainable branch (k>=2; documented k=0/k=1
semantic exceptions); 17. level/shape gradient separability; 18. synthetic
bilinear interaction task: ordering dies without the interaction branch
(3 seeds, 8 proteins, mean CI >= 0.70, mean gap >= 0.20, positive in >= 20/24
cells); 19. matched-wrong support clearly worse than correct support;
20. private task abstains to the level shift (2 seeds x 3 panels, mean full
<= mean level-only + 0.005).

## Part 3: not in scope here

Performance on `meta_val`, bootstraps, and any innovation claim — those
belong to Stage 2 (R6) screening and Stage 3 (R7) formal development.
`meta_test` remains sealed.
