# CIIP Contextual-Propagation Audit Preregistration

Date: 2026-08-20. This is a read-only successor to the finalized CIIP-1A
control-arm adjudication. It does not amend any earlier preregistration.

## Question and Scope

For the 49 ESM-covered, verified single-mutation Duong-Ly pairs, measure how a
single residue edit changes frozen ESM-2-150M per-residue states as a function
of distance from the verified edit. This is a representation audit, not a
predictive-model experiment. The assay endpoint remains centered functional
percent inhibition, never Ki, Kd, pK, or DeltaDeltaG.

## Frozen Inputs

- `stageCIIP_potential_bridge/DATA1A.json` and `DATA2X2.json` define pairs,
  split labels, and the 49 covered indices.
- `stageX_csc_signal/stageX0c_measurement_qualification_20260818/q1_esm_cache.npz`
  supplies original WT/variant ESM-2-150M hidden states.
- The locally cached `facebook/esm2_t30_150M_UR50D` tokenizer and model are
  used in evaluation mode only for the erasure control.

## Measurements

For each pair, `Delta h_i = h_var_i - h_wt_i` and `d_i = ||Delta h_i||_2`.
The audit reports mutation-site, radius-6 local mean, non-site-context mean,
and full-sequence mean delta norms, plus a distance-to-mutation curve. It also
stores the unreduced curves in `context_distance_curves.npz`.

The mutation-erasure control replaces the verified position in both WT and
variant sequences with the same `X` residue token. The program must assert
that the resulting input strings are exactly equal before inference; it then
reports every pair's maximum absolute embedding difference and fails if the
maximum exceeds `1e-5`.

## Constraints and Interpretation

- No model training, hyperparameter selection, target fitting, normalization,
  checkpoint selection, or test-label use is allowed.
- There is no context-only predictive claim in this stage. Evaluating whether a
  context-only representation predicts centered response requires a separately
  preregistered control matrix and train-parent-only fitting.
- This audit cannot authorize CIIP-1B, a deployable mutation-coordinate-free
  potential, BindingDB bridge work, or production changes.
- ESM contextual propagation, if observed, only explains why a distant window
  need not be a pure mutation-information null; it is not evidence of
  ligand-conditioned interaction.

## References Checked

- eSIG-Net (2026): residue-level mutation-site ESM features are used alongside
  interactor-conditioned discrepancy and original interaction objectives; its
  PPI task and labels are not equated with this functional kinase assay.
- ESM repository: residue-level and mean representations are distinct outputs;
  this audit uses residue-level states.
- PremPLI: mutation effects on protein-ligand binding are a structure-defined
  affinity-change task and are not conflated with the present endpoint.
- CS-DTA: cold-start claims require entity-disjoint evaluation; current split
  is pair-level only.
- Duong-Ly: the input endpoint is functional percent inhibition.
