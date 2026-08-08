# T-BASIS-R0 Fixed Radial Basis Recoverability

Registered: 2026-08-07

## Scientific Question

Can a fixed, continuous chemogeometric radial basis computed from privileged
holo coordinates be recovered from the existing sequence+2D P1B distance
distribution, and does recovery depend on the correct protein?

This is the minimum prerequisite for Universal Mechanistic Basis Distillation.
It does not test angular or many-body basis terms, affinity direction, few-shot
adaptation, production `z`, or the frozen mathematical operator.

## Forbidden Reads And Changes

- ChEMBL/BindingDB/Papyrus affinity values: forbidden.
- DAVIS/KIBA/recipient labels: forbidden.
- P1B, ESM, GINE, CSMO, Band, `K`, and theory changes: forbidden.
- Hyperparameter search, test-based model selection, and reuse of the T-DIR-P0
  panel: forbidden.

## Fresh Structural Panel

Select `192/64/64` train/validation/test complexes from the governed
`pilot20k_homology_split_v2` corpus by SHA-256 order of
`"T-BASIS-R0|source_entry_id"` before basis computation. Exclude all 40
T-DIR-P0 complexes. All 320 records must have globally distinct homology groups,
PDB IDs, exact sequences, and nonempty Murcko scaffolds. Validation/test
scaffolds must not occur in earlier selected splits. Validation/test records
must remain held out from P1B training. Protein coordinate-to-sequence mapping
coverage must be at least `0.999` to avoid treating unresolved slots as a
student error.

Create a one-to-one, same-split, score-blind wrong-protein map before basis
computation. Each wrong protein must have a different homology group, global
sequence identity `<0.40`, length ratio in `[0.5, 2.0]`, and no reuse.

## Fixed Teacher Basis

Use eight fixed, multi-hot ligand atom channels:

```text
hydrophobe, aromatic, donor, acceptor,
positive, negative, halogen, other
```

Use six fixed protein residue chemistry classes and six Gaussian radial
functions with centers `[2.0, 3.5, 5.0, 6.5, 8.0, 9.5] A`, sigma `1.0 A`, and
a cosine cutoff at `10 A`.

For canonical ligand atom `i` and P1B residue slot `s`, privileged coordinates
provide the minimum heavy-atom distance `d_is`. The fixed teacher is:

```text
Phi*[a,r,k] = (1 / N_atom) sum_i,s q[i,a] pi[s,r] R_k(d_is)
```

where `pi[s,r]` is the sequence-observable residue-class composition of slot
`s`. Coordinates are fixed and permutation invariant. They are descriptive
radial moments, not energy or uncertainty.

## Student Arms

1. `MEAN`: train teacher mean, with no pair information.
2. `P1B_RAW_CORRECT`: replace `R_k(d_is)` by its fixed expectation under the
   frozen five-bin P1B distance distribution for the correct protein.
3. `P1B_CAL_CORRECT`: apply one shared train-only `6 x 6` radial calibration
   plus intercept to the raw radial vectors. The deterministic solver is Ridge
   with `alpha=1e-3`; no model selection.
4. `P1B_CAL_DERANGED`: the same frozen calibrator with the score-blind wrong
   protein and the same ligand.
5. `TEACHER`: zero-error ceiling.

The bin-to-RBF expectation matrix is computed deterministically by uniform
quadrature within `[0,4)`, `[4,6)`, `[6,8)`, and `[8,10)`; the `>=10 A` bin is
zero because of the fixed cutoff.

## Metrics And Gate

Teacher coordinates are standardized using train-only mean and standard
deviation; coordinates with train standard deviation `<=1e-8` are excluded and
reported. Each complex is one independent unit.

For validation and test report complex-macro standardized MSE and:

```text
reconstruction_gain = (MSE_mean - MSE_correct) / MSE_mean
partner_gain        = (MSE_deranged - MSE_correct) / MSE_mean
```

Use 2,000 fixed-seed complex bootstrap replicates. No pair-level inference.

The terminal verdict is `RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED` only
if all conditions hold:

- validation reconstruction gain and partner gain are both positive;
- test reconstruction gain `>=0.10` with 95% bootstrap LCB `>0`;
- test partner gain `>=0.10` with 95% bootstrap LCB `>0`;
- all panel, mapping, hash, CUDA, and forbidden-read contracts pass.

Otherwise use one of:

```text
T_BASIS_DATA_OR_MAPPING_FAIL_CLOSED
RADIAL_BASIS_RECOVERY_NOT_IDENTIFIED
RADIAL_BASIS_PARTNER_DEPENDENCE_NOT_IDENTIFIED
```

A PASS authorizes only a separately registered angular/many-body basis
distillation study. It does not authorize real affinity, DAVIS, few-shot,
production integration, or claims of a universal interaction basis.
