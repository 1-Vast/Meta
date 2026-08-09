# Preregistration — P1R2B-PHASE2B-S0R

## Label-blind complete-panel replay of the synthetic objective audit

Stage identifier: `P1R2B-PHASE2B-S0R_COMPLETE_PANEL_REPLAY`

Written: 2026-08-10. Base repository commit: `fd30250`. This document is
frozen before any S0R teacher, trajectory, AP, BCE or bootstrap result exists.

## 1. Why S0 is not a contract-level result

The prior S0 result remains useful development evidence, but three violations
make its formal earliest verdict `SYNTHETIC_CONTRACT_INVALID`:

1. the decisive trajectory used `rows[:4000]`, representing only 2 closure
   components, instead of the registered complete held-out panel;
2. its pair universe was inherited from `build_pairs(..., masks)` and therefore
   depended on real MONN residue-edge symmetric differences despite the frozen
   synthetic-only label firewall;
3. canonical checkpoints, optimizer states and prediction tables behind every
   AP were not materialised and reloaded.

Consequently `SYNTHETIC_CONTROL_LOSS_MISALIGNED` is downgraded to a
two-component development diagnostic. It is not erased and its numbers are not
changed. The unexecuted S1 control repair (`4850c7d5...`) is superseded before
execution because it inherits the invalid panel and scales only one factor.

## 2. Scope

S0R is synthetic-only. It reads no real residue-edge label, affinity value,
ChEMBL, BindingDB, DAVIS, KIBA, recipient or metaval value. It does not train
the real Phase 2B candidate, modify `model/` or `scripts/`, add architecture,
change rank, lower `AP_bidir >= 0.50`, or touch

```text
A(F, z) = K(B(z)F(z)).
```

Seed `20260905` is used only to replay the historical development control and
is permanently forbidden as confirmation.

## 3. Frozen metadata-only panel

The census was performed before this document without constructing a teacher
or computing any prediction metric. The runner may read only the following
whitelisted fields:

```text
source_key, pdb_id, uniprot_id, uniprot_sequence, ligand_ccd,
derived seq_key, graph_key, Murcko scaffold, n_atoms, n_res, component
```

Pair eligibility is exact-sequence equality, distinct exact ligand graphs and
distinct non-empty Murcko scaffolds. It never uses `edges`,
`positive_binary_edges`, a residue mask, mask Jaccard or symmetric difference.

Frozen census:

| panel | pairs | constructs | protein closure components |
|---|---:|---:|---:|
| train universe | 228,845 | 797 | 587 |
| train hash panel | 14,333 | 367 | 298 |
| held-out A complete | 44,746 | 141 | 81 |

The train hash rule remains
`int(sha256("a|b")[:16], 16) % 16 == 0`. Component overlap across train and
held-out A is zero.

Frozen artifacts:

```text
metadata_only_records.jsonl
  dcc2a0cd0cd958640f6548ed7cd3e5076e6c0ae9f6e455e35ad3e82f534dc856
train_pairs.jsonl
  847a0770856e7e26903e56bc40d7249bb4bc21082392aa7edf1e3166014c1195
heldoutA_pairs.jsonl
  29925ed5139b1d054aa3b2a26d5a0b281336e717fc801dcfae553b5e6cc340ae
```

Any count or hash mismatch is `SYNTHETIC_CONTRACT_INVALID`.

## 4. Mandatory execution artifacts

Before training, materialise a deterministic batch stream and freeze its SHA.
Training must reload that file, verify the SHA and reject unknown IDs. At
updates `0, 1, 10, 100, 210`, save and hash:

* model checkpoint;
* optimizer state;
* complete held-out prediction table sufficient to recompute every pair AP;
* pair-, construct- and component-level metric tables.

Reload every checkpoint and prediction table and reproduce the reported AP and
BCE. The program, not a later manual script, writes the terminal verdict and
all `NOT_RUN` records before exit.

## 5. Candidate-path and label-firewall checks

The teacher copied through the production `Head` must satisfy the existing
tolerances: relative field error `<=1e-4`, AP difference `<=1e-3`, relative
product error `<=1e-4`, exact antisymmetry and identical-ligand zero. A sentinel
record whose `edges` accessor raises must pass the panel builder without raising.

Failure is `SYNTHETIC_CONTRACT_INVALID`, stop.

## 6. Full-panel replay

Run two otherwise identical trajectories under the original BCE, AdamW,
`lr=1e-3`, `weight_decay=1e-4`, gradient clip 5.0, sampler seed 20260902 and
210 updates.

### T1 — original teacher gauge

Initialise exactly at `(U*, V*)`.

### T2 — balanced ray-optimum gauge

Estimate the positive ray scale `a*` using only the complete 14,333-pair train
panel. Initialise

```text
U0 = sqrt(a*) U*
V0 = sqrt(a*) V*
```

and verify `U0^T V0 = a* U*^T V*`. Scaling only one factor is prohibited.
Report the factor-balance residual and product equality.

For both trajectories, report complete-train component-macro BCE and complete
held-out component-macro AP/BCE at every checkpoint. No prefix or subset is a
decision panel.

## 7. Component inference and decision rules

For each held-out component, first average pair AP inside construct and then
construct AP inside component. Bootstrap the 81 closure components with 2,000
replicates and seed 20260903. For `DeltaAP = AP_100 - AP_0`, report the point
estimate and one-sided 95% upper confidence bound.

A trajectory is `MISALIGNED` iff all hold:

```text
complete_train_BCE_100 < complete_train_BCE_0
heldout_component_macro_AP_100 <= AP_0 - 0.05
UCB95(DeltaAP) < 0
```

The old 2-component subset is never a decision unit.

## 8. Earliest terminal verdict

Exactly one:

```text
SYNTHETIC_CONTRACT_INVALID
  any label-firewall, count, hash, artifact, reload, determinism or numerical
  contract fails

SUBSET_SELECTION_ARTIFACT
  neither complete-panel trajectory is MISALIGNED

SCALE_PARAMETERIZATION_MISMATCH
  T1 is MISALIGNED but balanced T2 loses at most 0.05 AP or has UCB95 >= 0

SURROGATE_AP_MISALIGNMENT_FULL_PANEL
  balanced T2 is MISALIGNED

OPTIMIZER_OR_IMPLEMENTATION_DEFECT
  product equality, finite optimization or bit-exact replay fails after the
  contract checks
```

No budget scaling, margin filtering, continuous-field witness, full-W witness,
fresh teacher seed or student score is permitted in S0R. Those are authorized
only by a separately frozen next-stage contract selected by this verdict.

## 9. Required outputs

```text
S0R_INPUT_AND_FIREWALL_MANIFEST.json
S0R_FROZEN_STREAM_MANIFEST.json
S0R_CANDIDATE_PATH_WITNESS.json
S0R_ORIGINAL_GAUGE_TRAJECTORY.json
S0R_BALANCED_GAUGE_TRAJECTORY.json
S0R_COMPLETE_PANEL_COMPONENT_INFERENCE.json
S0R_ARTIFACT_RELOAD_AUDIT.json
PHASE2B_S0R_VERDICT.json
PHASE2B_S0R_COMPLETE_PANEL_REPORT.md
```

All code, inputs, checkpoints, prediction tables and outputs receive SHA-256
entries in a detached manifest. Self-hashes inside the file being hashed are
prohibited.

## 10. Frozen boundaries

No result from S0R identifies biology, ligand-conditioned localisation,
affinity, selectivity, energy, a few-shot section or a valid `z`. Real Phase 2B,
independent structural confirmation, source affinity, meta-adaptation and the
law operator all remain frozen.
