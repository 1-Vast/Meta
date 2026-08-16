# Tau transfer decision

**Date:** 2026-07-28

**Verdict:** `TAU_CONCEPT_TRANSFERABLE__PAIR_SPECIFIC_TEACHER_ABSENT__NO_TRAIN`

**Final category:** 3 - current data cannot identify the mechanism; a new
source-resolved substrate or prospective measurement condition is required.

## What tau establishes

The tau paper used four robot tasks with 100 trajectories per task. Current
tactile tokens and a horizon of 32 actions predicted changes in future
pretrained visual latents. Future tactile-change magnitude weighted the
training-only auxiliary loss.

The reported mean full-task success was `71.25%`, compared with `58.75%`
without action conditioning, `51.25%` without predictive self-supervision, and
`28.75%` without the tactile encoder. Each model-task condition had 20 physical
trials, and the paper did not report task- or seed-cluster confidence
intervals. Its limited object generalization was not a dual-cold analogue.

Future frames are real additional within-trajectory sensor observations, so
they add supervision density beyond algebraically duplicating one scalar.
They remain correlated within trajectories and do not increase the number of
independent tasks, sites, environments, or provenance lineages.

## Transferable and non-transferable parts

The transferable abstraction is:

```text
current interaction state + a verifiable intervention
    -> change in interaction state
```

For DTA, this requires canonical protein and ligand interventions plus a
genuinely pair-specific post-intervention teacher. A mutation string, matched
molecular pair, masked-token target, concatenated sequence/SMILES embedding,
or expanded affinity quartet does not supply an equivalent future
observation.

The proposed bilateral affinity differences, four-cell mixed differences,
counterfactual energy, and path consistency overlap REWIRE, MISO, DICE/AXIS,
WTPAIR, and R-MAON. If all predictions come from one scalar field `s(m,l)`,
cycle consistency is algebraic and an extra cycle loss adds no information.

## Executed label-free tests

The audit ran in the `drug` environment using only these TRAIN metadata
columns:

```text
target, conn, endpoint, scaffold, assays, docs, accession, hcluster,
dual_cold_split
```

It read 201,827 rows, 559 targets, 121,401 ligand parents, 48,234 scaffolds,
23,569 assay identifiers, 9,587 document identifiers, and 517 homology
components. It found 3,363 same-target, same-endpoint groups with at least 16
ligands and an identical aggregated assay-ID bundle. This is an optimistic
metadata-density count, not a count of verified same-assay/context rectangles.

That density is not a bilateral intervention graph. The registry contains
neither the required canonical protein fields (`base_protein`,
`construct_sequence`, `directed_substitution`) nor an explicit
`directed_ligand_edit`. Therefore protein transitions, directed chemical
transitions, and complete bilateral factorial contrasts are not identifiable.

Teacher admission was then tested one candidate at a time:

| candidate | result | decisive blocker |
| --- | --- | --- |
| PBCNet2.0 | fail | unavailable exact training membership, privileged pose contract, and no auditable overlap exclusion; released weights alone are insufficient |
| physical pose/contact teacher | fail | no single frozen rights-cleared teacher; prior docking/holo/native/physics framings failed matched controls |
| LEXOR-MC | fail | MC0 extraction contract failed; no quantitative pair-specific state and no pretraining-overlap exclusion |
| separable generic embeddings | negative control only | target-only plus ligand-only projection has no interaction residual by construction |
| random and constant teachers | negative controls only | deliberately contain no semantic pair information |

All preregistered negative controls were present. Zero pair-specific teachers
passed admission. `TAU-S0`, `TAU-I0`, and `TAU-M0` are therefore blocked.

## Reopening conditions

Reopen `TAU-G0` only with exact, source-resolved WT/single-substitution
constructs and canonical directed ligand edits that recur across independent
base proteins, broad families, assays, and provenance lineages.

Reopen `TAU-T0` only with a frozen teacher version, available
weights/deterministic engine, resolved training lineage, frozen rights, legal
fold-local inputs, pretraining-overlap exclusion, pair-specific output, and
executable wrong/shuffled/separable/random controls.

If both gates pass, synthetic calibration must reuse the existing direct
regular-null R-MAON carrier. Real-label association still requires OpenMut
`I0`; predictive integration additionally requires `OMUT-C0`. A failed teacher
cannot be rescued by changing auxiliary-loss weight, offsets, views, rank,
capacity, backbone, epochs, seeds, or coordinate concatenation.

That order applies only to a non-structural, independently measured teacher.
Any structure-, pose-, docking-, contact-, or complex-conditioned teacher
remains deferred until `OMUT-M1` passes and can enter only through `OMUT-R0`.

## Artifacts

- Preregistration: `reports/active/tau_transfer_preregistration_2026-07-28.md`
- Teacher manifest: `manifests/tau_teacher_admissibility.v1.json`
- Executed result: `reports/active/tau_feasibility_a0_t0_2026-07-28.json`
- Runner: `research/tau_feasibility.py`
- Focused tests: `tests/test_tau_feasibility.py`
- Frozen paper: `tmp/pdfs/paper_interaction_audit_2026_07_28/tau_future_visual_2607.24485v1.pdf`

The audit inspected the schema of an affinity-bearing TRAIN registry and
materialized only preregistered safe metadata columns. It materialized no
affinity column/value, development/confirmation value, Davis confirmation
asset, or sealed-test value. The result flags are explicitly scoped to this
audit and do not overwrite the permanent historical firewall ledger.
