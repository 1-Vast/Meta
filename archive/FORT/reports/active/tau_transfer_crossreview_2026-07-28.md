# Tau Transfer Cross-Review

**Date:** 2026-07-28
**Scope:** the Pair-JEPA proposal, the BCEL/physical-teacher/LLM proposal, the
original tau paper, the recent interaction-paper audit, and the current FORT
teacher-route evidence.
**Label firewall:** the schema of the affinity-bearing TRAIN registry was
inspected and only preregistered safe metadata columns were materialized. No
affinity column/value, development/confirmation value, Davis confirmation
asset, or sealed value was materialized in this audit.

## Decision

Canonical verdict:
`TAU_CONCEPT_TRANSFERABLE__PAIR_SPECIFIC_TEACHER_ABSENT__NO_TRAIN`.
The independent-review shorthand `TAU_ANALOGY_ONLY__NO_TRANSFER_TRAIN` has the
same scope and is not a separate state.

Tau supports one narrow design principle: an intervention-conditioned,
training-only delta target can increase supervision density and improve a
representation when the target is naturally co-observed. It does not supply
the missing target-specific ligand-reordering information in DTA. Neither
proposal identifies an admissible pair-specific auxiliary state on the current
substrate, so no teacher scoring, student fitting, or affinity training is
unlocked.

**Canonical stage-name mapping:** this independent review originally grouped
source/topology/teacher admission as `TAU-T0`, real teacher headroom as
`TAU-T1`, and predictive transfer as `TAU-T2`. The live `task.md` uses the
more granular names `TAU-A0` (concept), `TAU-G0` (topology), `TAU-T0`
(teacher admission), `TAU-S0` (synthetic calibration), `TAU-I0` (real
headroom), and `TAU-M0` (predictive transfer). The scientific order and stop
rules are identical.

## Established findings

1. Tau uses four robot tasks, 100 demonstrations per task, and 400 trajectories
   in total. Current tactile state plus the subsequent action sequence predicts
   the detached change in future visual features. A tactile-change-weighted
   cosine loss trains the auxiliary branch, which is removed at inference.
2. Each model was evaluated over 20 trials. Average full-task success was
   `71.25%`; the no-action-sequence, no-predictive-SSL, and no-tactile-module
   ablations were `58.75%`, `51.25%`, and `28.75%`. These results support a
   contribution from those components on the paper's four robot tasks; they do
   not establish a general mechanism or protein-ligand affinity transfer.
3. Future offsets and windows are correlated supervision from the same
   trajectory. Tau does not claim that they create independent experimental
   units. DTA ligand pairs, quartets, or teacher views likewise cannot increase
   the number of independent proteins, families, documents, sites, or
   provenance lineages.
4. The existing label-free FORT tau audit found 201,827 TRAIN registry rows,
   559 targets, and 517 homology components, but no explicit `base_protein`,
   `construct_sequence`, `directed_substitution`, or
   `directed_ligand_edit` fields. Bilateral intervention contrasts were
   therefore not identifiable. It also found zero admissible pair-specific
   teachers.
5. PBCNet2.0 is the closest relative-affinity teacher candidate, but its exact
   training membership is unavailable and its inputs require a reference
   co-crystal, Glide docking, and MCS-aligned poses. A clean
   pretraining-overlap-excluded evaluation and legal fold-local input contract
   cannot currently be constructed.
6. The physical-pose/contact route is not reopened by renaming it a future
   state. Docking, holo/native-complex, pocket, physics-increment, and
   structure-distillation framings have already failed or lack compliant
   target-reordering evidence.
7. Privileged completion also lacks demonstrated information headroom:
   `Delta_info = +0.0154`, 95% CI `[-0.0155, +0.0464]`, below the empirical
   MDE and frozen threshold `0.0452`. Correct-support specificity and public
   strict few-shot adaptation have also failed; FLOWR's positive evidence
   instead depends on dense private project SAR.
8. The current OpenMut program may recover a new source-resolved substrate,
   but it remains at source, firewall, topology, and power gates. This is an
   acquisition route, not authorization for a teacher or model.

## What Tau Transfers and Does Not Transfer

| Tau property | Transferable requirement for DTA | What does not transfer |
| --- | --- | --- |
| Naturally synchronized future visual state | A naturally observed or independently measured pair-specific auxiliary state | Docking scores, poses, SAR profiles, or LLM statements merely relabeled as a "future state" |
| Real action sequence conditions the future | A frozen, canonical intervention such as an explicit WT-to-substitution operation or directed ligand edit | An unordered ligand pair, random mutation, scaffold identity, or target ID treated as an intervention |
| Predicts a change rather than an absolute state | Target-centered ligand reordering or a mixed difference evaluated against destruction controls | Absolute potency, family, taxonomy, compatibility, pose confidence, or assay/document scale |
| Detached, training-only target | A teacher with frozen lineage and no inference-time privileged input | A teacher computed from the same affinity/SAR labels; that is supervised label transformation, not affinity-free SSL |
| Multiple future offsets per trajectory | More correlated constraints within an independent unit | More biological `n`, narrower component-level uncertainty, or greater provenance diversity |
| Auxiliary branch removed at deployment | Exact no-auxiliary and matched random-auxiliary nested controls | Evidence of strict family-cold and scaffold-cold performance without a separate prediction test |

The key non-isomorphism is causal: tau observes the same physical episode
before and after a real action. Observational DTA has no natural time axis,
action, or co-observed rich future modality. A constructed four-cell affinity
contrast can be a valid estimand, but it is an outcome, not automatically an
independent auxiliary state.

## Novelty Overlaps

| Proposed element | Existing overlap | Cross-review decision |
| --- | --- | --- |
| Pair-JEPA or BCEL four-cell difference | REWIRE, MISO/MISO-OR, DICE/AXIS, WTPAIR, KirHub DD, R-MAON, and the current OpenMut estimand | Valid algebra; not new by itself |
| Antisymmetric ligand-pair or mutation-pair head | DICE/AXIS, WTPAIR, and the direct centered R-MAON carrier | Carrier choice, not a new information source |
| Cycle consistency | Any scalar energy field already gives zero cycle sum | Diagnostic only; adds no identifying information |
| Relative-affinity teacher | PBCNet2.0 and prior pairwise/MMP ranking routes | Potentially useful only through new pair-specific teacher information; PBC is currently inadmissible |
| Physical or structural teacher | Atom-residue teachers, FIRE/BridgeFIRE, native-pose and privileged structure distillation | Closed on current evidence; requires a genuinely new measured state and clean lineage |
| Privileged completion distillation | RB-DR-QMAPD | Stopped by the failed information-headroom gate |
| Explicit observation/noise model | NARD, de-noising, STRATA, and UEL | Useful inference tool after information exists; cannot create target signal |
| Meta-adaptation or LoRA | Prior wrong-support audits and FLOWR | Public few-shot specificity failed; dense project SAR is a measurement condition |
| LLM mechanism constraints | K-LBP/LEXOR contamination and hidden-SAR concerns | Not tau transfer and not a quantitative affinity teacher |

No combination of these overlapping modules counts as novelty. A defensible
new contribution would be a separately measured, lineage-clean auxiliary state
that passes target/intervention corruption and predicts target-specific ligand
reordering across independent families and sources.

## Ordered Gates

### `TAU-T0` - auxiliary-source admissibility and natural correspondence

Affinity-free and student-free.

- Freeze source, version, checksum, license, deterministic engine or weights,
  training membership, deployment inputs, and fold-local construction.
- Require a naturally co-observed or independently measured pair-specific
  state. Affinity-derived pairs, SAR-derived profiles, docking scores, and LLM
  constraints do not qualify as self-supervision.
- On the registered auxiliary outcome, require the true target and true
  intervention to beat ligand-only projection, matched wrong target,
  within-family wrong target, intervention shuffle, and constant/random
  matched-capacity teachers.
- Count proteins, homology components, broad families, documents, sites, and
  provenance lineages before projecting a component-level MDE.

**Current status:** fail/stop. The earlier topology audit found no bilateral
intervention fields and the teacher audit admitted zero pair-specific
candidates. PBC/pose/physical candidates fail lineage, input, or prior-route
requirements.

### `TAU-T1` - frozen-teacher information-headroom gate

No student is fitted. This gate is blocked until the upstream OpenMut
source/firewall/topology/evidence/measurement gates establish a powered,
multi-family substrate. Structure- or pose-conditioned candidates remain
deferred to `OMUT-R0`, after `OMUT-M1`.

- Compare each frozen teacher separately with exact, TRAIN-only Ki or Kd
  target-centered reordering; analyze endpoints separately.
- Aggregate inference first by target/homology component and provenance.
- Beat ligand-only, matched 2D, target shuffle, within-family wrong target,
  ligand/intervention/pose shuffle, constant/random teacher, and the teacher's
  ligand-only projection.
- Require component-level LCB95 `> 0`, effect
  `>= max(0.03, empirical MDE)`, and no collapse in the worst prespecified
  family or source.

### `TAU-T2` - minimal training-only transfer

Allowed only after the same candidate passes `TAU-T1`.

- Remove the teacher at inference.
- Use an exact nested no-auxiliary null and a matched-capacity random-auxiliary
  control with identical data, basis, split, optimizer, and budget.
- Evaluate strict family-cold plus scaffold-cold validation with held
  provenance.
- Require both a prespecified predictive gain and mechanism-specific loss under
  target/intervention corruption. One favorable seed is insufficient.
- Keep a newly sealed, independent confirmation source separate from all
  teacher selection and development.

For the canonical live queue, a non-structural independently measured teacher
maps to `TAU-I0` and `TAU-M0`. Any structure-, pose-, docking-, contact-, or
complex-conditioned teacher remains blocked until `OMUT-M1` passes and may
enter only through `OMUT-R0`; tau terminology cannot move it earlier.

## Stop and Unlock Rules

- Failure at any gate closes that candidate. Do not rescue it by increasing
  offsets, views, teacher combinations, loss weight, latent width, rank,
  backbone capacity, epochs, seeds, pair expansion, or pose generation.
- `TAU-T0` does not imply affinity information, `TAU-T1` does not imply
  learnable transfer, and `TAU-T2` development success does not imply strict
  independent confirmation.
- The only current unlock is continued label-free OpenMut source recovery and
  power auditing. Real-label teacher testing remains blocked.
- A candidate may reopen only with a genuinely different admissible source:
  auditable lineage, legal non-privileged inputs, natural or independent
  pair-specific measurements, adequate multi-family/provenance units, and a
  preregistered MDE.

## Exact Evidence Paths

All hashes are SHA-256.

| Evidence | Exact path | SHA-256 |
| --- | --- | --- |
| Original tau paper | `tmp/pdfs/paper_interaction_audit_2026_07_28/tau_future_visual_2607.24485v1.pdf` | `b9dd13ee0fc69ce091c6f10b738caaca05aef95608b317e06b2c8021d4b2be4b` |
| Extracted tau text | `tmp/pdfs/paper_interaction_audit_2026_07_28/tau_future_visual_2607.24485v1.txt` | `9a9baa05d5efe896b80ee9cf81483029cfca2ee7d50bce24c02ec67de09ef3aa` |
| Pair-JEPA proposal | `C:/Users/59964/.codex/attachments/2894d3a6-c93a-4c7c-9312-71808d3cca55/pasted-text.txt` | `bd6e03265025114d32a5bde4df39d66642daf7cad3978c2861d2a8f0b8e179ad` |
| BCEL/teacher/LLM proposal | `C:/Users/59964/.codex/attachments/90eb8b94-7e66-40f4-aaaf-2d908cbb9df5/pasted-text.txt` | `69cb793a07d97d98a5f936ed18b0a0483025f2e70bfa01c0f382a0043c7c1109` |
| Recent paper audit | `reports/active/recent_interaction_paper_audit_2026-07-28.md` | `75464ac421f14909545e0622e61aa65856ebf72466e2e4aa31ea59a81999b026` |
| Recent paper gate | `reports/active/recent_interaction_paper_gate_2026-07-28.json` | `f232d55454aea0893f3df664dbcb48c9ef0e1f0f30633046aecd38d563eac51e` |
| Tau preregistration | `reports/active/tau_transfer_preregistration_2026-07-28.md` | `01afaf01cd8406f5bc8ca4cfcf148e1ffae444eadcbcb50089383cf5ec384ffa` |
| Tau feasibility result | `reports/active/tau_feasibility_a0_t0_2026-07-28.json` | `ad47f2f971f008adc50e7417bbd26b04b6ec40d1ce97c0aecfe15bf8c6dea63e` |
| Teacher manifest | `manifests/tau_teacher_admissibility.v1.json` | `2191cf76654a466ff92db6de783e42255408ffd182f6354d1be95e61070519cf` |
| OpenMut cross-review | `reports/active/openmut_delta2rank_feasibility_crossreview_2026-07-28.md` | `0f667e488aa7a60e22047d78ce72810ba5afca4f0429c06c7c7b60344e99480e` |
| Protein-conditioned signal audit | `reports/active/protein_conditioned_signal_audit_2026-07-27.md` | `8c56149da8375922b1efa153945fb6ed02ab9272f11c38cdd891d3f5f05e496b` |
| Privileged-completion result | `reports/active/rb_dr_qmapd_oracle_result.md` | `e4c6b903eb18976d64d47d414907c21ec7eb1f81399b8ac45a3f7e54e76c824d` |
| LOCK/CLOCK final audit | `reports/active/lock_clock_g0_final_audit_2026-07-28.md` | `4a0ff3f569c66d47ef2daf2a9050695c9c377de93fcd64c7ab2fa0ce8b762232` |

`task.md` and `history.md` are live coordination files and were not modified by
this cross-review; the immutable evidence artifacts above carry the
reproducible claims.
