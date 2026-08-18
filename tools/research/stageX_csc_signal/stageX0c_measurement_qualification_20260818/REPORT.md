# Stage X0c report — measurement-pipeline qualification (corrected successor)

Frozen preregistration: PREREGISTRATION.md, SHA-256
7de23c8131860ca4426e12c4e88de2b5453f47ca5b4d7b22754226e6309922cd (verified).
Governance: the original Stage X0 is ruled INVALID INSTRUMENT
(../X0_INVALID_INSTRUMENT_VERDICT.md); its distance-ratio capability gate is a
measurement-definition failure, and no original artifact or threshold was
modified.

## Gate summary

| Gate | Status | Evidence |
|---|---|---|
| Q0-A ProteinGym external validation | PASS | Q0A_PROTEINGYM_VALIDATION.json |
| Q0-B historical/construct mapping | PASS | Q0B_MAPPING_AUDIT.json, Q0B_ALIAS_LEDGER.md |
| Q1 representation capability (probe selectivity) | PASS | Q1_SELECTIVITY.json |
| Q2 fully synthetic planted harness | FAIL | Q2_PLANTED.json |
| Q3 Saifudeen panel qualification | PASS (census delivered) | Q3_SAIFUDEEN_CENSUS.json |
| I6 production-dataflow integrity | PASS (23 tests) | tests/test_x0c_integrity.py |

## Q0 — variant-coordinate layer

### Q0-A external validation (ProteinGym)

PASS — 45,623 sampled records across 217 DMS assays (SHA-256-seeded
sample from the official ProteinGym v1.3 zip, SHA-256
3a83766254ac9ac9984ec25cb73c6e010ea4418f5e35f143933e6b6e6473b921, reference
sequences from the official reference file): old-residue agreement 100%,
mutated-sequence agreement 100% against the frozen threshold of 99.5%.

### Q0-B historical/construct mapping

BRAF V599E -> V600E is an explicit, cited historical alias with sequence
evidence: the 1992 reference M95712.1 CDS is 2298 nt (765 aa) and lacks
exactly 3 nt (one codon) relative to the current NM_004333.4 CDS (2301 nt,
766 aa; its translation equals UniProt P15056). The alias is NOT generalized
to any other protein. PDGFRalpha rows use human P16234 (S1's Q9DE49 is Danio
rerio pdgfra); D842V is quarantined because its reported construct range
668-1210 exceeds the canonical length (1089) and is never silently repaired.
All 76 Duong-Ly variant rows are typed: 65 admitted single-point pairs,
11 quarantined with reasons. KLIFS pocket numbering validated: the gatekeeper
maps to pocket index 45 for 10 known gatekeeper mutations. Davis construct
census recorded (mutant/phosphorylation-state semantics).

## Q1 — representation capability (probe + control-task selectivity)

Frozen gate: at least one representation outside {edit_descriptor, random}
with selectivity >= 0.10 and cluster-bootstrap 2.5% lower bound > 0 on task
T-A (pocket membership) under leave-one-parent-out. Result: PASS.

Passing representations: pair_centered_local_esm (selectivity +0.189,
CI [0.033, 0.363]), mutation_position_only (+0.110, CI [0.021, 0.230]),
substitution_type_only (+0.209, CI [0.007, 0.420]). Not significant:
residue_identity_context (+0.054), local_onehot_window (+0.027),
klifs_pocket (-0.086), and the ID/composition/global representations
(0.000). Global pooled ESM and composition read no pocket-membership
information; the ESM local-window representation reads it.

Capacity curve (MLP-8, same LOO-parent split): pair_centered_local_esm
0.919 vs linear 0.954 - a larger probe does not inflate the result.
Random-label curve: linear probes on globally shuffled labels reach 0.451 /
0.446 / 0.644 for the three passing representations, i.e. no gain from
label structure for the ESM and position representations (the type
representation's 0.644 is small-sample overfit; its real-label selectivity
0.209 remains the admissible claim).

Interpretation: Q1 certifies readable information only, not downstream
causal use. Distance ratios remain diagnostics.

## Q2 — fully synthetic planted-signal harness

Gate point (tau*=1.0, rank 4, dense, eval cells, median of 3 seeds):
correct-arm Spearman 0.033 (need >= 0.30), dead-zone sign accuracy 0.504
(need >= 0.70), gap vs ligand_only -0.018 (need >= 0.05) -> FAIL.
Detection ladder tau* {0.125,0.25,0.5,1.0,2.0} x rank {1,4,16} is in
Q2_PLANTED.json. Negative controls behaved as designed: label permutation
dead-zone sign accuracy 0.52-0.57 (chance), tau*=0 interaction head
recovers nothing, floor-clamp imputation induces the expected spurious
recovery (dz 0.588), no-interaction head and shuffled/family-shuffled/
random protein ~0.50-0.55, free-target-id upper bound ~0.50-0.55.
Diagnosis (separate frozen-seed runs, recorded in the artifact): the
oracle arm (P@U input) recovers the centred interaction at dz 0.68-0.76,
so the information exists; the correct arm (one-hot pocket input) tops
out at dz ~0.58 across 8+6 restarts and two training protocols, so the
failure is representation-learning/optimization capacity at this sample
size, not information absence.

## Endpoint / censoring ladder (Duong-Ly + Saifudeen)

Duong-Ly: continuous % remaining activity with observed range -12.5..191.3;
4,023/17,934 cells > 100; no hard floor. Ladder used throughout: (1) logit
transform with interval censoring when the rounded endpoint hits 0/100;
(2) raw % remaining; (3) sign-only comparisons for dead-zone evidence.
Saifudeen: 21.3% of values exactly 100, 0.5% exactly 0; only 103/349
variants have >= 25% of inhibitor values strictly inside (20, 80), so any
B1-stage reuse must model per-row censoring, never impute at the floor.

## Q3 — Saifudeen 2026 panel qualification

Verified against the paper and its supplementary file: 92 inhibitors (86
approved of ~100), 409 wild-type activity columns, 349 variant columns,
duplicate measurements at 1 uM, per-kinase Km ATP (supplementary methods),
KIRHub portal kirhub.fredhutch.org, license CC BY-NC-ND 4.0.

Pairability census (Q3_SAIFUDEEN_CENSUS.json): 313/349 variants have a
matched WT gene; 272 have a matched WT gene AND substrate; 34 have no
matched WT; 21.3% of the 27,088 activity values are exactly 100 (saturated),
only 103/349 variants have a responsive window (>= 25% of inhibitor values
strictly inside 20-80%). Construct-background equality between WT and
variant rows is NOT assumed (full-length vs domain fragments differ); it
remains an unresolved per-pair dimension. The panel is a functional
inhibition positive control, never pK/pIC50/DTA.

## I6 — production-dataflow integrity

23 contract tests pass against the production objects (csc.py, q2.py
generator/training, x0_common): antisymmetry, identity-pair zero,
reference-term sign flip, train-only reference statistics, stable seeds
across processes, no Python hash(), gradient coverage for every trainable
branch, dead-branch capture, permutation controls destroy target
information, matched arms share cells/masks, no parent/scaffold crosses
blocks, unique cells, interval-bound directions, sign-only target direction,
planted truth bitwise recomputable, cluster bootstrap resamples clusters,
restricted data not committed, old-residue hard rule, BRAF alias not
generalized.

## Stage result

Q0-A, Q0-B, Q1, Q3 and I6 pass; Q2 FAILED its frozen gate. Therefore
B1/B2/C/D are NOT authorized and no real-data biological inference may be
drawn from this harness. The Q2 failure is optimization-limited (the
oracle arm recovers the planted signal at dz 0.68-0.76, the correct arm
does not), which defines the single highest-information next step: a new
preregistration for the representation-learning fix or a revised gate
point — the current frozen gate may not be moved retroactively.
