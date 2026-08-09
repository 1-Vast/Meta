# MetaSieve experimental evidence ledger

Updated: 2026-08-09.

This is the canonical human-readable summary of completed experiments. It
does not replace immutable manifests, raw JSON results, or the full historical
ledger in `history.md`.

## Evidence precedence

When records disagree, use this order:

1. immutable input/output manifests and raw machine-readable Gate artifacts;
2. `project_state.json` for current authorization and freeze state;
3. this ledger for the current interpretation of completed experiments;
4. `history.md` for the complete chronological failure and decision record;
5. prose supplied from outside the repository, which is not evidence until its
   commits, manifests, predictions and label-access audit are recovered.

Historical test counts and status statements describe the repository at the
time of each experiment. The current consolidated regression count is 75.

## Current conclusion

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_IDENTIFIED
PAIR_COMPATIBILITY_IDENTIFIED
AGGREGATE_ESM_ECFP_PROBE_NOT_PROTEIN_SPECIFIC
PAIR_LOCAL_P1B_MECHANISM_OBSERVABILITY_NOT_TESTED
POSE_FREE_CLASS_NOT_CLOSED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

The earliest unresolved boundary is not whether protein information reaches
the geometry bridge. P1B established that it does. The unresolved question is
whether the actual frozen P1B atom-local, residue-local and pair-local state
contains a correct-protein mechanism that can later add affinity value beyond
ligand-only and wrong-protein controls.

## Retained PASS evidence

| Stage | Verdict | Verified result | Scope |
|---|---|---|---|
| P0 | PASS | Canonical DTA, sealing, label-firewall and frozen-operator contracts pass regression. | Software/data contract only. |
| P1A | PASS | 14,906 governed holo complexes, 14,906 receptor sequences and 2,869 chemotypes after 739 protected-homology exclusions. | Open structural corpus; no affinity claim. |
| P1B | `PASS_GEOMETRY_IDENTIFIED` | On 1,477 controlled test complexes, correct contact AUPRC 0.43885 versus wrong-protein 0.05149 and wrong-ligand 0.23895. Correct distance MAE 1.97541 A versus wrong-protein 2.66531 A. | Correct-protein contact/distance geometry only. |
| D0-C | PASS | Official ChEMBL37 SQLite archive verified by SHA-256; 343,562 canonical Ki/Kd rows and 41,619 tasks constructed deterministically. | Immutable source corpus provenance; no model trained. |
| D1 | PASS | 3,817 governed tasks, 697 targets and 253 homology/document closure components; fixed fold sizes 1467/588/588/587/587 with zero closure crossing. | Independence governance; no affinity identification. |
| E0S/E0R | Synthetic diagnostic PASS | Teacher statistic was reconstructable from frozen geometry (maximum error 2.19e-7). Objective mismatch was found; a Moore-Penrose witness gave train RMSE 3.18e-8 and holdout RMSE 0.01269, with correct CI 0.99737. | Synthetic realization and numerical diagnosis only. |
| T-BASIS-R0 | Structural PASS | Fixed 288D radial chemistry basis was held-out recoverable and partner-dependent: reconstruction gain 0.5312 [0.4433, 0.5962], partner gain 0.1561 [0.1070, 0.2007]; all 288 coordinates active. | Structural statistic only; affinity reads zero. |

The retained P1B raw Gate is
`report/mechanism_refactor/p1b_gate_pilot20k_seed17_v4/gate_report.json`.
The retained D0/D1 report is
`report/mechanism_refactor/p1r2b_d0_chembl37_v1/STAGE_REPORT.md`.
Later synthetic and research implementations were removed after consolidation;
their results remain in `history.md` and their source trees are recoverable from
Git history.

## Mechanism-to-affinity experiments

| Stage | Main result | Verdict and interpretation |
|---|---|---|
| P1C | Ligand CI 0.71110; correct mechanism 0.60629; deranged 0.60715. Correct-minus-ligand -0.10481 and correct-minus-deranged -0.00086. | FAIL. The legacy readout had no usable correct-protein affinity increment. |
| P1R0 | Legacy statistic was atom-order dependent; PCA32 retained 48.78% of correct-vs-deranged energy. | Readout contract defect confirmed; PCA compression alone did not explain failure. |
| P1R1 | Invariant 288D MIF: ligand 0.71110, correct 0.71874, deranged 0.69516. Correct-minus-ligand +0.00764 [-0.00179, 0.01685]; correct-minus-deranged +0.02359 [0.00990, 0.03761]. | FAIL below frozen +0.03 Gate. Protein contrast recovered, but not ligand-baseline affinity value. |
| P1R2A | Variance decomposition: ligand 77.58%, protein 20.76%, non-additive interaction 1.67%. Correct-minus-ligand -0.00596; correct-minus-deranged +0.02373. | FAIL. The interaction residual remained partner-sensitive but was not affinity-incremental. |
| P1R2B0 | Source OOF correct-minus-ligand: Ridge -0.00263, spline -0.01980, MLP -0.03026. Metaval Ridge +0.00441; nonlinear arms remained non-positive. | FAIL. Greater nonlinear capacity after global pooling did not rescue affinity semantics. |
| P1R2B1 | Strongest metaval MLP correct-minus-deranged +0.03895 [0.01791, 0.06010], but correct-minus-ligand +0.01915 [-0.00296, 0.04096]. | FAIL. Compatibility/wrong-protein penalty without a stable positive correct-protein increment. |
| E-AFF-P0/H0A/H0C | No population-shared radial affinity direction; task-local headroom was ligand/series SAR; centered radial interaction residual did not recover partner affinity. | Negative evidence against the tested radial affinity mappings, not against all pair-local biology. |

These experiments establish the recurring failure mode:

```text
correct protein changes geometry or compatibility
    !=
correct protein adds transferable affinity-ranking information
```

## Data and estimand feasibility

| Stage | Status | Meaning |
|---|---|---|
| F0R | Historical failure closed | Live ChEMBL API rehydration could not reproduce legacy JSONL bytes. This does not invalidate ChEMBL37; D0 replaced the live API with a release-pinned dump. |
| E-AFF-X0 | `STOP_SOURCE_INTERACTION_UNDERDETERMINED` | The frozen crossed-rectangle independence requirement was unattainable in the available source design. This is a data/estimand limitation, not a model failure. |
| E-AFF-L0/L0R | NOT RUN scientifically | The positive control failed. The 195-task/78-component audit therefore did not test protein-specific affinity location. |
| XP1/XP2 | Development evidence only | Consumed kinase panels contained interaction signal, but did not meet the `k<=5`, double-held-out and fresh external-admission requirements. |
| XP3/XP4/XP5 | Public-data boundary | Low-noise panels had too few independent protein groups; broader BindingDB panels had estimated assay noise (0.7774) above interaction SD (0.4058), and the tested radial basis did not generalize. |

A blocked or underdetermined estimand is not recorded as a failed biological
model. New data may reopen it only under a separately frozen acquisition,
independence and power contract.

## Structural self-supervision evidence

S0-S4 used 1,118 RCSB complexes with 621 protein clusters and 586 ligand
scaffolds, disjoint from the 10,468 P1B-exposed PDB IDs. The deterministic
six-channel 3D teacher passed rotation, translation, atom-permutation and
determinism checks at machine precision.

The aggregate mean-ESM plus ECFP Ridge probe recovered two teacher totals:

| Channel | R2 versus mean | R2 versus random | Correct minus deranged |
|---|---:|---:|---:|
| Directional H-bond | +0.268 [0.166, 0.378] | +0.366 [0.222, 0.505] | +0.037 [-0.015, 0.084] |
| Hydrophobic burial | +0.299 [0.162, 0.454] | +0.307 [0.167, 0.431] | -0.006 [-0.035, 0.024] |

Attribution localized this aggregate signal to the ligand. For H-bond totals,
ligand-only R2 was 0.266, protein-only 0.009 and joint-minus-ligand 0.0015. For
hydrophobic burial, ligand-only was 0.331, protein-only -0.017 and
joint-minus-ligand -0.032.

The valid conclusion is therefore
`AGGREGATE_ESM_ECFP_PROBE_NOT_PROTEIN_SPECIFIC`. S4 did not evaluate P1B's
atom-local, residue-local or atom-by-slot tensors and is not an upper bound on
the whole sequence-plus-2D class. No GPU training was authorized by this
aggregate result.

## Claims not admitted as evidence

The externally supplied S5-S9 SSL report has no matching commits, immutable
manifests, checkpoints, prediction files or result JSON locally or on the
audited remote branches. Its numerical claims remain
`EXTERNAL_CLAIM_NOT_REPRODUCED`. Independent code-contract findings from that
report were separately audited and fixed; those fixes do not validate its
experimental claims.

## Failure taxonomy

The 94 historical failure entries in `history.md` are grouped as follows:

| Range | Theme | Durable lesson |
|---|---|---|
| F-01--F-19 | Early meta-architecture and representation attempts | More attention, deeper encoders and generic adaptation did not establish support or correct-protein causality. |
| F-20--F-60 | Phase-Z data, context, acquisition and recovery | Provenance, endpoint semantics, target independence and immutable releases must be treated as model prerequisites. |
| F-61--F-85 | P1 mechanism and E0 realization | Geometry is identifiable; global summaries and tested affinity readouts are not sufficient. Synthetic objective defects must not be confused with biological failure. |
| F-86--F-94 | Affinity estimand and public-data feasibility | Some questions are underdetermined by available independent panels; a blocked Gate is not permission to relax independence. |

## Active authorization

Only `P1R2B-S5_LOCAL_MECHANISM_OBSERVABILITY` is active. It must test the
actual frozen P1B local contract in this order:

1. CCD atom and canonical residue mapping;
2. single-chain versus oligomer/interface deployment audit;
3. chemistry-faithful deterministic pseudo-teacher;
4. exact-residue to 128-slot retention ceiling;
5. ligand-only, random, deranged-protein and pair-shuffle controls;
6. synthetic trainability control;
7. conditional lightweight pair-local GPU distillation.

The following remain frozen:

- real ChEMBL/BindingDB affinity training;
- DAVIS and recipient labels;
- typed-interaction production integration;
- production biological `z`;
- CSMO, Band, mesh and frozen theory;
- P2-P4.

Even a structural S5 PASS would only authorize a separately registered source
affinity Gate. Production admission requires closure-component OOF
`correct-ligand >= 0.03` and `correct-deranged >= 0.03`, both with 95% lower
confidence bounds above zero, followed by a sealed transfer Gate.

## Recovery and reproducibility

Terminal-negative implementations and duplicate reports were removed after
their evidence was recorded. They remain recoverable from commits `3281780`,
`12a2765`, and `608decf`. Large releases, embedding banks and caches are not
redistributed; see `DATA_AVAILABILITY.md`.

For the current executable state run:

```powershell
conda run -n drug python -m pytest -q
```

