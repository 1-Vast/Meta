# DCST-R14 to R17 Complete Decision Record

Date: 2026-07-29  
Status: frozen summary of completed label-blind audits  
Scope: PLINDER-to-ChEMBL transport, SISMT, DTIOD, and KLIFS source replacement

The root-level all-dataset index is
[`DATASET_RECORD_SUMMARY.md`](../../DATASET_RECORD_SUMMARY.md). This file
retains the detailed DCST route record.

## Executive decision

The four routes below were executed without loading a ChEMBL affinity column,
without consuming confirmation labels, and without consuming the sealed test.
All routes stopped before student training, support transport, or new affinity
fitting.

| Route | Information object | Formal decision | Main failure |
| --- | --- | --- | --- |
| DCST-R14 | Absolute R6 privileged mechanism state | `STOP_PLINDER_SOURCE_EXPAND_STAGE1` | Target and ligand support are too sparse; domain separation is strong |
| SISMT/R15 | Stable source/target spectral intersection | `STOP_SISMT_LABEL_BLIND_SUPPORT_GATE` | No privileged direction survived bootstrap stability |
| DTIOD/R16 T1 | Local protein-segment by active-Morgan mixed difference | `STOP_DTIOD_T1_NO_PRIVILEGED_TANGENT` | Wrong-ligand destruction increased the response |
| KLIFS/R17 | Exact KLIFS interaction-fingerprint bridge | `STOP_KLIFS_BRIDGE_INADEQUATE` | Bridge and rectangular core are too small and provenance-sensitive |

The historical substitution-geometry conclusion remains category 3. These
audits do not show that protein information is useless. They show that the
tested information objects are not yet identified as transferable,
ligand-specific, strict dual-cold interaction mechanisms.

## 1. Frozen decision packets

### R14

- Preregistration:
  `dcst_r14_transport_support_preregistration_2026-07-28.md`
- Result: `dcst_r14_transport_support_seed1729.json`
- Decision: `dcst_r14_transport_support_decision_2026-07-29.md`

Frozen estimands and gates:

```text
delta_m(t,d) = m(t,d) - mean_anchor m(t,a)
```

- target support: at least 20% of ChEMBL TRAIN targets with maximum source
  4-mer containment at least 0.40;
- ligand support: at least 20% of audited ChEMBL TRAIN ligands with maximum
  source Morgan Tanimoto at least 0.40;
- pair responsiveness: ChEMBL TRAIN centered-moment RMS at least 20% of the
  PLINDER source value;
- strong domain separation: any five-fold out-of-fold domain AUC at least
  0.80;
- route selection was frozen before execution and used train support only.

Observed values:

| Quantity | Value |
| --- | ---: |
| exact-target firewalled PLINDER source rows | 2,106 |
| PLINDER exact targets | 767 |
| ChEMBL TRAIN targets | 559 |
| ChEMBL TRAIN ligand audit | 20,000 |
| pair responsiveness ratio | 0.946823 |
| target support >= 0.40 | 5.0089% |
| ligand support >= 0.40 | 9.8500% |
| target domain AUC | 0.908641 |
| ligand domain AUC | 0.882452 |
| privileged centered mechanism AUC | 0.941130 |

Pair responsiveness passed, but both entity-support gates failed. Therefore
the absolute PLINDER representation was stopped. Global MMD/OT, nonlinear
fusion, and source-wide importance weighting are not valid rescues.

### SISMT/R15

- Preregistration: `dcst_r15_sismt_preregistration_2026-07-29.md`
- Result: `dcst_r15_sismt_seed1729.json`
- Decision: `dcst_r15_sismt_decision_2026-07-29.md`

Frozen generalized eigensystem:

```text
Sigma_T v = lambda (Sigma_S + epsilon I) v
```

Frozen gates:

- support eigenvalue in `[0.25, 4.0]`;
- privileged mechanism overlap at least `0.05`;
- 20 target-block bootstrap median squared projection at least `0.50`;
- retained dimension between 1 and 16;
- partial transport ESS at least 88 for both source and target;
- transported mass at least 0.20;
- target and ligand coverage at least 20%;
- matched wrong-target and wrong-ligand destruction controls.

Observed values:

| Quantity | Value |
| --- | ---: |
| PLINDER source pairs | 1,613 |
| PLINDER source targets | 767 |
| ChEMBL TRAIN pairs | 4,472 |
| ChEMBL TRAIN targets | 559 |
| privileged support-compatible directions | 143 |
| directions passing initial mechanism overlap | 1 |
| surviving direction eigenvalue | 0.619565 |
| surviving mechanism overlap | 0.050474 |
| bootstrap median projection | 0.0 |
| stable retained dimension | 0 |

The only borderline direction was not stable. Partial transport and Stage-2
affinity fitting were therefore not run. Lowering the overlap or stability
threshold would be post-result selection and is prohibited.

### DTIOD/R16 T1

- Preregistration: `dcst_r16_dtiod_preregistration_2026-07-29.md`
- Result: `dcst_r16_dtiod_t1_seed1729.json`
- Complete target-block table:
  `dcst_r16_dtiod_t1_target_blocks_seed1729.csv`
- Decision: `dcst_r16_dtiod_t1_decision_2026-07-29.md`

Frozen operator:

```text
Omega(r,a) = P - P_target_mask - P_ligand_mask + P_both_mask
```

Intervention contract:

- 32 protein segments;
- four active-Morgan masks per development pair;
- protein mask replaces one segment by the mean of the other 31 raw ESM
  segments;
- ligand mask clears one active Morgan bit while preserving descriptors;
- no coordinate-dependent raw embedding Jacobian is used.

The CSV has 84 exact target blocks and 215 pairs. Each row contains:

- `target_block`;
- `pairs`;
- `PrivSemantic`;
- `NoPrivSemantic`;
- `RandomTeacherSemantic`;
- `PrivRandomSegment`;
- `PrivRandomBit`;
- `PrivWrongTarget`;
- `PrivWrongLigand`.

Aggregate target-macro medians:

| Arm | Median RMS |
| --- | ---: |
| PrivSemantic | 2.452507e-4 |
| NoPrivSemantic | 3.784427e-8 |
| RandomTeacherSemantic | 1.720295e-7 |
| PrivRandomSegment | 1.150935e-4 |
| PrivRandomBit | 1.906539e-5 |
| PrivWrongTarget | 4.873898e-5 |
| PrivWrongLigand | 2.719599e-4 |

Destruction results:

```text
D_T = 1 - R_wrong_target / R_true = 0.8012688
D_L = 1 - R_wrong_ligand / R_true = -0.1089054
```

The privileged tangent is large and target-sensitive, but it is not
ligand-specific. The wrong-ligand response is 10.89% larger than the true
pair response. T2 student distillation, T3 support testing, and Stage-2
affinity fitting were correctly skipped.

### KLIFS/R17

- Preregistration: `dcst_r17_klifs_bridge_preregistration_2026-07-29.md`
- Result: `dcst_r17_klifs_bridge_seed1729.json`
- Decision: `dcst_r17_klifs_bridge_decision_2026-07-29.md`
- Correction: `dcst_r17_target_support_correction_2026-07-29.md`

Observed values:

| Quantity | Value |
| --- | ---: |
| valid human KLIFS mechanisms | 11,250 |
| retained complexes | 5,161 |
| accessions | 181 |
| ligands | 1,681 |
| scaffolds | 1,340 |
| repeated-pair IFP Jaccard | 0.972953 |
| rectangular 2-core edges | 393 |
| rectangular 2-core connected components | 1 |
| exact ChEMBL TRAIN bridge pairs | 300 |
| exact bridge targets | 50 |
| exact bridge homology components | 45 |
| ChEMBL TRAIN ligand support >= 0.40 | 7.41% |

KLIFS improves pocket alignment and structural quality, but the bridge remains
small and is dominated by one promiscuous-ligand component. No RIFB model was
constructed.

## 2. Controls and claim boundary

The four routes used destructive or matched controls rather than relying on
one positive score:

- NoPrivileged teacher;
- uniform and random-frozen representations;
- random segment and random Morgan-bit interventions;
- wrong-target and wrong-ligand derangements;
- target, ligand, and mechanism domain classifiers;
- target-block bootstrap stability;
- rectangular target/ligand degree core;
- exact source and checkpoint hashes.

The following claims are not licensed:

- R14 proves that structural information is useless;
- R15 proves that no spectral transport is possible in any source;
- R16 proves that local protein responses do not exist;
- R17 proves that KLIFS is unusable in general;
- any route improves strict dual-cold affinity prediction;
- any failed route can be rescued by rank, capacity, epochs, or a larger PLM.

The defensible joint claim is narrower: these particular information objects
failed their preregistered label-blind transfer gates on the current source
and therefore cannot be promoted to affinity modeling.

## 3. Frozen local assets and Git boundary

Protected local assets remain outside Git:

- PLINDER and ChEMBL parquet registries;
- ESM-2 segment and accession caches;
- Morgan ligand feature caches;
- R6 checkpoint files;
- structural privileged `.npz` caches;
- all historical `.pt` model snapshots.

The Git record contains the code, tests, preregistrations, decisions, JSON
results, and the DTIOD target-block CSV. It does not contain protected data,
checkpoints, temporary downloads, or sealed labels.

## 4. Current route and reopening conditions

The DCST sequence is closed at R14-R17. The only active route is outcome-free
`OMUT-X7` / `DCST-R18`, restricted to BindingDB native-article records with
explicit PDB complex identifiers, exact single-chain sequence identity, and
non-ChEMBL provenance. A safe reader may record metadata and Ki/Kd field
presence only; numeric affinity and PDB-coordinate acquisition require a
separate authorization.

Any future DCST-like route requires a new source or information object that
passes, before affinity access:

1. source rights and immutable version checks;
2. exact target, ligand, scaffold, family, assay, document, and provenance
   firewalls;
3. nontrivial factorial target-by-ligand topology;
4. effective sample size and MDE planning;
5. a frozen mechanism-specific coordinate;
6. target and ligand destruction controls;
7. worst-family and worst-source stability;
8. a preregistered exact-null downstream operator.

No threshold may be relaxed after observing a failed route.

## 5. Verification

Executed with `D:\anaconda\envs\drug\python.exe`:

```text
622 passed, 1 warning
```

The DTIOD target-block CSV has 84 rows and 215 total pairs. Its arm medians
match the JSON artifact exactly. No new affinity label was read by R14, R15,
R16, or R17.
