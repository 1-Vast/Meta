# MetaSieve L2B provenance and recovery handoff

**Date:** 2026-08-09  
**Purpose:** identify the actual sources behind the five missing L2B artifacts and provide a reproducible restart point for a local agent.  
**Governing verdict:** `FIVE_ARTIFACTS_NOT_MATERIALIZED`.

## 1. Correction to the consolidated report

The following five files named by `METASIEVE_L2B_FIVE_ARTIFACTS_CONSOLIDATED_REPORT_2026-08-09.md` are not present in the supplied workspace or Git evidence:

1. `METASIEVE_L2B_NEXT_STAGE_ANALYSIS_2026-08-09.md`
2. `L2A_ORACLE_FACTOR_BUDGET.json`
3. `L2B_INDEPENDENT_DATA_GATE.json`
4. `L2B_RESIDUE_STUDENT_FEASIBILITY.json`
5. `PREREG_L2B_PLM_LOCALIZER.md`

The consolidated report is therefore a secondary narrative reconstruction. It is not evidence that the five underlying artifacts, their scripts, checkpoints, row-level predictions, bootstrap indices, or closure assignments were ever materialized.

The five missing objects must not be recreated by copying numbers out of the consolidated report. Any replacement must be generated from newly acquired inputs under new artifact names and a newly frozen protocol.

## 2. Where the interaction-labelled data actually came from

The only identified public source for the L2/L2B residue–atom corpus is:

```text
Repository: https://github.com/lishuya17/MONN
Commit:     f2b62ccf49c18a9502aa0eb0d582c6e0735ef200
Paper:      Li et al., Cell Systems 10, 308–322 (2020)
DOI:        10.1016/j.cels.2020.03.002
License:    algorithm and data restricted to NON-COMMERCIAL use by the repository README
```

The repository documents that its benchmark was built from PDBbind v2018/PDB structures and PLIP non-covalent-interaction output. The load-bearing files are:

| Fixed-commit path | Role | Bytes | SHA-256 |
|---|---|---:|---|
| `data/out7_final_pairwise_interaction_dict` | PDBbind-v2018 development interaction dictionaries | 53,017,418 | `9e7d1128a79139cb3a43d077ba5d19cce6376ddc9cf35d65db925e7f5e7e9d82` |
| `data/independent_dataset_interaction_dict` | Additional-PDB interaction dictionaries | 4,234,610 | `377b83080190e56a5ceea09101b73234596fbf069acc4be872806c92d4d68598` |
| `data/pdbbind_all_datafile.tsv` | PDB/UniProt/ligand/sequence/assay table; **must not be read in localization-only R0** | 10,792,015 | `4234a3fe8675b733a886d8f46cf8e871a91995cf81950625bc897de36f223c30` |
| `data/independent_dataset_datafile.tsv` | Additional-PDB metadata table | 1,082,601 | `d4e0bc8e3cb16c96fbdef59d33724dac03162e6e1c431d722e29c678e1c47209` |
| `data/mol_dict` | Serialized RDKit ligand objects for the development set | 59,740,405 | `973b1edccd74dfc329e825d19fc7a31197e2c49cf9ef9237bc5276de80c89d12` |
| `data/independent_dataset_mol_dict` | Serialized RDKit ligand objects for the additional-PDB set | 2,259,203 | `fe44626723fb3bee5fa8cb26e6702158f05b5c2c9438737b0883b006eef54253` |

The MONN construction sources that explain the schema are:

- `create_dataset/Dataset_construction_protocol.txt`
- `create_dataset/step4_get_interaction.py`
- `create_dataset/step7_final_interaction.py`

`step4_get_interaction.py` gives every PLIP event an identifier of the form `<interaction class>_<event number>`. `step7_final_interaction.py` stores the same full identifier in:

```text
atom_bond_type:    (ligand_atom_name, full_event_id)
residue_bond_type: (zero_based_uniprot_residue_index, full_event_id)
```

Consequently, a positive residue–atom edge is reconstructed only by joining atom and residue records that share the **complete event identifier**. Joining only on a base class such as `Hydrogen Bonds` would create false Cartesian-product edges.

## 3. Results reproduced directly from the fixed public inputs

The included `rebuild_monn_edge_corpus.py` was run against the fixed commit above without reading either TSV affinity table. The following counts are directly reproducible from the public pickle schema:

| Quantity | Development | Additional PDB |
|---|---:|---:|
| Raw interaction dictionaries | 12,987 | 1,853 |
| Dictionaries with at least one mapped positive pair | 12,738 | 1,851 |
| Unique binary `(residue_index, atom_slot)` positives | 195,798 | 9,832 |
| Unique typed `(residue_index, atom_slot, interaction_type)` positives | 202,766 | 9,832 |
| Missing atom-name references | 0 | 0 |

This verifies the origin of the consolidated report's `12,987`, `12,738`, and `195,798` values.

It does **not** verify any split, model, metric, or confidence interval.

## 4. Claim-by-claim provenance status

| Claim in the consolidated report | Actual provenance | Current evidence status |
|---|---|---|
| 12,987 raw development dictionaries | MONN `out7_final_pairwise_interaction_dict` | **Directly reproduced** |
| 12,738 mapped development complexes | Full-event join in the same MONN pickle | **Directly reproduced** |
| 195,798 binary positive edges | Full-event join and pair deduplication | **Directly reproduced** |
| 8,646 anchor complexes | Missing custom giant-component/closure script | **Unverified external claim** |
| 4,067 satellites / 701 components | Missing exact-graph, homology and component assignment artifacts | **Unverified external claim** |
| Full pair AP 0.01786 | Missing B4/L2 checkpoint and per-complex predictions | **Unverified external claim** |
| Oracle residue AP 0.25262 and all oracle intervals | Missing oracle program, predictions and bootstrap indices | **Unverified external claim** |
| 72,226 same-protein/different-ligand pairs | Missing exact-ligand-graph corpus and analysis program | **Unverified external claim** |
| Student AP 0.08664 vs 0.06798 | Missing feature generator, fitted student, B4 predictions and bootstrap indices | **Unverified external claim** |
| 524 complexes / 157 components / 2,555 edges | Missing cross-corpus closure program and retained-row manifest | **Unverified external claim** |
| ESM2 preregistration and thresholds | Only summarized inside the consolidated report; original preregistration absent | **Not formally frozen** |
| `esm2_t33_650M_UR50D` was not executed | No PLM checkpoint or prediction artifact is present | **Supported by absence, not a performance result** |

All unverified numerical claims must retain external-claim status. They are not asserted false; they are inadmissible as baselines until reproduced.

## 5. Why the supplied theory/model files are not the source of the five artifacts

The supplied theory files define constraints and boundaries:

- the sole production operator is `A(F,z) = K(B(z)F(z))`;
- pairwise/listwise/metric ranking is not claimed;
- a biological statistic cannot enter `z` before identification;
- support-derived information must respect section/identifiability constraints.

The supplied implementation contains a low-rank atom–residue contact/distance bridge. That bridge predicts contact and distance-bin geometry. It is not the missing exact-residue PLIP localizer checkpoint used as the B4 incumbent.

The supplied README also states that there is no assembled production DTA pipeline and that no biological statistic has passed admission. Thus these files constrain interpretation but cannot generate the reported L2A/L2B AP values.

## 6. Independent blockers at R0

### 6.1 Missing labelled corpus and derivation code

This part is recoverable from the fixed MONN source, subject to its non-commercial-use condition. The included script recovers only the raw interaction corpus. Exact ligand graph identities, homology closure, component assignment and all statistical results still require new implementation.

### 6.2 Missing B4 incumbent

The comparison `B5 - B4 >= 0.02 AP` is undefined because no B4 exact-residue localizer checkpoint, prediction table, or preregistration is present. P1B `best.pt` cannot substitute: it targets contact/distance geometry over a different representation and estimand.

Before B5/ESM2 is registered, a new B4 must be independently preregistered, trained, frozen and evaluated. The B4 artifact bundle must contain at least:

- model source and configuration;
- checkpoint plus SHA-256;
- training-row manifest;
- negative-sampling manifest;
- exact development/satellite component assignments;
- per-complex complete-matrix predictions;
- metric implementation and deterministic tie policy;
- bootstrap seed and sampled component indices.

## 7. Correct restart sequence for a local agent

Use new stage and metric identities. Do not reuse `4,067/701`, `B4`, or `B5` as if the missing run existed.

1. **Authorization and license gate.** Record that MONN data will be used for non-commercial research and obtain project authorization to acquire it.
2. **Acquire and pin source.** Clone the MONN repository and check out exactly `f2b62ccf49c18a9502aa0eb0d582c6e0735ef200`.
3. **Raw-corpus gate.** Run `rebuild_monn_edge_corpus.py --strict-hashes`; require exact agreement with the six reproduced counts above.
4. **Ligand identity gate.** In a separately pinned RDKit environment, deserialize the two molecule dictionaries and generate canonical exact-graph hashes. Record RDKit version, sanitization policy and failures.
5. **Closure preregistration.** Freeze exact PDB, exact sequence, UniProt, exact ligand graph and 40% homology rules. If CD-HIT-2D is unavailable, either install it or preregister a replacement and accept a new split identity. The old report's command was `-c 0.4 -n 2 -G 0 -aS 0.8 -g 1`, but the absent orchestration script means the old closure cannot be assumed reproduced.
6. **Build and freeze B4.** Register a matched non-PLM exact-residue baseline before looking at B5 improvement.
7. **Register B5.** Only after B4 exists, freeze the PLM, head, controls, thresholds and evaluation protocol.
8. **Run development.** Keep any independent confirmation cohort sealed until development passes.
9. **Register downstream affinity testing separately.** Localization does not authorize affinity or few-shot support adaptation.

Minimal acquisition commands:

```bash
git clone https://github.com/lishuya17/MONN.git external/MONN
git -C external/MONN checkout f2b62ccf49c18a9502aa0eb0d582c6e0735ef200

python research/Residue_locator/R0_PROVENANCE_RECOVERY/rebuild_monn_edge_corpus.py \
  --monn-root external/MONN \
  --output-dir report/residue_locator/r0_raw_corpus \
  --strict-hashes
```

Adapt the script path to the repository location in which this recovery package is placed.

## 8. Formal state after this recovery

```text
PUBLIC_MONN_SOURCE: IDENTIFIED_AND_HASHED
RAW_EVENT_JOIN: REPRODUCED
RAW_DEVELOPMENT_COUNTS: REPRODUCED
FIVE_ORIGINAL_ARTIFACTS: NOT_MATERIALIZED
OLD_DERIVED_SPLITS_AND_METRICS: UNVERIFIED_EXTERNAL_CLAIMS
B4_INCUMBENT: ABSENT
B5_PLM_GATE: NOT_REGISTERABLE_UNTIL_B4_EXISTS
AFFINITY_TRANSFER: NOT_AUTHORIZED
K_LE_5_SUPPORT_ADAPTATION: FORBIDDEN
FROZEN_PRODUCTION_OPERATOR: UNCHANGED
```

## 9. Primary sources

- MONN repository: <https://github.com/lishuya17/MONN>
- MONN paper: <https://doi.org/10.1016/j.cels.2020.03.002>
- PLIP repository: <https://github.com/pharmai/plip>
- ESM repository and model names: <https://github.com/facebookresearch/esm>
- RCSB PDB Data API: <https://data.rcsb.org/>

These sources establish data and software provenance. They do not validate the missing custom splits, models or results.

