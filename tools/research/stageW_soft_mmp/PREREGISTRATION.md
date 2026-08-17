# Stage W preregistration — soft controlled chemical-change surfaces on independent open datasets

Frozen **2026-08-17, before any Stage W coverage statistic was read.** Stage
U/V closed the exact-MMP route on BindingDB-Ki double-cold. This stage does
**not** continue or relax that route. It opens the redirected research goal:
first establish whether independent open datasets (Davis and KIBA separately)
provide an identifiability surface for a *soft but controlled* local
chemical-change class, then — only if the surface passes — build the local
protein × ligand representation and its training curriculum.

## 0. Redirected goal

The exact-MMP branch is closed as not estimable on BindingDB-Ki. The new
research question is unchanged in spirit but moved to surfaces that can support
it:

* Davis: independent training, independent split, independent pKd reporting.
* KIBA: independent training, independent split, independent reporting; KIBA
  scores are never compared on the pK MSE scale.
* ChEMBL Ki: optional, only after Davis/KIBA W0 passes and only with
  assay/document governance. Not part of W0.

No cross-dataset training, normalization, support, labels or checkpoints.

## 1. Data provenance

* Davis rows: `dataset/raw/dta/davis.tab` (standard DeepDTA Davis tabular
  snapshot; 30,056 rows expected). Recorded SHA-256 in the W0 artifact.
* KIBA rows: `dataset/raw/dta/kiba.tab` (standard DeepDTA KIBA tabular
  snapshot; 118,254 rows expected). Recorded SHA-256 in the W0 artifact.
* The AdaMBind public snapshot is **not used for training**: its own manifest
  records `training_authorized: false`.
* Protein-component splits are new per dataset: CD-HIT 40% on the unique
  target sequences, built with the local `tools/runtime/cdhit/4.8.1/cd-hit.exe`;
  the split seed and assignment SHA-256 are recorded. Components are the
  protein-cold unit.
* No BindingDB-Ki labels or artifacts enter this stage. `meta_test` of the
  BindingDB lineage is not opened.

## 2. W0 — soft MMP surface census and frozen admission gate

### 2.1 Construction (label-blind)

* Ligand parsing: RDKit, isomeric SMILES. Single-cut MMP fragments with
  `rdMMPA.FragmentMol` as in Stage U (`mmp.py`).
* Same-target pairs only; core = larger fragment; ties by canonical SMILES.
* **Soft family key** (structure-only, frozen):
  `sha256( murcko_core | attachment_element | attachment_aromatic |
  attachment_in_ring | category(R_a) >> category(R_b) )`
  where
  * `murcko_core` = canonical isomeric SMILES of the core after replacing the
    `[*:1]` dummy with `[H]`, passed through RDKit Murcko scaffold;
  * `category(R)` = quantized pharmacophore/change class:
    `(heavy_atoms_bucket{0,1-3,4-7,>=8}, aromatic_bool, ring_bool,
      HBD_bucket{0,1,>=2}, HBA_bucket{0,1,>=2}, charge_sign{-,0,+})`.
* Canonical direction of a family: sort `(category(R_a), R_a canonical SMILES)`
  vs `(category(R_b), R_b canonical SMILES)`; `delta_y = y(R_b ligand) -
  y(R_a ligand)`.
* Deduplication: one observation per `(target, exact core, R_a, R_b)`; the
  family row then aggregates observations with the median (expected to be
  mostly one).
* A soft family aggregates observations across different exact cores and R
  identities. It is therefore **not** an exact transformation; the W0 artifact
  must report the within-target, within-family across-exact-core residual
  (median and p95) as the chemical-control price of the softened key.

### 2.2 W0 statistics (per dataset, separately)

Report: observations, targets, protein components, exact keys, soft family
keys; family target/component degree; families spanning >=3 targets and >=3
components; connected components; top-1/top-10 family and target observation
shares; cross-component double-difference row counts for families repeated
across components; same-core residual; random-protein null for the family
statistic (permute protein identity across components, 500 draws, stable
seed).

### 2.3 Frozen W0 admission gate (each dataset independently)

1. same-target MMP observations >= **1,000**;
2. targets >= **30**;
3. protein components (CD-HIT40) >= **10**;
4. soft family keys spanning >=3 targets and >=3 components >= **20**;
5. cross-component D rows on repeated soft families >= **500**;
6. top-1 soft-family observation share <= **0.10**; top-10 <= **0.30**;
   top-1 target share <= **0.20**; top-5 <= **0.60**;
7. within-target within-family across-core residual median <= **1.00 pK**
   (Davis) / <= **1.00 KIBA score unit** (KIBA).

If a dataset fails any gate: that dataset is recorded as **insufficient for
this surface**, no W1/W2 model is trained on it, and no threshold is moved.
Both datasets failing = stop the soft-MMP route and write the negative report.

## 3. W1 (only for datasets passing W0) — local interaction representation

Deferred until W0 is read. The design is fixed in outline by the active goal:
residue/region tokens (sequence-derived), ligand atom/subgraph tokens with
explicit pharmacophore features, multiple latent pocket states, local
cross-attention interaction slots, no early pooling, independent level/shape
heads. Exact hyperparameters are frozen in a W1 amendment **before any W1
training metric is read**.

## 4. Training curriculum and stop rules (inherited from the redirecting instruction)

Only after W1 representation exists:
1. new surface insufficient -> stop, no training;
2. local representation not better than ligand-only, or shuffled protein
   reproduces the gain -> close the representation;
3. correct-protein gain only from assay/target level -> no meta-learning;
4. representation passes, then meta-learning must beat fixed Tanimoto;
5. single seed screens; promotion requires multi-seed, component bootstrap and
   independent-dataset replication.

Wrong-protein, shuffled, residue-permuted, protein-blind, capacity-matched
random branch and label-shuffle controls are mandatory in W1 evaluation.
No correct-vs-wrong hinge loss is trained (Stage S degenerate solution).

## 5. Artifacts

`W0_SOFT_MMP_CENSUS.json`, `W0_DATA_PROVENANCE.json`, `COMPONENTS_*.json`,
`tests/`, `COMMANDS.md`, `ENVIRONMENT.json`, and later W1 artifacts.
The preregistration SHA-256 is recorded in every W0 artifact.
