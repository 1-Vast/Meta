# OPEN_DATA_ONLY amendment

Adopted 2026-07-25. Binds every stage of the FORT dual-cold DTA program from this point.

## 1. Admissibility rule

A dataset may be used for training, representation learning, model selection, external evaluation or
confirmation only if all of the following hold:

1. it downloads from a stable public repository without a private agreement, NDA or author-only
   access;
2. its licence or public data-use terms are recorded in `manifests/open_sources.json`;
3. raw measurements, compound structures, target identifiers and assay metadata are recoverable;
4. the exact version, download URL, download date and SHA-256 are recorded;
5. another researcher can reconstruct the same registry from public materials;
6. no proprietary API, commercial database, private supplement or undisclosed internal dataset is
   required.

**A publicly readable paper whose underlying data is not downloadable does not satisfy this rule.**

Open availability never establishes statistical independence. Every candidate is additionally
audited for overlap at publication/document, assay, target construct, protein-homology component,
compound connectivity, scaffold and measurement-cell level. Two database names are not two
independent sources.

## 2. Dataset-role separation

Each dataset is assigned exactly one role *before* its labels are used, and a role is never changed
after results are observed:

| dataset | role | status |
|---|---|---|
| ChEMBL-37 local extract (CC BY-SA 3.0) | train-only mathematical development | active |
| `panel_metz` (ChEMBL doc `CHEMBL1201862`) | train-only development; its development rows are **spent** | spent by Gates PA/PB/PC |
| `panel_davis` (ChEMBL doc `CHEMBL1908390`) | single-use external panel | registered, sealed, **underpowered**, unconsumed |
| BindingDB native articles 202607 (CC BY 3.0 US) | candidate external panel | **rejected on power** (section 4) |
| Klaeger / PXD005336 `Kdapp` | candidate external panel | **not openly recoverable** (section 3) |
| PKIS2 / KCGS percent-inhibition | representation pretraining only | never continuous-affinity confirmation |

One external dataset may not be used for model selection and later described as independent
confirmation. If only one powered external open panel ever exists, every module choice, seed,
hyper-parameter, support draw, query row, baseline, metric and statistical test is frozen before its
first model score.

## 3. OPEN-S0 result — availability and provenance

Probed from this environment on 2026-07-25; exact URLs and HTTP statuses are recorded in
`reports/active/open_s_registry.json` so the finding can be verified or refuted independently.

| source | endpoint probed | status | outcome |
|---|---|---:|---|
| ChEMBL API | `/chembl/api/data/status.json` | 200 | open; ChEMBL_37, release 2026-05-01 |
| Klaeger supplement (publisher) | `science.org/doi/suppl/10.1126/science.aan4368/...` | **403** | not downloadable without publisher access |
| Klaeger supplement (Europe PMC) | `europepmc/.../PMC6542668/supplementaryFiles` | 200 | body is an error: *"Article with id PMC6542668 is not open access one"* |
| Klaeger deposition | PRIDE `projects/PXD005336` | 200 | *"Target Landscape of Clinical Kinase Inhibitors"*, submitted 2017-08-28 |
| Klaeger deposition files | PRIDE `projects/PXD005336/files` | 200 | inventory is mass-spectrometry **RAW** files only |
| ProteomicsDB dose-response | `logic/api/dose_response.xsodata/$metadata` | **404** | no public service exposing drug dose-response found |
| BindingDB native articles | `bindingdb.org/rwd/bind/downloads/BindingDB_BindingDB_Articles_202607_tsv.zip` | 206 | open, 18,114,757 bytes |
| PKIS2 archives | Zenodo search | 200 | open (KCGS v2.0, CC BY 4.0) but percent-inhibition |

**Klaeger `Kdapp` is therefore not admissible here.** The quantitative matrix sits behind a
publisher paywall; the public deposition contains raw spectra, not `Kdapp`; Europe PMC explicitly
declares the article not open access; and no public ProteomicsDB API service exposing drug
dose-response was found at any probed path. Reconstructing `Kdapp` from PXD005336 would require
re-running the full proteomics pipeline, which is not registry reconstruction from public materials
in any usable sense. This is a statement about what is retrievable from the probed public endpoints,
not a claim that the data does not exist.

Downloaded and hashed: `BindingDB_BindingDB_Articles_202607_tsv.zip`, SHA-256
`d2584d1519318d00ab5f46289da5ab3549affe732d598a5072f8777b6b3b5262`, 18,114,757 bytes, retrieved
2026-07-25T13:32:14Z, licence CC BY 3.0 US.

## 4. OPEN-S1 result — endpoint, structure and power shape

The power requirement is measured, not assumed. `reports/active/panel_davis_registration.md`
established that 102 independent components with 12 query ligands per target yields
`MDE80 = 0.1596` against the `0.0614` arm-heterogeneity reference. An admissible external panel
therefore needs **both** ~100 independent homology components **and** ~40 query ligands per target
after the dual-cold firewalls.

BindingDB native subset (93,712 rows; 93,023 curated by BindingDB itself, only 429 imported from
ChEMBL, which is what makes it a candidate independent source), restricted to human single-chain
targets with exact relations:

| endpoint | assay-controlled blocks (target x article) | distinct targets | blocks >= 20 ligands | blocks >= 40 | **targets with a block >= 40** |
|---|---:|---:|---:|---:|---:|
| pKi | 1,175 | 346 | 222 | 53 | **38** |
| pKd | 218 | 163 | 14 | 6 | **6** |

Median block depth is 7 ligands (pKi) and 3 (pKd). Pooling a target's ligands across articles would
raise pKi to 81 targets, but that pools different assay protocols into one within-target ranking and
is inadmissible under this program's own endpoint/assay contract; it is reported only as an upper
bound and is marked inadmissible in the registry.

**38 is far below 100, and every remaining filter can only reduce it** — homology clustering of the
38 accessions, the target-side firewall against ChEMBL train and the Metz panel, and the ligand-side
scaffold/connectivity/Tanimoto firewall which also cuts per-target query depth. The shape of
BindingDB-native curation is structurally wrong for this purpose: it is ligand-deep and
target-shallow, one medicinal-chemistry article at a time, whereas an external dual-cold panel needs
target-breadth.

## 5. Verdict

```text
NO_OPEN_POWERED_INDEPENDENT_PANEL
```

Consequences, all mandatory:

* the power threshold is **not** relaxed, and no private, paywalled or author-only source is
  substituted;
* `OPEN-P`, the one-seed external predictive gate, the three-seed stage and final confirmation are
  **all blocked**; running them would require an evaluation source that does not exist;
* the sealed `panel_davis` stays sealed with `consumed=false`. It is underpowered, not spent;
* `PD-M2` and `PD-H` remain admissible because they are train-only mathematical development on the
  open ChEMBL substrate, but **no result from either may be described as predictive evidence**, and
  neither may authorise a later stage;
* any future attempt requires a genuinely open, endpoint-consistent panel with ~100 components and
  ~40 query ligands per target. The most likely such object is a dense kinome or family-wide
  profiling panel released under an open licence with recoverable continuous values — the Klaeger
  matrix would qualify on shape and fails only on access.
