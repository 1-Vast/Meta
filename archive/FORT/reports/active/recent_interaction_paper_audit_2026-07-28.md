# Recent protein-ligand interaction paper audit

**Date:** 2026-07-28  
**Formal verdict:** `NO_NEW_IDENTIFIABLE_INTERACTION_SOURCE__NO_TRAIN`  
**Final category:** 3 - current data cannot identify the proposed mechanism; new prospective
measurement conditions are required.

## Scope and evidence boundary

This audit asked whether eleven recent papers add a new source of information, stronger supervision,
different estimand, data condition, or independently validated target-conditioned prior for strict
target-ligand dual-cold affinity prediction. A new architecture, larger pretraining corpus, pairwise
loss, LoRA, flow matching, or Bayesian regularizer was not sufficient by itself.

The following claims were kept separate:

1. a paper succeeds on its own task;
2. an abstract mechanism from the paper is relevant to DTA;
3. the mechanism supplies target-specific ligand reordering on the present data;
4. it improves affinity prediction under the complete dual-cold and provenance firewalls.

No candidate reached claim 3. No affinity model or student was trained. The low-cost local audit read
only TRAIN topology columns `conn`, `scaffold`, `docs`, `target`, and `dual_cold_split`; it did not
read `affinity`. No development, Davis, confirmation, or sealed label was read in this round.

## Frozen literature and code

Primary files were downloaded to
`tmp/pdfs/paper_interaction_audit_2026_07_28`; extracted text, rendered title pages, rendered key
method/result pages, supporting information, and repository snapshots are retained there.

| paper | stable identifier | primary PDF SHA-256 |
| --- | --- | --- |
| AAMFM | arXiv `2607.20057v1` | `b94df85679ee7ccd52e5b0468931f3f06bc52caf0d50692b83445f046f4e000c` |
| tau: Future Visual Supervision | arXiv `2607.24485v1` | `b9dd13ee0fc69ce091c6f10b738caaca05aef95608b317e06b2c8021d4b2be4b` |
| EnsembleEGNN | arXiv `2607.21561v1` | `25eb6b2fccf7ee2c5dad301dbc379e6cfb987f0c55c3bce35cb7d979192b0ef9` |
| Gaussian Molecular Shape Overlap | arXiv `2607.20766v1` | `b6ec53ebeb369d32134b8ec10ef5caec678a5ba3e9682c4f14e25aea54d7940b` |
| PAC-DP | arXiv `2607.24296v1` | `e133d6ae374b12ac64b378e65f287479bb60397ef6a3e5e8dd3b362d0fe96dd9` |
| Hybrid UQ | DOI `10.1021/acs.jcim.5c01597` | `1299e63ba3979eafa15f55eba736ba1f2c2df92bfbc98b1e61f37f25f1cab526` |
| NISE | DOI `10.1101/2025.04.22.649862` | `b1194f6a974067aabd38fc572073ff3df52e061a632a5617bcc5d422ede0caf9` |
| Protein-ligand Neosurface Design | DOI `10.1038/s41586-024-08435-4` | `22a680fa8f814091fdcb9ab5f55050018bfbaa1a87372350f652a56f77631d42` |
| PBCNet2.0 | DOI `10.1101/2025.06.04.657800` | `b3e5176e0528a970f7694a564c9cfb4c4c09accd5d03a516daa77ae85caf2abd` |
| FLOWR.ROOT | arXiv `2510.02578v6` | `0eb5ba75b6dbac3ae7eb002f14c302b8828ea4cffdf263f9ed9ea465d4188d9f` |
| Vision Transformer for affinity | DOI `10.1038/s41598-026-62430-5` | `f9b54769965550a2e2adb4e9b729d4e4d03c655f1c2e2b134b64221634d562b1` |

Supporting archives were also frozen: Hybrid UQ
`36e86afb6991dbf60e637bb815b12341d100b073083239ecd8fb40d2c366869e` and Neosurface
`b90efa9d9372eecb15eddf68831bc1a7a52369e35f02d2af8f34187ee3150169`.
The title-page and key-page contact sheets have SHA-256
`bf9aad18b8b8e603ec73901e1760763582e7403f162462e2a5f290f538ceb5a9` and
`f5f826c8687d214e7615073b32fe7fe815621019b3f9b9a4c376376527881ce6`.

| repository | frozen commit |
| --- | --- |
| `YuJie-0202/PBCNet2.0` | `3d46e6e594531c5692376e242b606641979e8550` |
| `jule-c/flowr_root` | `b2263e2516ad798d0119f8c4b531698860cc846e` |
| `AaronFeller/EnsembleEGNN` | `660e09a03deaef3e33a9c8fcf4284dc6f1d3cb2f` |
| `polizzilab/NISE` | `236c5cabe3dc8b3aa71663ff4f00ef8ac6d6d17e` |
| `XL-S224/AAMFM` | `36036588de55f9bd1675dbbf92d726e3ce8db941` |
| `JPoziemski/VIT_for_affinity` | `a60dbfa174f85e43630be527e58f59b8ca627a55` |
| `LPDI-EPFL/masif-neosurf` | `b522fdc537e0724e7c574753f9068f88e687c736` |

The Gaussian-shape paper names an open-source `shape3d` library but gives no resolvable repository URL
in the paper; an author/title GitHub search on 2026-07-28 returned no repository. No public code
repository was identified for tau or PAC-DP. Hybrid UQ code was read from the frozen supporting
archive and its public `CDDLeiden/uqdd` description.

## Paper-by-paper decision

| paper | what the paper actually establishes | transferable abstraction | boundary for strict dual-cold DTA | disposition |
| --- | --- | --- | --- | --- |
| AAMFM | Antibody-antigen CDR generation conditioned on antigen geometry; about 30k preferences are labeled by Protenix/AF3 score plus AntiBERTy PLL and evaluated mainly by the same oracle family, pTM/ipTM, PLL, RMSD, and AAR. | Training-only privileged pair conditioning and preference alignment. | Antibody design is not small-molecule affinity. The preference oracle is circular with the main functional evaluation and supplies no independently measured ligand reordering. | Read-only example of privileged supervision; no teacher or model route. |
| tau | A robot VLA learns tactile state from action-conditioned future visual latent targets available only during training. | A privileged future-state teacher can be removed at deployment. | It contains no protein, ligand, affinity, biological prior, or DTA split. It supplies an analogy, not a teacher state. | Conceptual citation only; grouped with AAMFM/Neosurface, not a separate innovation. |
| EnsembleEGNN | A conformer-set EGNN improves cyclic-peptide PAMPA permeability after CREMP self-supervised pretraining. | Permutation-invariant conformer pooling is a plausible ligand-only representation. | No protein input or affinity estimand. Released evaluation is row-level random 5-fold; the default path does not isolate sequence/scaffold/source, and the downstream CREMP molecules are exposed during representation pretraining. | Hypothesis for B0 only; no evidence of DTA interaction or compliant generalization. |
| Gaussian Shape Overlap | An importance-sampling estimator gives Gaussian shape overlap/union volume, analytic standard errors, adaptive sampling, and rigid alignment. | Conformer diversity, alignment, Monte Carlo error, and adaptive compute diagnostics. | Shape-only DUD-E/LIT-PCBA screening is not quantitative affinity or target-specific reordering. | Numerical tool only. |
| PAC-DP | A PAC-Bayes KL regularizer improves low-data robot diffusion-policy learning. | A complexity penalty can regularize an already meaningful prior. | No protein, ligand, affinity, biological prior, or mechanism-identification result. A KL term cannot create missing target information. | Theory analogy only; blocked until a biological prior passes independently. |
| Hybrid UQ | Evidential/deep-ensemble hybrids improve uncertainty metrics on Papyrus++ under random and chemical scaffold shifts. | Aleatoric/epistemic reporting and rejection/calibration metrics. | There is no target- or provenance-disjoint evaluation. In released multitask code, `n_targets > 0` sets the protein descriptor to `None`; the model uses chemical features with target-specific outputs. | UQ tool after a valid predictor exists; not an interaction source or performance module. |
| NISE | Reciprocal sequence/structure networks design exatecan-binding proteins; four designs bind, and two local mutations improve one designed binder to `Kd = 1.2 +/- 0.2 nM`. | Co-structure self-consistency can be useful in local protein design. | The original default objective is ligand pLDDT. Current code can log Boltz `affinity_pred_value`, but its `pbind` objective maximizes `affinity_probability_binary`. These are confidence/compatibility variables, not calibrated multi-target quantitative affinity or shared-ligand reordering. | Strong result on its own design task; no quantitative DTA teacher. |
| Neosurface | A learned surface fingerprint retrieves/designs protein partners for ligand-induced ternary interfaces and experimentally validates three binder systems. | A ligand can alter a protein surface in a biologically meaningful way. | Fourteen ternary complexes yield 28 partner-search directions against 8,879 decoys. The output is protein-partner retrieval/design, not small-molecule affinity. Inference requires a ligand-bound PDB surface; full reproduction uses TensorFlow 1.9, MSMS/APBS/PyMesh and several TB of data. | External structural biology evidence only; wrong estimand and privileged input for DTA. |
| PBCNet2.0 | A frozen Siamese Cartesian-tensor model predicts relative activity for two similar ligands in one protein pocket and reports strong same-series benchmarks/prospective examples. | A pair-conditioned teacher is the closest candidate to a target-specific relative-affinity signal. | The 8.6M count is paired expansion of about 2.81M BindingDB 2023.12 measurements mixing IC50/EC50/Ki/Kd. Poses depend on a similar co-crystal ligand, Glide docking, MCS alignment, and the same-pocket/same-series assumption. Exact training files are restricted, so overlap cannot be excluded. | Strongest teacher candidate, but inadmissible for confirmation and impossible on current inputs. No student/distillation. |
| FLOWR.ROOT | A joint 3D generation/affinity model uses broad pretraining, curated structural refinement, and project-specific LoRA; the paper reports that zero-shot models fail on four unseen proprietary SAR projects and LoRA recovers them. | Dense project-specific measurements can support adaptation. | The paper says the FEP+/OpenFE benchmark likely substantially overlaps training ligand/target space and is not evidence of substantial generalization. Four adaptation projects are proprietary. The repository is an early release without final converged weights, requires protein PDB/CIF plus ligand/pocket inputs, and recommends at least 40 GB VRAM; this host has 8,188 MiB. | Measurement-design evidence, not a runnable strict few-shot route. |
| Vision Transformer | A ViT predicts affinity from a known 3D protein-ligand complex on PDBbind/CASF/CoreSet and gives architecture/augmentation ablations. | Matched known-complex structural baseline. | Remaining PDBbind complexes are randomly divided 90/10 for train/validation. The paper reports no test protein with zero or only one structurally similar training protein; its hardest 18-target subset still permits up to three `TM-score > 0.75` neighbors. | Architecture baseline only; no new information or strict target-cold evidence. |

## External-teacher contamination and input audit

### PBCNet2.0

The paper's 8.6M examples are not 8.6M independent target-ligand measurements. Molecules are paired
within chemical series after endpoint aggregation and pose filtering. Pair count is therefore a
quadratic training construction, not an independent sample size.

Zenodo record `15656365` was queried through its public API on 2026-07-28. It is published with
`access_right=restricted` and exposes zero downloadable file entries. The repository includes weights
and test inputs, but the exact training identities cannot be intersected with FORT targets, homology
groups, ligand parents, scaffolds, structures, documents, or provenance families. BindingDB and ChEMBL
also share source literature, so database-name separation is not a provenance firewall.

The required inference object is not `(sequence, SMILES)`. It is two same-pocket complexes constructed
using a co-crystal reference, a selected docking pose, and MCS-to-reference alignment. The current
registry has no query pose, reference-complex, or reference-ligand field. Generating them would reopen
the previously closed physical-pose route and would still not solve pretraining contamination.

Consequently no pretraining-overlap-excluded subset, component-level effect, or MDE can be constructed.
PBCNet2.0 is exploratory at most. It cannot be a confirmation teacher, and student training is
forbidden.

### FLOWR.ROOT

Stage-1/2 sources include ZINC3D, PubChem3D, Enamine REAL, OMol25, Plinder, BindingMOAD, BindingNet,
SAIR, KIBA-3D, Davis-3D, Kinodata-3D, SPINDR, and HiQBind. This broad mixture directly intersects or
derives from sources already used in FORT. The paper itself states that the FEP+/OpenFE ligand and
target spaces likely substantially overlap training and that those results do not demonstrate
substantial generalization.

The four decisive LoRA cases are private project datasets. Their dense SAR condition cannot be
reconstructed, and no accession/homology/ligand/scaffold/document/provenance intersection can be
audited. The public repository labels itself an early release and promises final converged weights
later. Its inference contract requires protein and ligand structure files and recommends 40 GB VRAM,
versus the local RTX 4060 Laptop GPU's 8,188 MiB. A resource workaround would not fix the data or
identification failures.

### Boltz/NISE, AAMFM, and Neosurface

These systems produce co-folding confidence, binary binder probability, structure-oracle preferences,
or protein-partner surface complementarity. Without an independent pretraining-overlap exclusion and
a measured correct-target versus wrong-target reordering gate, none can be relabeled as quantitative
affinity. AAMFM and NISE also optimize or rank with the same oracle family later used as evidence of
functional quality, which is unsuitable as independent confirmation.

## Frozen candidate definitions

### C1: frozen pair-conditioned teacher

For a legal same-target pair, define

```text
s_tij = teacher(t, d_i, d_j, pose_i, pose_j)
r_tij = (y_ti - y_tj) - ligand_only_difference(i, j)
```

The independent vector must be aggregated first by target, full-sequence homology component,
document/provenance family, and any entity-disjoint packing. Raw pair count is never the inference
unit. The frozen identifying gate requires true teacher residual association to beat ligand-only,
matched 2D, constant/random teacher, target shuffle, family-matched wrong target, ligand shuffle,
pose/interface shuffle, teacher ligand-only projection, and a pretraining-overlap-excluded sensitivity
set, with a positive component interval and adequate MDE.

Current status: not executable. PBC exact overlap and required poses are unavailable; FLOWR weights and
inputs are not compliant; NISE/Boltz/Neosurface outputs have the wrong estimand. No teacher score was
generated and no student was trained.

### C2: ligand conformer ensemble B0

For conformers `c_dm` and frozen weights `w_dm`, the only first-stage estimand is ligand-only:

```text
b_ens(d) = h(Pool_m(phi(c_dm), w_dm))
```

It must separately beat the original 2D B0 and a matched-capacity 2D neural baseline under parent,
scaffold, high-similarity, equivalent-molecule, document, and provenance-family component splits.
Required controls are single conformer, uniform/Boltzmann pooling, shuffled conformer set, duplicated
conformer, and parameter-matched random ensemble. Target-conditioned conformer weighting is not
authorized until ligand-only B0 passes.

Current status: topology stop before labels or training.

### C3: FLOWR-inspired support adaptation

The minimal nested sequence is:

```text
f0(t,d)                              no adaptation
f0(t,d) + a_t                        intercept only
f0(t,d) + beta_t^T x_d               ligand-only adaptation
f0(t,d) + u_t^T V x_d                low-dimensional interaction adaptation
```

Correct-target support must beat wrong-target, label-permuted, cross-target, and matched-capacity LoRA
support on identical episodes and query rows. Historical matched-episode audits already found correct
support indistinguishable from cross-target or label-permuted support; later support-conditioned
posterior, SCGD, and QACO routes also failed wrong-support or label-specificity gates. FLOWR adds no
new support information. Its positive result instead identifies dense project-specific SAR as a data
condition.

Current status: historical null remains binding; no LoRA or other capacity rescue.

## Label-free TRAIN topology audit

The reproducible artifact is
`reports/active/recent_interaction_paper_gate_2026-07-28.json`, SHA-256
`f232d55454aea0893f3df664dbcb48c9ef0e1f0f30633046aecd38d563eac51e`.
The runner and test hashes are:

- `research/recent_paper_gate_audit.py`:
  `bb7068344b56b357fd4bdc61725ba52469d4ae132176522ce427630122cb0b21`;
- `tests/test_recent_paper_gate_audit.py`:
  `f7a46c0746f6a1a7bf80d29639eb395beb734b5f3d568078732628032365fb6b`.

Only TRAIN rows and the five declared non-affinity columns were projected. Ligand parents were joined
when they shared either a Bemis-Murcko scaffold or an atomic ChEMBL document ID.

| topology quantity | result |
| --- | ---: |
| TRAIN rows | 201,827 |
| ligand parents | 121,401 |
| scaffolds | 48,234 |
| atomic document IDs | 9,587 |
| scaffold-or-document components | 2,197 |
| largest component ligands | 90,288 (74.37%) |
| largest component rows | 163,117 (80.82%) |
| largest component targets | 527 |
| components with at least 40 ligands | 171 |
| singleton components | 751 |

This is an optimistic lower bound on component merging: high-similarity chemical edges were not added.
Adding them can only merge components. More importantly, the registry contains document IDs but not
reliable source-lineage/provenance-family metadata. RECRO previously showed that document-disjointness
does not imply provenance-family independence. Therefore a compliant conformer-B0 confirmation split,
its independent statistical units, and its MDE are not identifiable on this registry.

Formal topology verdict: `NO_COMPLIANT_TRAINING_SUBSTRATE`.

## Historical conflict audit

| proposed idea | binding prior evidence | decision |
| --- | --- | --- |
| AAMFM/tau/Neosurface privileged supervision | Physical pose, pocket, native-complex and structural distillation routes already failed to establish leakage-free target reordering. | New examples of the same supervision family; no route reopening. |
| NISE/Boltz compatibility | Pair-compatibility and confidence variables have not predicted quantitative target-specific ordering under destruction controls. | Keep as compatibility diagnostics only. |
| PBC pairwise relative objective | Ordinary MMP/pairwise ranking cannot create target information; PBC would be different only through its external pose-conditioned teacher. | The potentially new information is inadmissible because input and contamination gates fail. |
| FLOWR LoRA | Matched support episodes, posterior adaptation, SCGD, and QACO failed correct-support specificity; increasing adapter capacity is prohibited. | Dense private SAR is a measurement condition, not evidence for strict few-shot adaptation. |
| ViT | Replacing the interaction head or backbone without new information is closed. | Baseline only. |
| PAC-DP | Bayesian regularization cannot rescue an unidentified or non-regular biological prior. | Blocked until a prior passes independently. |
| EnsembleEGNN | Ligand-only B0 remains valid, but conformer improvement has never passed a provenance-family-disjoint split. | A distinct B0 hypothesis, stopped by the current topology gate rather than declared biologically false. |

The PBC, EnsembleEGNN, and FLOWR papers therefore do not overturn earlier negative mechanisms. They
clarify three future data requirements: auditable frozen-teacher lineage plus legal poses, a
provenance-independent ligand graph for B0 representation tests, and dense correct-target SAR for
adaptation.

## Power and stopping decision

For C1, the number of legal target/homology/provenance units is unknown because exact teacher training
membership is inaccessible. For C2, source families are absent and the scaffold/document graph has a
74.37% ligand giant component. For C3, the paper's only convincing adaptation units are proprietary
projects and existing public matched episodes failed support specificity. In all three cases, an
empirical component-level MDE for the requested claim cannot be computed before violating a firewall.

The stop is therefore not a null estimate of biological effect. It is a preregistered observability and
power failure. Tuning a kernel, rank, LoRA width, teacher combination, model capacity, seed count, or
training duration cannot repair it.

Zero candidates entered affinity training. Frozen-teacher distillation, target-conditioned conformer
weighting, FLOWR-style LoRA, and any combined paper module are not authorized.

## Reopening conditions

The active next step remains a prospective factorial reliability panel:

- 12 targets, exactly two from each of at least six protein families;
- at least 16 shared, scaffold-diverse ligands against every target;
- two operationally independent sites or provenance lineages;
- one preregistered pKi or pKd endpoint;
- complete randomized `12 x 16 x 2 = 384` inclusion before technical replication;
- inactive, censored, failed, and out-of-quantification outcomes retained;
- frozen target, homology, binding-profile, ligand, scaffold, chemical-neighbor, assay, document, and
  provenance firewalls.

This A0 panel estimates cross-site reliability, target-specific mixed-difference variance, and empirical
MDE. It does not by itself validate a `0.03` predictive gain. A predictive route still requires the
broader powered component count already frozen in `task.md`.

PBC-like teachers additionally require a downloadable exact pretraining manifest, an
overlap-excluded subset, legal label-blind poses/reference complexes, and the complete C1 destruction
suite. Conformer B0 additionally requires provenance-family metadata and a component split that
survives high-similarity closure. Adaptation additionally requires dense correct-target support that
beats wrong and permuted support before any LoRA comparison.

## Verification

The topology runner compiled, reran deterministically, and its focused test returned `3 passed` in the
`drug` environment. The project-owned suite, explicitly scoped as `pytest -q tests`, returned
`321 passed`. The artifact parses as strict JSON with `affinity_column_read=false`,
`development_labels_read=false`, `confirmation_labels_read=false`, and
`sealed_test_consumed=false`.

The rendered title-page contact sheet was visually inspected and contains all eleven papers. The
rendered key method/result pages have legible text, tables, captions, and page boundaries with no
clipping, black blocks, or missing panels. This visual pass was used alongside, not replaced by, text
extraction and repository inspection.

A bare repository-root `pytest` is not a valid FORT test boundary after freezing third-party
repositories under `tmp/`: it traversed AAMFM, FLOWR.ROOT, and MaSIF-neosurf test files and stopped at
six collection errors for their uninstalled BioPython, FLOWR package, TensorFlow, IPython, and local
module environments. No FORT test failed, and no external dependency was installed to disguise that
boundary.
