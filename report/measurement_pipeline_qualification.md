# Qualifying the Measurement Pipeline Before Testing for a Transferable Protein-Conditioned Interaction Signal

**Independent research assessment. Design phase only — no code, no repository changes, no training.**
17 August 2026

Verified facts are marked ✔ (checked this session against a primary or official source). Items marked ⚠ come from secondary sources or require reading a supplementary file, and are flagged for verification. Hypotheses are labelled as such.

---

## Contents

1. [Executive scientific judgment](#1-executive-scientific-judgment)
2. [Diagnosis of the invalid capability and planted-signal designs](#2-diagnosis-of-the-invalid-designs)
3. [Literature map](#3-literature-map)
4. [Verified mutation / variant dataset comparison](#4-verified-mutation--variant-dataset-comparison)
5. [Recommended representation-capability metric](#5-recommended-representation-capability-metric)
6. [Corrected planted-signal benchmark design](#6-corrected-planted-signal-benchmark-design)
7. [Ranked biological positive-control panels](#7-ranked-biological-positive-control-panels)
8. [Two-level capability / localization protocol](#8-two-level-capability--localization-protocol)
9. [Censoring and endpoint analysis plan](#9-censoring-and-endpoint-analysis-plan)
10. [Full control matrix](#10-full-control-matrix)
11. [Integrity-test specification](#11-integrity-test-specification)
12. [Power and cluster-bootstrap plan](#12-power-and-cluster-bootstrap-plan)
13. [Ranked post-gate model directions](#13-ranked-post-gate-model-directions)
14. [Staged decision tree](#14-staged-decision-tree)
15. [Outcome definitions](#15-outcome-definitions)
16. [Data access, time and compute](#16-data-access-time-and-compute)
17. [The single best next action](#17-the-single-best-next-action)

---

## 1. Executive scientific judgment

### 1.1 The cycle is in a healthy state, and that should be said plainly

A programme that invalidates its own instruments before reporting a result is doing science correctly. Sixteen specific defects were found across two designs, and **not one biological claim was made on top of them**. That is the difference between a slow programme and a wrong one. The recommendations below assume this standard continues.

### 1.2 Four judgments that change what should be built next

**(1) The audit's own conclusions are correct, and one of them is independently confirmable.** The BRAF annotation mismatch is not a data-entry error; it is a documented historical renumbering. The 2002 Sanger report of the common melanoma BRAF mutation as V599E was later renumbered to V600E because of a reference-sequence difference — the classic one-residue shift from counting the mature protein after loss of the initiator methionine, the same phenomenon that makes sickle-cell HbS E6V into p.E7V under HGVS. Any variant annotation drawn from pre-~2004 literature or from vendor construct catalogues carries this risk. This is a **systematic** class of defect, not four bad rows.

**(2) The "local representation looks sensitive" result was not merely unproven — it was measuring something with no defined meaning.** Three of the seven defects (mutant windows centred on the mutation while WT windows were centred at the sequence midpoint; the zero-valued inter-protein denominator; the mutation token being an edit descriptor) are not independent bugs. They are one conceptual error: **a sensitivity ratio was computed between objects that do not live in a common coordinate system, normalised by a quantity that is zero by construction.** The correct repair is not to fix the centring and rerun. It is to replace the metric family entirely (§5), because *distance ratios are the wrong instrument for this question*.

**(3) The planted-signal design had one structural flaw that subsumes most of the nine listed.** Adding synthetic interaction to **real biological endpoint values** makes the estimand unidentifiable before any code runs. If `y_obs = y_real + τ·I_synth`, and `y_real` already contains an unknown real interaction `I_real`, then no evaluation can attribute recovered structure to `I_synth`. The planted-signal benchmark must have **fully synthetic labels**; realism must come from reusing the *observation graph* — which cells exist, which are censored, which proteins and scaffolds — never from reusing the *values*. Get that right and the "raw endpoint vs interaction truth," "main effects not removed," and "endpoint offsets dominate sign accuracy" defects all dissolve, because there are no endpoint offsets left to dominate.

**(4) A dataset published four months ago supersedes the planned positive control.** Saifudeen et al., *Nature Biotechnology*, 20 April 2026: <cite index="285-1">86 of the ~100 approved kinase inhibitors were profiled against 758 kinases, including 409 wild-type and 349 oncogenic variants using a biochemical kinase assay</cite>, and <cite index="285-1">all compounds in the full panel were screened in duplicate at 1 μM with Km ATP, generating ~290,000 kinase–drug interaction measurements across 758 kinases (409 wild-type and 349 mutants)</cite> ✔ ([doi:10.1038/s41587-026-03090-8](https://www.nature.com/articles/s41587-026-03090-8)). This is roughly **4.6× the variant count** of the Duong-Ly panel, with wild types and variants measured **inside the same study**.

### 1.3 The Duong-Ly problem the audit has not yet surfaced

The Duong-Ly panel does not contain its own wild-type measurements. The paper states that the heatmap and dataset <cite index="277-1">include data reported previously for the corresponding wild-type forms for comparison (Anastassiadis et al., 2011)</cite> ✔. The 76 mutants derive from <cite index="277-1">21 cognate wild-type kinases</cite> and the screen produced <cite index="277-1">over 13,000 mutant kinase-compound pairs</cite>, with <cite index="277-1">each measurement performed in duplicate</cite> ✔.

Three consequences the programme must plan around:

- **Every WT→mutant difference computed from Duong-Ly + Anastassiadis is a cross-study contrast spanning roughly five years**, not a within-experiment contrast. A batch term rides on every Δ and does not cancel.
- **The effective independent-parent count is 21, not 76.** Leave-one-parent-out cross-validation therefore has at most 21 folds, and fewer once compound coverage is required. Variants are unevenly distributed across parents (EGFR alone carries eleven ⚠), so a result can be one kinase's result.
- **The compound sets are described as "overlapping," not identical** (183 vs 178). Overlap must be established compound-by-compound on standardised structures, never by name.

### 1.4 An assay-physics point that cuts in the programme's favour

My earlier framing treated fixed-ATP activity panels as a general obstacle. That was too coarse. For a purely ATP-competitive inhibitor, Cheng–Prusoff gives `IC50 = Ki·(1 + [ATP]/Km)`. Two regimes behave completely differently:

- **Fixed [ATP] across kinases** → the factor `(1 + [ATP]/Km_i)` varies per kinase. It contributes a protein main effect (which cancels in any within-protein-pair centred contrast) **but also interacts with inhibitor mechanism**, since type-II and allosteric inhibitors do not obey the competitive form. That residual does *not* cancel.
- **[ATP] set at each kinase's own Km** → the factor is `2` for every ATP-competitive inhibitor, **independent of the kinase**. Cross-kinase log-differences are then Ki-comparable up to a constant.

Saifudeen et al. ran at Km ATP ✔. Reaction Biology's platform is also run at per-kinase Km ATP in some protocols (e.g. a published protocol specifies Km ATP of 15, 20, 20 and 10 µM for ARAF, BRAF, BRAF-V600E and RAF1 respectively ⚠), while other HotSpot runs use a flat 10 µM or 100 µM ⚠. **ATP concentration must therefore be read per kinase from each study's supplementary methods table and treated as a first-class covariate — never assumed from the vendor name.** The brief's instruction not to infer comparability from a shared commercial provider is exactly right, and this is the concrete reason.

**Hypothesis, flagged as such:** for a WT/mutant pair where the mutation changes Km(ATP) — common for activating mutations — a Km-ATP protocol assays the two constructs at *different* ATP concentrations. For ATP-competitive inhibitors this is self-correcting; for type-II and allosteric inhibitors it is not. WT→mutant Δ values should therefore be stratified by inhibitor binding mode before interpretation.

### 1.5 Calibrated priors, recorded so they can be wrong

| Claim | P |
|---|---|
| Local, coordinate-comparable representations pass a properly constructed capability test | 0.85 |
| Global mean-pooled PLM embeddings fail the same test | 0.85 |
| A corrected planted-signal harness recovers τ* ≥ 0.5 (in noise-SD units) at protein-component level | 0.80 |
| Biological positive control (WT→mutant response prediction, unseen parent) passes on Saifudeen | 0.55 |
| Same, on Duong-Ly + Anastassiadis alone | 0.35 (batch confound + 21 parents) |
| Localization to the correct residue succeeds where capability succeeds | 0.30 |
| Cold-parent transfer of the interaction signal at useful effect size | 0.45 |

---

## 2. Diagnosis of the invalid designs

### 2.1 The sixteen defects, grouped by archetype

Grouping matters more than the list, because archetypes recur and individual fixes do not prevent recurrence.

| Archetype | Defects it explains | Why it happens | Structural fix |
|---|---|---|---|
| **A. Coordinate incomparability** | mutant windows centred on mutation, WT on midpoint; mutation token as edit descriptor; construct/isoform offsets | the position index was derived from the *name* rather than being a property of the *pair* | a typed variant record carrying `(parent_seq, variant_seq, index_in_each, provenance)`; both members of a pair are rendered by one function |
| **B. Degenerate normaliser** | zero inter-protein denominator for mutation tokens | a scale was chosen without checking it is non-zero on the actual objects | forbid ratio metrics; require a denominator that is a spread over *distinct, non-identical* comparators (§5) |
| **C. Truth/prediction mismatch** | synthetic interaction added to real endpoints; raw predictions compared to interaction truth; main effects not removed from fitted predictions; endpoint offsets dominating sign accuracy | the object being evaluated was not the object generated | fully synthetic labels; a single projection operator applied identically to truth and prediction (§6.5) |
| **D. No held-out** | training and evaluation on the same rows | evaluation was written as a sanity print, then promoted | split declared before any model is instantiated; assertion on index-set disjointness |
| **E. Unmatched control** | ligand-only control was not a trained matched model | controls were computed analytically while the treatment was trained | every control is the *same* training procedure with one input ablated; identical optimiser, schedule, capacity, seeds |
| **F. Nondeterminism** | process-dependent Python hashing in the random control | `hash()` on `str` is salted per process unless `PYTHONHASHSEED` is fixed | seeds derived from content hashes (`blake2b` of a canonical key), never `hash()`; determinism test across two processes |
| **G. Unversioned biological annotation** | 4/76 mutation/residue mismatches; historical BRAF numbering; PDGFRα entries | annotations were treated as free text | versioned resolution against a pinned UniProt release, with an explicit historical-alias table and a hard failure on mismatch |
| **H. Untested code path** | tensor equation did not run; generated main effects unused; KLIFS features not implemented | components were written but never exercised end-to-end | a "generator round-trip" test: recover planted parameters from generated data before any model is involved |

### 2.2 Why the capability result cannot simply be re-run after fixing the centring

Fixing the window centring produces a *comparable* pair of vectors, but the reported quantity would still be a **distance ratio**, and distance ratios answer the wrong question. Three reasons:

1. **Scale is not capability.** A representation can move a lot under a perturbation while carrying no recoverable information about *which* perturbation occurred — e.g. if the movement direction is dominated by a nuisance axis. Conversely a small but consistent displacement in a low-variance direction is perfectly readable by a linear probe.
2. **Anisotropy makes Euclidean distance uninterpretable.** PLM embedding spaces are strongly anisotropic; raw ‖·‖ mixes high-variance nuisance directions with the informative subspace. Whitening helps but introduces its own estimation problem.
3. **The edit-descriptor pathology has no distance-based defence.** A one-hot "V600E" descriptor has enormous WT–mutant distance and *zero* protein information. Any metric that ranks it highly is measuring the label, not the representation.

The metric family must therefore be **probe-based with a control task**, which disqualifies the edit descriptor automatically (§5.3).

### 2.3 The planted-signal design's identifiability failure, stated formally

Let the real panel be `y_real(p,l) = μ + a(p) + b(l) + I_real(p,l) + ε`. The draft added a synthetic term:

```
y_obs = y_real + τ·I_synth
      = μ + a(p) + b(l) + [ I_real(p,l) + τ·I_synth(p,l) ] + ε
```

The bracketed sum is the only interaction-like quantity in the data, and `I_real` is unknown. Therefore:

- Recovery of `I_synth` cannot be separated from recovery of `I_real` — the very quantity whose existence is in question.
- A τ = 0 negative control is not a null: it is the real-data experiment, which may legitimately show signal.
- Worse, a *successful* recovery at small τ could be entirely `I_real`, and would be read as pipeline validation when it is actually the biological result the programme has not yet earned.

**This is the most consequential single defect in the list**, because it would have produced a confidently wrong "pipeline qualified" verdict. Fully synthetic labels are non-negotiable.

---

## 3. Literature map

### 3.1 Probing methodology — the correct formalism for RQ1

| Work | Contribution used here | Link |
|---|---|---|
| Hewitt & Liang, EMNLP 2019, "Designing and Interpreting Probes with Control Tasks" | **Selectivity** = task accuracy − control-task accuracy, where the control task assigns each *type* a random but consistent label. Isolates representation content from probe memorisation. Also shows complex probes memorise, and that **dropout does not improve selectivity** while other regularisation does; low-dimensional MLP probes are far more selective ✔ | [ACL Anthology D19-1275](https://aclanthology.org/D19-1275/) · [author's summary](https://www.cs.columbia.edu/~johnhew/interpreting-probes.html) |
| Hewitt & Manning, NAACL 2019, structural probe | distance-based probes done properly (learned metric, not raw ‖·‖) | doi:10.18653/v1/N19-1419 |
| Elazar et al., "Amnesic probing" (TACL 2021) | information *removal* as a test of causal use, complementing read-out probing | doi:10.1162/tacl_a_00359 |
| Ravfogel et al., INLP | iterative nullspace projection for removing a property from a representation | doi:10.18653/v1/2020.acl-main.647 |
| Belinkov, *Computational Linguistics* 2022, "Probing classifiers: promises, shortcomings, advances" | survey of probe pitfalls; why probe accuracy ≠ model use | doi:10.1162/coli_a_00422 |

### 3.2 Protein language models, variants and benchmarks

| Work | Contribution | Link |
|---|---|---|
| **ProteinGym** (Notin et al., NeurIPS 2023 D&B) | substitution benchmark of ~2.7M missense variants across **217 DMS assays** plus 2,525 clinical proteins; indel benchmark ~300k mutants across 74 DMS assays ✔. The external validation set for any mutation representation | [GitHub](https://github.com/OATML-Markslab/ProteinGym) · [NeurIPS](https://papers.nips.cc/paper_files/paper/2023/hash/cac723e5ff29f65e3fcbb0739ae91bee-Abstract-Datasets_and_Benchmarks.html) |
| Meier et al., NeurIPS 2021, ESM-1v | zero-shot variant effect from masked-LM likelihood — the standard baseline for "does a PLM know about this mutation" | arXiv:2107.05340 |
| Rao et al., ICLR 2021 | PLM attention maps as unsupervised contact predictors — motivates contact-aware local neighbourhoods | [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.12.15.422761v1) · [ESM](https://github.com/facebookresearch/esm) |
| Lin et al., *Science* 2023 | ESM-2 / ESMFold | doi:10.1126/science.ade2574 |
| **MaveDB** / Atlas of Variant Effects | open repository of DMS assays; the upstream source for much of ProteinGym | [mavedb.org](https://www.mavedb.org) |
| Marquet et al., *Human Genetics* 2022 | embeddings from PLMs predict conservation and variant effects — evidence local embeddings carry variant information | doi:10.1007/s00439-021-02411-y |
| Fine-tuning PLMs with DMS (2024) | notes that **DMS score scales are highly assay-dependent and require rescaling before pooling across assays** — directly relevant to endpoint governance | [arXiv:2405.06729](https://arxiv.org/pdf/2405.06729) |

### 3.3 Kinase resistance, mutation and pocket structure

| Work | Contribution | Link |
|---|---|---|
| Persky et al., *Nat. Struct. Mol. Biol.* 2020 | DMS across multiple kinases identifies **generalisable** ATP-competitive-inhibitor resistance residues; predictions transferred to TBK1, CSNK2A1, BRAF; an activation site confirmed in BRAF, EGFR, HER2, MEK1 | [doi:10.1038/s41594-019-0358-z](https://www.nature.com/articles/s41594-019-0358-z) |
| Azam, Seeliger, Gray, Kuriyan & Daley, *NSMB* 2008 | activation of tyrosine kinases by gatekeeper-threonine mutation — the mechanistic basis of the gatekeeper positive control | doi:10.1038/nsmb.1486 |
| KLIFS (Kooistra et al., *NAR* 2016) | 85-position kinase pocket alignment; reported 0.8 ± 0.1 Å superposition RMSD for aligned residues ✔ | [doi:10.1093/nar/gkv1082](https://academic.oup.com/nar/article/44/D1/D365/2502606) · [klifs.net](https://klifs.net/) |
| van Westen et al., HIV antivirogram PCM | prior positive result on WT→mutant response prediction with residue-aware descriptors, generalising to unseen mutants | [PMC3578754](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3578754/) |

### 3.4 Planted-signal, identifiability and experimental design

| Topic | Why it matters | Anchor |
|---|---|---|
| Spiked / planted low-rank models | detectability of a rank-r spike depends on amplitude **and** rank **and** aspect ratio, not amplitude alone — hence the τ grid must be two-dimensional (§6.4) | Baik–Ben Arous–Péché phase transition literature |
| Simulation-based calibration (Talts et al. 2018) | the discipline of validating an inference pipeline by recovering known generated parameters — the correct template for a planted-signal harness | arXiv:1804.06788 |
| Persistent excitation (system identification) | an interaction operator is identifiable only if the observation graph excites it; must be audited before training | Ljung, *System Identification* |
| Multimodal collapse: Wang, Tran & Feiszli (CVPR 2020); Peng et al. OGM-GE (CVPR 2022) | a null with a collapsed protein branch is an optimisation artifact, not a falsification | — |
| Hewitt–Liang control tasks (above) | the general principle that every capability claim needs a matched control that *cannot* succeed for the right reason | — |

### 3.5 Censoring

- **SparseChem** (MELLODDY): censored regression with explicit ±1/0 masks and one-sided squared loss — [arXiv:2203.04676](https://arxiv.org/pdf/2203.04676)
- **AMPL** (ATOM): maximum-likelihood mean estimation for partially censored pIC50 — [arXiv:2002.12541](https://arxiv.org/pdf/2002.12541) · [code](https://github.com/ATOMconsortium/AMPL)
- Landrum & Riniker, *JCIM* 2024: cross-assay noise floor — minimal curation gives ~65% of same-target IC50 pairs differing >0.3 log and 27% >1 log; maximal curation 48%/13%, Kendall τ ≈ 0.71 ✔ — [doi:10.1021/acs.jcim.4c00049](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049)

---

## 4. Verified mutation / variant dataset comparison

### 4.1 The panels

| Field | **Saifudeen 2026** | **Duong-Ly 2016** | **Anastassiadis 2011** | **Davis 2011 (KINOMEscan)** | **PKIS2 (Drewry 2017)** |
|---|---|---|---|---|---|
| Citation | Saifudeen, Zhu, Liang et al., *Nat. Biotechnol.* (2026) ✔ | Duong-Ly et al., *Cell Rep.* 14:772 ✔ | Anastassiadis et al., *Nat. Biotechnol.* 29:1039 ✔ | Davis et al., *Nat. Biotechnol.* 29:1046 ✔ | Drewry et al., *PLOS ONE* 12:e0181585 ✔ |
| DOI / URL | [10.1038/s41587-026-03090-8](https://www.nature.com/articles/s41587-026-03090-8) | [10.1016/j.celrep.2015.12.080](https://www.cell.com/cell-reports/fulltext/S2211-1247(15)01536-3) | [nbt.2017](https://www.nature.com/articles/nbt.2017) | [nbt.1990](https://www.nature.com/articles/nbt.1990) | [PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0181585) |
| Data location | Supplementary Table 4; portal at `kirhub.fredhutch.org` ⚠ | Table S2; `kir.fccc.edu` ✔ | supplementary tables | supplementary XLS ✔ | supplementary, public-domain release ✔ |
| Licence | **verify** ⚠ | **CC BY-NC-ND 4.0** ✔ (no derivative redistribution) | verify ⚠ | verify ⚠ | **CC BY 4.0** ✔ |
| Parent proteins | 409 wild-type ✔ | 21 cognate WT ✔ | 300 WT ✔ | ~363 distinct domains ⚠ | ~250–406 assays ⚠ |
| Variants | **349 oncogenic variants** ✔ | 76 mutants ✔ | minimal | ~21–27 disease-relevant mutants + phospho-state variants ⚠ | ~21 disease-relevant mutants ⚠ |
| Compounds | 86 approved inhibitors ✔ | 183 ✔ | 178 ✔ | 72 ✔ | 645 ✔ |
| Measurements | **~290,000** ✔ | >13,000 mutant–compound pairs ✔ | ~53,400 ⚠ | 31,824 cells ✔ | ~250k ⚠ |
| Replication | **duplicate** ✔ | **duplicate**, 33 discrepant pairs (0.2%) removed ✔ | duplicate ⚠ | single ⚠ | ⚠ |
| Assay | biochemical kinase assay ✔ | Reaction Biology HotSpot ✔ | Reaction Biology HotSpot ✔ | KINOMEscan competition **binding** ✔ | KINOMEscan ✔ |
| Compound conc. | **1 µM** ✔ | 0.5 µM ⚠ | 0.5 µM ⚠ | primary screen at 10 µM; Kd follow-up ✔ | single-conc. ⚠ |
| ATP | **Km ATP (per kinase)** ✔ | per-kinase, read from supplementary ⚠ | per-kinase, read from supplementary ⚠ | n/a (binding assay) ✔ | n/a ✔ |
| Endpoint | % inhibition / remaining activity ⚠ | % remaining activity ✔ | % remaining activity ✔ | **Kd (nM)** ✔ | % control ✔ |
| Orientation | lower % remaining = stronger | lower = stronger | lower = stronger | **lower Kd = stronger** | lower %ctrl = stronger |
| Censoring | floor 0% / ceiling 100% ⚠ | floor/ceiling ⚠ | floor/ceiling ⚠ | **blanks = Kd > 10 µM or undetected at 10 µM** ✔ | floor/ceiling ⚠ |
| **WT and variant in the same experiment?** | **YES** ✔ | **NO** — WT imported from Anastassiadis 2011 ✔ | n/a | partially (variants in-panel) ⚠ | partially ⚠ |
| Leave-one-parent-out feasible? | yes, large | yes, **≤21 folds** ✔ | n/a | limited | limited |

### 4.2 Fields that must be read from supplementary files before any power calculation

For every panel: median shared ligands per WT–variant pair; exact censored fraction; per-kinase ATP concentration and substrate; construct boundaries (residue ranges) and tags; whether the "variant" is a point mutation, deletion, fusion, ITD or phosphorylation state; compound identity as standardised structures (InChIKey), not names; UniProt accession **and isoform**; KLIFS pocket mapping availability.

**None of these should be inferred.** In particular, the Saifudeen panel contains fusions as well as point variants ⚠ — fusions are a different perturbation class and must be typed separately, not pooled with substitutions.

### 4.3 A trap specific to KINOMEscan panels

KINOMEscan lists phosphorylated and non-phosphorylated forms of the same kinase as separate assays, and includes autoinhibited and truncated constructs. These have **identical amino-acid sequences** but genuinely different affinities. A sequence-only representation assigns them identical features, which (a) puts an irreducible ceiling on achievable performance and (b) violates the functional form `y = f(sequence, ligand)` that the whole programme assumes. They must be excluded or given an explicit activation-state covariate.

---

## 5. Recommended representation-capability metric

### 5.1 Four requirements any valid metric must satisfy

| # | Requirement | What the previous design violated |
|---|---|---|
| R1 | **Coordinate comparability** — WT and variant features produced by one function differing only in the sequence argument | WT windows centred at sequence midpoint, mutant windows at the mutation |
| R2 | **Non-degenerate reference scale** — the normaliser must be a spread over distinct, non-identical comparators | inter-protein denominator was zero by construction |
| R3 | **Capability separated from correctness** — the metric must not require the representation to predict affinity | conflated throughout |
| R4 | **Evaluation on unseen parents** — capability on seen-protein variants does not imply capability on cold parents | not addressed |

A fifth, implicit requirement: the metric must **fail** for an edit descriptor. Any metric that rewards a one-hot "V600E" token is measuring the label, not the protein.

### 5.2 The core object: the contrast vector

For a variant record `v = (parent p, position i, wt residue A, mut residue B)`, after coordinate resolution (§11.2):

```
Δφ(v) = φ_i( seq_p^WT ) − φ_i( seq_p^v )
```

Both terms evaluated at **the same resolved index** in their respective sequences. The index is a property of the *pair*, never derived from a name string. For length-changing edits, `Δφ` is undefined at the token level and the variant goes to a separate track (§5.5).

### 5.3 Primary metric: Local Selectivity, S

Following Hewitt & Liang's control-task construction, define three probes, all **low-capacity** (linear or ≤10-dimensional MLP, since larger probes memorise and regularisation aimed at generalisation gap does not restore selectivity):

| Probe | Input | Target | Symbol |
|---|---|---|---|
| Task probe | `Δφ(v)` | an **external** biological label: the ProteinGym DMS score for that exact variant, or (fallback) KLIFS pocket membership of position `i` | `A_task` |
| Control probe | `Δφ(v)` | a **random-but-consistent** label assigned per *substitution type* `(A→B)`, sampled from the empirical label distribution | `A_ctrl` |
| Identity probe | `Δφ(v)` | parent protein identity | `A_id` |

Then:

```
Local Selectivity      S = A_task − A_ctrl
Identity Leakage       L = A_id − chance
```

All three evaluated on **held-out parent proteins** with a parent-level cluster bootstrap.

**Verdict rule.** A representation is *capable* iff `S > 0` with the 95% cluster-bootstrap CI excluding 0, **and** `L` is small (pre-register a ceiling, e.g. `L ≤ 0.10` in normalised accuracy units).

**Why this construction is the right one:**

- It disqualifies the edit descriptor **automatically**. An edit descriptor is a deterministic function of the substitution type, so it can solve the control task exactly as well as the real task: `A_task = A_ctrl`, hence `S = 0`. No special-casing needed.
- It disqualifies mean-pooled global embeddings without appealing to distances: `Δφ` is dominated by nuisance variation, both probes land near chance, `S ≈ 0`.
- It has **no denominator**, so R2 is satisfied structurally rather than by inspection.
- It separates capability from correctness (R3): the task label is a *biological property of the substitution*, not the affinity change the programme is ultimately after.
- The identity probe catches the failure mode where `Δφ` encodes "which protein this is" rather than "what changed" — the local analogue of the target-ID shortcut.

### 5.4 Diagnostics (report alongside S; never as the headline)

| ID | Diagnostic | Definition | Reads out |
|---|---|---|---|
| **D1** | Matched-substitution displacement z-score (**MDZ**) — the correct replacement for the broken ratio | compare `‖Δφ(v)‖` against `‖Δφ(v')‖` for counterfactual substitutions `v'` at **other positions in the same protein**, matched on BLOSUM62 score (±1) and burial class; `MDZ = (‖Δφ(v)‖ − median) / MAD` | magnitude anomaly. Denominator is a spread over distinct real objects, so it cannot be zero by construction |
| **D2** | Locality profile | `‖φ_j(WT) − φ_j(mut)‖` as a function of `\|j − i\|` | whether the representation is actually local; a global pooling gives a flat, noise-dominated profile |
| **D3** | Aligned-position cross-parent comparability | can a probe trained on parents 1..n recover residue identity at KLIFS position k for **held-out** parents from `φ_k` alone? | whether aligned coordinates are shared across parents — a precondition for cold-parent transfer |
| **D4** | Controlled-perturbation calibration | `S` as a function of BLOSUM severity for randomly injected substitutions | a dose–response curve and an empirical detection threshold |
| **D5** | Effective rank / anisotropy of the `Δφ` cloud | participation ratio of the eigenvalue spectrum | if rank ≈ 1, the representation encodes only "something changed" |
| **D6** | Local Jacobian | finite-difference `∂φ_i/∂(one-hot at i)` over all 19 alternatives | the full local sensitivity tensor; a sanity check that the encoder is sensitive at all |

**Hierarchy: S is primary; D1–D6 are diagnostic.** Only `S` has a matched control that can fail for the right reason. D1 in particular must never be promoted to primary — it was the previous design's estimand class and it measures magnitude, not information.

### 5.5 The thirteen representation families

| Representation | Biological content | WT/mut comparable coordinates? | Non-degenerate reference? | Expected `S` | Expected `L` | Cold-parent suitability |
|---|---|---|---|---|---|---|
| Whole-sequence mean pooling | global fold/family | yes | yes (protein-pair spread) | ≈ 0 | high | poor — dominated by identity |
| Residue token at mutation site | local chemistry + context | yes | yes | **high** | low–moderate | good |
| WT−mut residue-token difference | the edit in context | yes | yes | **high** | **low** | **best** |
| Mutation-centred sequence window (one-hot) | local primary sequence | yes, if both centred at `i` | yes | moderate | low | good |
| Local PLM window (contextual) | local chemistry + evolutionary context | yes, if both centred at `i` **and** the same window bounds | yes | **high** | moderate | good |
| KLIFS aligned pocket positions | pocket chemistry in shared coordinates | **yes across parents** | yes | high | low | **best for cold parents** |
| Structure-derived pocket neighbourhood | 3-D environment | yes if a common frame exists | yes | high | low | good; needs structures |
| Contact-aware local neighbourhood | coupled positions | yes | yes | high | low | good |
| MSA / co-evolution features | constraint, coupling | yes | yes | moderate | moderate (MSA depth is family-correlated) | moderate |
| Learned variant embedding | whatever supervision put there | yes | yes | uninterpretable until probed | **potentially very high** | risky — must pass `L` |
| **Explicit mutation edit descriptor** | **none about the protein** | n/a — not a protein representation | n/a | **≈ 0 by construction** | 0 | **none** |
| Global + local decomposition | both | yes | yes | high | moderate | good; the practical default |
| Random / target-ID controls | none / identity only | n/a | n/a | ≈ 0 / ≈ 0 | 0 / maximal | none |

### 5.6 The seven governance answers the brief requires

**Historical numbering, isoforms, construct offsets.** Every variant record carries a `numbering_provenance` field with an explicit enum: `{uniprot_canonical, uniprot_isoform:<id>, construct_local, legacy_pre2004, vendor_catalogue}`. Resolution runs against a **pinned UniProt release**, with a curated historical-alias table seeded by the BRAF V599E→V600E case (a documented reference-sequence shift arising from counting the mature protein after loss of the initiator methionine) and the analogous mature-protein offsets in receptor tyrosine kinases with cleaved signal peptides — the likely source of the PDGFRα mismatches. **Hard rule: if the residue at the resolved index does not equal the annotated wild-type residue, the record fails and is quarantined. It is never silently coerced.** Four failures out of 76 was the pipeline working; the fix is to make the check mandatory and the alias table explicit, not to reduce the failure count.

**Insertions, deletions, multiple mutations, fusions, ITDs.** A typed schema with separate tracks:
`substitution` (single, token-aligned — the main track) · `multi_substitution` (all edits in one parent; contrast vector is the sum of per-position contrasts only if positions are non-adjacent, otherwise treat jointly) · `indel` (token alignment breaks; requires alignment-based contrast, e.g. profile over the aligned region) · `fusion` / `ITD` (different perturbation class entirely — **must not be pooled with substitutions**, and the Saifudeen panel contains them).

**PLM token limits and special-token indexing.** ESM prepends `<cls>` and appends `<eos>`, so a 1-based residue `r` sits at token index `r`; the off-by-one here is the single most common silent bug in this area. **Assertion: decode the token at the computed index and require it to equal the annotated wild-type residue letter.** Common ESM configurations also cap positional context near 1,022 residues. Several kinases in these panels exceed that by a wide margin — LRRK2 is ~2,527 residues and its clinically central G2019S sits far outside a 1,022-token window; MTOR and the ATM/DNA-PK family are worse. Policy required: sliding-window inference with overlap and recorded offsets, or domain-restricted sequences with the offset stored in the record. Silent truncation must raise, not warn.

**Distinguishing capability from correctness.** Two numbers, reported separately and never combined: `S` (capability — can a low-capacity probe read the perturbation out of the representation, beyond substitution-type memorisation) and, later, the affinity-change prediction accuracy (correctness). A representation may be capable and biologically wrong; it may not be incapable and biologically right.

**Detecting target-ID and family shortcuts.** The identity probe `L`; plus a family probe (can Manning group be recovered from `Δφ`?); plus the requirement that `S` survives on **held-out parents**, which no identity-based solution can do.

**Which public benchmarks validate it.** ProteinGym's substitution benchmark — ~2.7M missense variants across 217 DMS assays ✔ — is the right external validator, because it is large, labelled, and entirely independent of the kinase-affinity question. MaveDB is the upstream repository. Note the standing caveat that DMS score scales are assay-dependent and require per-assay rescaling before pooling.

**Suitability for unseen parents.** Only representations that place *aligned positions in shared coordinates across parents* (KLIFS pocket, structure-derived neighbourhoods, contact-aware windows) can support cold-parent transfer. Residue-token differences are excellent for capability but carry no guarantee of cross-parent comparability — which is exactly what diagnostic D3 measures, and why D3 should be run before committing to a representation for Stage 2.

---

## 6. Corrected planted-signal benchmark design

### 6.1 The governing principle

**Reuse the observation graph. Never reuse the values.**

Realism comes from the real panel's structure — which cells exist, which are censored, which proteins and scaffolds, what the replicate noise is. Labels are generated entirely from the synthetic model, so the truth is known exactly and no unknown real interaction contaminates the estimand.

### 6.2 Generation procedure

```
1.  Extract from the real panel (values discarded):
      G          observation mask over (protein, ligand) cells
      C_lo,C_hi  per-assay censoring thresholds
      P          protein components with their pocket sequences
      L          ligands with scaffold cluster ids
      sigma      noise SD, estimated from real replicate disagreement

2.  Build generative features that are computable for UNSEEN entities:
      x_local(p) = Z-scale encoding of KLIFS pocket positions of p
      z(l)       = ligand descriptors / fingerprint projection
    (Using real FEATURES is required. Using real LABELS is forbidden.)

3.  Main effects, actually used:
      a(p) ~ N(0, sigma_a^2),  sigma_a matched to the real per-protein mean spread
      b(l) ~ N(0, sigma_b^2),  sigma_b matched to the real per-ligand mean spread
    ASSERT: Var(a) and Var(b) recovered from the generated y match targets
            within tolerance.   <-- catches "generated but unused"

4.  Interaction truth, from LOCAL features only:
      I(p,l) = sum_r  lambda_r * <u_r, x_local(p)> * <v_r, z(l)>
      u_r supported on a sparse subset of pocket positions (locality parameter)
      I standardised to unit variance over the observed cells

5.  Labels:
      y_true = mu + a(p) + b(l) + tau * I(p,l) + eps,   eps ~ N(0, sigma^2)

6.  Censor:
      emit (y_obs, lo, hi, censor_flag) by applying C_lo / C_hi
      Both floor and ceiling for percent-activity endpoints.
```

### 6.3 What each step fixes

| Previous defect | Fix |
|---|---|
| tensor equation did not run | step 3 assertion plus a **generator round-trip test**: recover `a`, `b`, `λ`, `u`, `v` from noiseless generated data before any model exists |
| main effects generated but unused | the same assertion |
| synthetic interaction added to real endpoints | step 5 — labels are fully synthetic |
| train and eval on the same rows | splits declared in §6.6 before generation; index-disjointness asserted |
| raw predictions compared to interaction truth | the projection operator, §6.5 |
| main effects not removed from predictions | the same operator, applied identically to both sides |
| ligand-only control not a matched trained model | §6.7 — every control is the identical training procedure with one input ablated |
| endpoint offsets dominating sign accuracy | no real endpoint offsets exist |
| invalid local representation reused | gated behind §5 passing first |

### 6.4 The effect-size grid — two-dimensional, not one

Amplitude alone is the wrong axis. Detectability of a planted low-rank structure depends jointly on amplitude, **rank**, and **sparsity**, and on the aspect ratio and density of the observation graph. A one-dimensional τ grid cannot distinguish "the signal was too weak" from "the signal was too complex for the sample size."

Express amplitude in **noise-SD units**: `τ* = SD(τ·I) / σ`. This makes the grid portable across datasets and endpoints.

| Axis | Grid | Rationale |
|---|---|---|
| Amplitude `τ*` | **{0, 0.125, 0.25, 0.5, 1.0, 2.0}** | `τ* = 0` is the mandatory must-fail arm; the rest straddle the computed detection threshold |
| Rank `R` | {1, 4, 16} | detectability degrades with rank at fixed amplitude |
| Locality (driving pocket positions) | {3, 10, 85} | tests whether sparse, biologically realistic drivers are harder than dense ones |

Set the centre of the amplitude grid by computing the graph's approximate detection threshold, `τ* ≳ √( R·(n_p + n_l) / N_obs )`, for the actual mask `G`, and confirm the grid brackets it. The proposed `{0.2, 0.4, 0.8, 1.6}` is a reasonable *shape* but is unanchored to noise and silent on rank; the version above supersedes it.

### 6.5 Extracting the fitted interaction — one operator, applied to both sides

Define an ANOVA projection `Π` that removes fitted global, protein and ligand main effects, estimated by least squares **on the training cells only** (the observation graph is unbalanced, so simple row/column double-centring is not sufficient):

```
Π(M) = M − (mu_hat + a_hat(p) + b_hat(l))
```

Then evaluate `corr( Π(ŷ), Π(τ·I) )` on the held-out cells, where `Π` is **the same operator with the same fitted main effects** on both sides.

This single construction resolves three of the nine defects at once: predictions are no longer compared to interaction truth on mismatched scales; main effects are removed from the fitted side; and endpoint offsets cannot dominate sign accuracy because they are projected out of both.

**Preferred alternative where the model supports it:** have the model emit `(â, b̂, Î)` explicitly and evaluate `Î` directly. Report both routes; disagreement between them is itself a diagnostic.

### 6.6 Splits

Nested, declared before generation, asserted disjoint:

`parent protein` ⊃ `pocket component` (KLIFS pocket-identity clusters at ≤50%) ⊃ `mutation family` (variants of the same parent stay together) — crossed with `ligand scaffold cluster`. Four evaluation regimes: warm/warm, cold-ligand, cold-protein, cold-both. The truth generator must use only features available for held-out entities, or cold transfer is impossible by construction and the benchmark is unfair rather than hard.

### 6.7 Controls — each a matched trained model

| Arm | Construction | Must |
|---|---|---|
| **Oracle** | given the true generative features `u_r·x_local`, `v_r·z` | **pass** — this is the discriminator between optimisation failure and representational incapacity |
| Correct representation | the candidate | the thing under test |
| Ligand-only | identical procedure, protein input ablated | ≈ 0 |
| Target-ID | free per-protein embedding, matched dimension | ≈ 0 on cold-protein splits |
| Shuffled protein | protein→feature assignment permuted | ≈ 0 |
| Family-preserving shuffle | permuted within Manning group | ≈ 0 |
| Random capacity-matched | iid features, same dimension, **content-hash seeded** | ≈ 0 |
| `τ* = 0` | no interaction planted | ≈ 0 for **every** arm including the correct one |

**The optimisation-vs-incapacity rule:** if the oracle passes and the candidate representation fails, that is representational incapacity. If the oracle also fails, the failure is in optimisation, the loss, the split, or the harness — and nothing about the representation has been learned. Log per-branch gradient norms throughout to catch modality collapse independently.

### 6.8 Evaluation

Sign accuracy on `|Π(truth)|` above a pre-set threshold; Spearman; **slope** of the OLS regression of `Π(ŷ)` on `Π(truth)` (scale recovery; target 1.0); calibration by quantile–quantile comparison; and R² expressed as a **fraction of the oracle's R²** rather than in absolute terms.

### 6.9 Avoiding architecture favouritism

Generate truth under **three generative families** — bilinear low-rank; additive-in-positions (GAM-like); and a thresholded nonlinearity — and report per-family results. Require the candidate to pass on at least bilinear and sparse-position families. Operationally: the generator seed is sealed and the candidate model is specified in writing before the seed is revealed, or generation and modelling are done by different people. This is the cheapest available protection against tuning the model to the simulator.

### 6.10 Power at the protein-component level

Subsample `n_p ∈ {10, 20, 40, 80, 160}` protein components and plot detection rate against `n_p` at each `τ*`. The output is a **detection surface**, and the operationally important number extracted from it is the **minimum detectable `τ*` at the real panel's `n_p`** — which becomes the pre-registered minimum detectable effect for every subsequent biological null.

---

## 7. Ranked biological positive-control panels

Ranked on: biological purity (does the contrast isolate the protein change?) × statistical power (effective independent parents) × licensing practicality × cold-parent transfer relevance.

| Rank | Panel | Biological purity | Power | Licensing | Cold-parent relevance | Verdict |
|---|---|---|---|---|---|---|
| **1** | **Saifudeen 2026** — 86 approved inhibitors × 758 kinases (409 WT + 349 variants), ~290k measurements, duplicate, 1 µM at **Km ATP** ✔ | **Highest available**: WT and variants in the **same study**, same assay, ATP normalised per kinase | **Highest**: 349 variants over a large parent set | **verify** ⚠ | high — many parents enable leave-one-parent-out | **Primary.** Acquire and qualify first. |
| **2** | **KINOMEscan in-panel variants** (Davis / Karaman / PKIS2) | high — **binding** assay, so no ATP-competition distortion at all; Kd is a clean thermodynamic endpoint | moderate — ~21–27 disease-relevant variants ⚠ | mixed; PKIS2 is CC BY 4.0 ✔ | moderate | **Co-primary.** The orthogonal-assay-physics complement to #1. |
| **3** | **Duong-Ly + Anastassiadis** — 183 × 76 variants from **21 parents** ✔ | **compromised**: WT values imported from a 2011 study ✔ ⇒ every Δ carries a cross-study batch term | limited — ≤21 leave-one-parent folds | **CC BY-NC-ND 4.0** ✔, no derivative redistribution | moderate | Demote from primary. Useful as a replication cohort with an explicit batch covariate. |
| **4** | **Gatekeeper substitutions across parents** (ABL1 T315I, EGFR T790M, KIT T670I, ALK L1196M) | high — the *same positional change* in different kinases | low n, high information | mixed | **highest** — this is the direct transferability test | **The localization test** (§8), not a capability test. |
| **5** | **Stanford HIVdb genotype–phenotype** | moderate — fold-change vs WT is already a Δ, but readout is cell-based | high rows, **low effective n** (mutations co-occur) | open | low for kinases; high as an independent label system | Independent-label replication. |
| **6** | **Kinase DMS under small-molecule selection** (Persky NSMB 2020; FGFR1–4 saturation; EGFR L858R) | high for residue effects | very high positions, **very few ligands** | supplementary | high — Persky specifically shows resistance residues generalise across the kinome | Supporting evidence; **different label system, do not merge**. |
| **7** | Same-assay ortholog panels | high | very low n | mixed | moderate | Supporting. |
| **8** | GPCR mutation panels (GPCRdb) | moderate | low | open | **high** — independent superfamily replication | Stage-2 replication. |
| **9** | Antibiotic-resistance mutation × compound matrices | moderate | varies | varies | low (different fold space) | Exploratory. |

**Acquisition order:** Saifudeen (verify licence and supplementary structure first) → KINOMEscan variants → Duong-Ly/Anastassiadis with batch covariate → ProteinGym (for §5, not for §7) → HIVdb.

---

## 8. Two-level capability / localization protocol

### 8.1 The two questions, kept apart by construction

| | **Level A — protein-conditioning capability** | **Level B — biological localization** |
|---|---|---|
| Question | Can *any* valid representation predict the ligand-response change between WT and variant, or between two proteins, on an **unseen parent / pocket component**? | Can the response be attributed to the correct mutation site, aligned pocket position, or mechanism? |
| Unit of analysis | (parent, variant) pair | (parent, position) |
| Passing evidence | signed prediction of Δ on cold parents beating every control in §10 | positional attribution beating position-level controls **and** agreeing with independent mutational data |
| A global representation | **may** pass | **cannot** count as passing |
| Failure semantics | closes the tested representation on the tested surface | **does not** imply failure of Level A |

**Stated in advance:** B failing while A passes is expected and acceptable. Localization is strictly harder and is sometimes formally unidentifiable — when variants of a parent co-occur in correlated patterns (HIVdb), or when a parent contributes only one variant position, no estimator can attribute the effect to a position. Reporting a Level-B null as a failure of protein conditioning would be an error of the same class the audit has already caught twice.

### 8.2 Level-A controls

correct variant · matched-wrong variant (a different variant of the same parent, matched on substitution severity) · family-preserving protein shuffle · target-ID embedding · nearest-pocket-neighbour retrieval · capacity-matched random protein · ligand-only · majority-sign baseline · per-pair constant baseline.

The last two matter more than they look: with imbalanced Δ distributions, "always predict the majority sign" and "predict the pair's mean Δ" are strong and are frequently the true source of an apparent effect.

### 8.3 Level-B controls

correct mutation position · random position in the same protein · **BLOSUM/physicochemical-matched unrelated position** · surface residue matched by substitution severity · aligned-pocket position vs non-pocket position · and the gatekeeper cross-parent transfer test (train on some gatekeeper-bearing kinases, predict the gatekeeper effect on a held-out kinase).

### 8.4 Split units and statistics

**The split unit is the parent protein**, nested inside pocket-identity clusters at ≤50% over the 85 KLIFS positions. Ligand rows sharing a single mutation pair are **not independent** and must never be bootstrapped as such. Crossed random effects over `{parent, variant-within-parent, ligand scaffold}`; cluster bootstrap resamples **parents**, taking all their variants and ligands together. With ~21 parents (Duong-Ly) the effective N is 21; this must appear in every CI.

---

## 9. Censoring and endpoint analysis plan

### 9.1 Endpoint semantics, per platform

| Endpoint | Orientation | Natural scale | Censoring | Correct treatment |
|---|---|---|---|---|
| **Kd (KINOMEscan)** | lower = stronger | pKd (log) | blanks = Kd > 10 µM or undetected at 10 µM ✔ → **left-censored in pKd** | interval-censored likelihood; never impute pKd = 5 |
| **Ki panels (Metz)** | lower = stronger | pKi | detection limits | as above |
| **Fixed-ATP % activity** | lower % remaining = stronger | **logit** (see 9.2) | floor 0%, ceiling 100% | logit transform + double interval censoring; **stratify by inhibitor binding mode** |
| **Km-ATP % activity (Saifudeen, some RBC runs)** | lower = stronger | logit | floor/ceiling | as above, but the Cheng–Prusoff factor is constant (9.3) |
| **Apparent Kd (kinobeads)** | lower = stronger | log | detection limits | separate label system; lysate equilibrium ≠ purified enzyme |
| **Composite KIBA** | — | — | — | **excluded** — merged endpoints, no defensible harmonisation |

**No merging across rows of this table into one regression target.** Cross-platform agreement may be *measured*; it may not be *assumed*.

### 9.2 Percent activity: use the logit scale, and quantify the ceiling problem

For a single-dose measurement at concentration `C` with Hill slope 1 under competitive equilibrium:

```
fraction remaining  f = 1 / (1 + C/IC50)
  ⇒  logit(1 − f) = log C − log IC50
```

So **the logit of fractional inhibition is linear in pIC50**, and differences on the logit scale are proportional to log-potency differences. Differences in raw percentage points are not. This single transform makes percent-activity panels usable as a Δ estimand; the Hill-slope assumption should be checked against the dose–response follow-ups both Anastassiadis and Duong-Ly report ⚠.

**The ceiling problem, which is a power issue, not just a nuisance.** At 0.5–1 µM, a compound far more potent than the dose gives ~0% remaining on **both** WT and variant, so a 10-fold resistance shift can be entirely invisible. Single-dose panels have their **narrowest dynamic range exactly where WT→variant resistance contrasts live.** Required Stage-0 output: the fraction of WT/variant compound cells where at least one member lies in a responsive window (say 10–90% remaining). Cells outside it become interval bounds with known direction, or are excluded — and the exclusion rate must be reported, because it is a direct multiplier on effective sample size.

### 9.3 The ATP-concentration analysis

For an ATP-competitive inhibitor, `IC50 = Ki·(1 + [ATP]/Km)`.

- **Fixed [ATP] across kinases** → the factor varies per kinase. It contributes a protein main effect (which cancels in centred within-pair contrasts) **and** interacts with binding mode (type II, allosteric), which does **not** cancel.
- **[ATP] = Km per kinase** → the factor is exactly 2 for every ATP-competitive inhibitor, independent of kinase. Cross-kinase log-differences become Ki-comparable up to a constant. **This is the better protocol for the programme's purposes.**
- **WT vs variant under a Km-ATP protocol**: if the mutation shifts Km(ATP), the two constructs are assayed at different ATP concentrations. Self-correcting for competitive inhibitors; not for type-II or allosteric ones. **Stratify every WT→variant Δ by inhibitor binding mode before interpretation**, and record per-kinase ATP as a covariate read from each study's methods table — never inferred from the vendor name.

### 9.4 Four analyses, all reported

1. **Interval-censored likelihood** (primary) — SparseChem-style ±1/0 masks or AMPL-style ML estimation.
2. **Rank / partial-order** (co-primary) — restrict to pairs whose order is *determinate*; a value below the floor is unambiguously weaker than a potent one even when its magnitude is unknown. This endpoint survives both the Davis floor and the percent-activity ceiling.
3. **Determinate-only subset** (sensitivity) — report the induced selection bias toward potent, promiscuous compounds.
4. **Floor-imputed** (diagnostic only, never evidence) — run the field-standard imputation purely to quantify how much it changes the answer.

**Decision rule:** if 1–3 disagree in direction, the result is censoring-driven and must be reported as such.

### 9.5 The three-level separation, restated

**within-platform existence** → **cross-platform reproducibility** → **cross-platform transfer**. A failure at level 2 or 3 is weak evidence about biology and strong evidence about assay physics (§9.3). No experiment in this programme can establish universal biological absence; every negative must name the dataset, estimand, representation and platform it applies to.

---

## 10. Full control matrix

Every arm is a **matched trained model**: identical training procedure, optimiser, schedule, capacity and seeds, with exactly one input ablated or corrupted. Analytically computed "controls" are not controls.

### 10.1 Capability stage (§5)

| Arm | Construction | Expected | Detects |
|---|---|---|---|
| Task probe | `Δφ` → external biological label | above control | representation content |
| **Control task** | labels randomised per substitution *type* | at chance for a good representation | probe memorisation of substitution type |
| Identity probe | `Δφ` → parent identity | near chance | identity encoded in the contrast |
| Family probe | `Δφ` → Manning group | near chance | family shortcut |
| Edit-descriptor arm | one-hot `A→B` | `S ≈ 0` **by construction** | that the metric behaves correctly |
| Shuffled-pairing arm | `Δφ` computed between mismatched WT/variant | `S ≈ 0` | pairing integrity |
| Held-out-parent arm | probe trained on parents 1..n, tested on n+1.. | effect persists | seen-variant memorisation |

### 10.2 Planted-signal stage (§6)

Oracle (must pass) · correct representation · ligand-only · target-ID free embedding · shuffled protein · family-preserving shuffle · capacity-matched random (content-hash seeded) · `τ* = 0` arm (all arms must fail) · rank sweep · locality sweep · generative-family sweep.

### 10.3 Biological stage (§7–8)

**Level A:** correct variant · matched-wrong variant · family-preserving protein shuffle · target-ID · nearest-pocket-neighbour retrieval · capacity-matched random protein · ligand-only · majority-sign · per-pair constant.

**Level B:** correct mutation position · random position same protein · BLOSUM-matched unrelated position · severity-matched surface residue · pocket vs non-pocket position · gatekeeper cross-parent transfer.

**Data-structure controls throughout:** ligand-identity exclusion (strict arm: both members novel) · scaffold-novel subset · pocket-identity-cold at ≤50/40/30% · cross-study batch covariate (mandatory for Duong-Ly+Anastassiadis) · binding-mode stratification · responsive-window stratification · censoring strata.

### 10.4 Integrity assertions that void a run if they fail

antisymmetry of any pairwise contrast · protein-identity-zero (`Δφ(p,p) = 0` exactly) · reference-term sign reversal · every regulariser has non-zero gradient · every permutation control provably destroys the intended information · encoder is permutation-sensitive before permutation controls are run · matched evaluation rows across all arms · index-set disjointness of splits · determinism across two processes.

---

## 11. Integrity-test specification

### 11.1 Runtime assertion table

| # | Assertion | Test | Fails on |
|---|---|---|---|
| I1 | Label orientation | named-anchor check: staurosporine broadly potent; imatinib potent on ABL1/KIT/PDGFR, weak on CDK2; lapatinib EGFR/ERBB2-selective | sign flips, pKd/Kd confusion |
| I2 | Units | distribution shape; a hard floor at exactly 10⁴ nM signals imputed censoring, not measurement | µM/nM errors |
| I3 | Residue identity | resolved index residue **equals** annotated WT residue, else quarantine | BRAF-class numbering, PDGFRα offsets |
| I4 | Parent/variant sequence pairing | Hamming distance between paired sequences equals the number of annotated substitutions | mis-paired constructs |
| I5 | Seed stability | identical outputs across two OS processes | `hash()` salting (`PYTHONHASHSEED`) |
| I6 | Mutation-centred coordinate equality | WT and variant windows have identical `(start, end)` bounds | the midpoint-centring defect |
| I7 | ESM token indexing | decode token at index `r`; must equal expected residue | off-by-one from `<cls>` |
| I8 | Truncation / out-of-range | raise, never warn, if position exceeds the model's positional window | LRRK2 G2019S, MTOR, ATM |
| I9 | Contrast antisymmetry | `f(l\|p,q) + f(l\|q,p) = 0` within tolerance | ordering bugs |
| I10 | Protein-identity zero | `Δφ(p,p) = 0` exactly | representation nondeterminism |
| I11 | Reference-term sign | flipping the reference reverses the contrast sign | centring bugs |
| I12 | No test-label centring in inputs | main effects fitted on **training cells only**; assert the fitting index set excludes test | subtle label leakage through the centring term |
| I13 | Interval-bound preservation | after every transform (e.g. logit), bounds remain valid and ordered | censoring corrupted by transforms |
| I14 | Decomposition consistency | `â + b̂ + Î` reconstructs `ŷ` to numerical tolerance | broken decomposition |
| I15 | Gradient coverage | every parameter tensor receives non-zero gradient in the first N steps | dead branches, unused modules |
| I16 | Live regulariser | gradient of each penalty w.r.t. parameters is non-zero on random inputs | the softmax-L1 no-op |
| I17 | Permutation efficacy | the permutation measurably changes the encoder output | vacuous permutation controls |
| I18 | Matched evaluation rows | all arms scored on the identical row index set | silent arm-specific filtering |
| I19 | Bootstrap unit | resampling is over protein components, not rows | CI inflation by 10× or more |
| I20 | Split isolation | train/val/test index sets pairwise disjoint at every nesting level | the same-rows defect |
| I21 | Ligand/scaffold novelty | test ligands and scaffolds absent from training as declared | recall leakage |
| I22 | Generator round-trip | planted parameters recovered from noiseless generated data | untested generator code paths |
| I23 | Main-effect realisation | `Var(a)`, `Var(b)` recovered from generated `y` match targets | "generated but unused" |

### 11.2 The variant-coordinate layer

A typed record, produced by one function, consumed everywhere:

```
VariantRecord {
  parent_uniprot, uniprot_release, isoform_id
  parent_sequence, variant_sequence
  edits: [ {position_source, position_canonical, wt_aa, mut_aa} ]
  edit_type: substitution | multi_substitution | indel | fusion | ITD | phospho_state
  numbering_provenance: uniprot_canonical | uniprot_isoform | construct_local
                      | legacy_pre2004 | vendor_catalogue
  construct_start, construct_end, construct_tags
  plm_token_index, plm_window_bounds, truncation_policy
  assay: platform, atp_conc, atp_basis(fixed|km), substrate, compound_conc, replicate_n
  provenance: study_id, table_id, row_id
  qc_status: pass | quarantine, qc_reason
}
```

Every downstream artefact references a `VariantRecord` by id. Nothing parses a name string.

### 11.3 Tests that appear to pass but are vacuous

| Vacuous test | Why it cannot fail | Replacement |
|---|---|---|
| L1 penalty on softmax attention | `Σ\|aᵢ\| = Σaᵢ = 1` identically; gradient is zero | entropy penalty, L1 on pre-softmax logits, sparsemax/α-entmax, top-k, L0 gates |
| Residue permutation on a permutation-invariant encoder | mean/sum pooling is invariant by construction | assert permutation-sensitivity first (I17) |
| A permutation applied *consistently* across all proteins | a globally consistent relabelling preserves every learnable relation | permute **within** each protein, or permute the protein→feature assignment |
| Toy-only identity checks (`f(0) = 0` on a 2×2 example) | passes for broken code at real scale | property tests on real-shaped, real-sparsity inputs |
| Distance-ratio sensitivity with a same-object denominator | denominator is zero by construction | matched-substitution spread (D1) or, better, selectivity (§5.3) |
| "Control performs worse" without absolute numbers | satisfied by a model that is worse than ligand-only | require correct > ligand-only in **absolute** terms |
| Row-level bootstrap CIs | ignores within-pair correlation | protein-component cluster bootstrap |

---

## 12. Power and cluster-bootstrap plan

### 12.1 Effective sample size

```
N_eff = (# independent parent proteins / pocket components at <=50% pocket identity)
        x (# independent ligand scaffold clusters)
        restricted to determinate, in-responsive-window cells
```

Counted **after** censoring exclusions, not before. For Duong-Ly the parent term is ≤21 ✔. For Saifudeen it is materially larger and should be computed from Supplementary Table 4 ⚠.

### 12.2 Minimum requirements before any biological claim

| Quantity | Minimum |
|---|---|
| Independent parents (Level A) | ≥ 15 |
| Independent parents contributing ≥ 3 variants (Level B) | ≥ 8 |
| Ligand scaffold clusters | ≥ 20 |
| Determinate, in-window cells per cold parent | ≥ 50 |
| Power to detect parent-level ρ = 0.20 | ≥ 80% |
| Minimum detectable `τ*`, taken from §6.10 | pre-registered before the biological run |

If these are not met, the correct report is **"underpowered — inconclusive,"** and the run must not be described as a falsification.

### 12.3 Bootstrap and inference

Crossed random effects over `{parent, variant-within-parent, ligand scaffold}`; **cluster bootstrap resampling parents** (taking all their variants and ligands together), 2,000 resamples; leave-one-parent-out and leave-one-family-out jackknife on every headline number; per-parent effect distributions reported, never only the pooled value; 5 seeds with per-seed values shown; all results expressed as a fraction of the empirically estimated noise ceiling.

---

## 13. Ranked post-gate model directions

Not authorised until §14's gates pass. Scored on the brief's ten criteria (H/M/L).

| Rank | Direction | Bio. plausibility | ID-shortcut resistance | Cold-parent | k=0 | k∈{1,2,3,5} | Data need | Compute | Novelty | Interpretability | Repeat-failure risk |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **Censoring-aware interaction decomposition** (`μ + a(p) + b(l) + interaction`, with orthogonality) | H | **H** — identity is given its own named term | H | H | **H** — `a(p)` is exactly what anchors calibrate | L | L | M | **H** | **L** |
| 2 | **Aligned-pocket residue–fragment bilinear** | H | H | **H** — shared coordinates across parents | H | H | M | L | M | H | L |
| 3 | **Global + local protein representation** | H | M | H | H | H | L | L | L | M | L |
| 4 | **Mutation-aware protein operator** (protein → function over ligand space) | H | **H** if the operator subspace dim ≪ #targets | H | H | H | M | M | **H** | M | M |
| 5 | **Sparse / entmax cross-attention** | M | M | M | H | H | M | M | M | **H** | **M — this is where the softmax-L1 no-op lives** |
| 6 | **Measured-selectivity contrastive learning** | H | H — negatives matched on ligand identity | M | H | M | M | M | M | L | M |
| 7 | **Conditional neural processes** | M | **H** — the support-only vs support+protein gap as a function of k is a built-in identifiability test | M | **H** | **H** | H | M | L | L | M |
| 8 | **Task-conditioned operators** | M | H | M | H | H | H | M | H | L | M |
| 9 | **Gradient-conflict / modality-collapse control** | n/a (infrastructure) | n/a | n/a | n/a | n/a | L | L | L | M | **L — mandatory regardless** |
| 10 | **MSA / co-evolution priors** | M | M | M | M | M | M | M | L | M | M |
| 11 | **Structure / pocket priors under leakage discipline** | H | M | M | M | M | H | **H** | M | M | **H** — co-folding features can re-import PDB memorisation as leakage disguised as physics |

**Recommended post-gate sequence:** 1 → 2 → 3, with 9 running throughout as instrumentation, and 7 reserved for the k-shot stage where its support-set ablation doubles as an identifiability test.

---

## 14. Staged decision tree

```
════════════════════════════════════════════════════════════════════════
STAGE Q0 — VARIANT-COORDINATE LAYER                    (2 weeks, no model)
════════════════════════════════════════════════════════════════════════
  Build the typed VariantRecord layer (§11.2) + assertions I1-I8, I10
  Validate against EXTERNAL ground truth: resolve a ProteinGym / MaveDB
  variant set and require >=99.5% residue-identity agreement
  ── GATE Q0 ── external agreement >= 99.5%; all quarantines explained;
     BRAF-class and construct-offset aliases covered; LRRK2/MTOR-class
     truncation policy exercised.
     FAIL -> fix. Nothing downstream is interpretable without this.

════════════════════════════════════════════════════════════════════════
STAGE Q1 — REPRESENTATION CAPABILITY                   (2 weeks)
════════════════════════════════════════════════════════════════════════
  Local Selectivity S on ProteinGym (external, labelled), held-out parents
  Diagnostics D1-D6; edit-descriptor arm must score S ~ 0
  ── GATE Q1 ── at least one representation with S > 0 (parent-level CI
     excluding 0) AND L below the pre-registered ceiling AND D3 showing
     cross-parent coordinate comparability.
     FAIL -> representation research. NOT a biological conclusion.

════════════════════════════════════════════════════════════════════════
STAGE Q2 — PLANTED-SIGNAL HARNESS                      (3 weeks)
════════════════════════════════════════════════════════════════════════
  Fully synthetic labels on the real observation graph (§6)
  tau* x rank x locality grid; three generative families
  ── GATE Q2 ── oracle passes; tau*=0 arm fails for ALL arms; every
     corruption control ~ 0; minimum detectable tau* recorded and
     pre-registered as the MDE for all later nulls.
     Oracle fails -> optimisation/harness defect, not incapacity.

════════════════════════════════════════════════════════════════════════
STAGE Q3 — DATA QUALIFICATION                          (2 weeks, parallel)
════════════════════════════════════════════════════════════════════════
  Saifudeen licence + Supplementary Table 4 structure; per-kinase ATP
  and basis; compound identity by InChIKey; censoring and responsive-
  window fractions; binding-mode annotation; N_eff computation
  ── GATE Q3 ── N_eff meets §12.2; endpoint semantics documented per
     panel; no label systems merged.

════════════════════════════════════════════════════════════════════════
STAGE B1 — LEVEL-A BIOLOGICAL POSITIVE CONTROL         (3 weeks)
════════════════════════════════════════════════════════════════════════
  WT -> variant response prediction on HELD-OUT PARENTS, Saifudeen
  primary, KINOMEscan variants as orthogonal-assay replication,
  Duong-Ly+Anastassiadis with explicit batch covariate
  ── GATE B1 ── signed Delta prediction beats ligand-only, majority-sign,
     per-pair-constant, target-ID, nearest-pocket-neighbour and every
     corruption arm, in ABSOLUTE terms, on cold parents, 5 seeds,
     replicated on a second panel with a different assay physics.

════════════════════════════════════════════════════════════════════════
STAGE B2 — LEVEL-B LOCALIZATION                        (3 weeks)
════════════════════════════════════════════════════════════════════════
  Position attribution; gatekeeper cross-parent transfer
  ── GATE B2 ── beats position-level controls AND agrees with independent
     mutational data. FAIL -> proceed with the dense model; record as a
     Level-B null only. Do NOT scale in response.

════════════════════════════════════════════════════════════════════════
STAGE C — COLD-PROTEIN INTERACTION TEST                (4 weeks)
════════════════════════════════════════════════════════════════════════
  Centred protein contrasts on cold pocket components; full §10 matrix
  ── GATE C ── the actual scientific question. Only this authorises
     model-scale training.

════════════════════════════════════════════════════════════════════════
STAGE D — ZERO-SHOT AND k-SHOT DTA                     (8 weeks)
════════════════════════════════════════════════════════════════════════
  Directions 1-3 (§13); CNP support-only vs support+protein gap at
  k = 0,1,2,3,5; strong fine-tuning baselines included
  ── GATE D ── non-zero k=0 skill with controls flat AND a persistent
     (non-converging) gap out to k=5.
```

**Closure rule.** A failed gate closes only the tested {dataset × estimand × representation × platform} combination. It never closes biology. Every negative report names all four in its title.

---

## 15. Outcome definitions

### 15.1 PIPELINE-QUALIFIED

All of: Gate Q0 (≥99.5% external residue agreement) · Gate Q1 (at least one representation with `S > 0`, CI excluding 0, `L` under ceiling, D3 comparability) · Gate Q2 (oracle passes, `τ* = 0` fails everywhere, controls flat, MDE recorded) · Gate Q3 (`N_eff` sufficient, endpoint semantics documented, no merged label systems) · all §11.1 assertions passing · determinism verified across processes.

**This is a statement about instruments only.** It licenses biological experiments; it asserts nothing biological.

### 15.2 BIOLOGICALLY POSITIVE

Pipeline-qualified, **and** Gate B1 passed on held-out parents, **and** replicated on a second panel with different assay physics (binding vs activity), **and** the effect is not carried by ≤2 parent proteins under leave-one-parent-out, **and** signed alignment holds (a prediction change under protein substitution that does not align with signed truth counts as zero evidence), **and** the correct-protein model improves its own absolute predictions rather than only degrading controls, **and** 5 seeds agree.

Level-B localization is **not** required for this verdict and must be reported separately.

### 15.3 UNRESOLVED

Any of: Q0/Q1/Q2 failing (instrument defect — no biological information) · Q3 showing `N_eff` below §12.2 (underpowered) · B1 null with a failed oracle, a collapsed protein branch, a dead regulariser, a vacuous permutation control, or an effect size below the pre-registered MDE · B2 null with B1 passing (capability without localization — expected, acceptable, and **not** grounds for scaling the model) · effect confined to one parent family.

### 15.4 FALSIFIED-AS-TESTED

Requires **all** of: pipeline-qualified · MDE from Q2 at or below the effect size of interest · censoring handled by all four analyses in §9.4, agreeing · binding-mode and responsive-window stratification applied · batch covariates included where WT and variant come from different studies · no branch collapse · every §11.1 assertion passing · **and then** no Level-A signal across ≥2 panels, ≥2 representations, ≥5 seeds.

Even then the claim is bounded: *"no transferable protein-conditioned interaction signal was recoverable for {estimand} on {panels} using {representations} at {platform} under {power}."* Per §9.5, no experiment in this programme can license a claim of universal biological absence.

---

## 16. Data access, time and compute

| Stage | Duration | People | Compute | Access | Risk |
|---|---|---|---|---|---|
| Q0 coordinate layer | 2 wks | 1 | CPU | UniProt (pinned release), ProteinGym/MaveDB — open | none |
| Q1 capability | 2 wks | 1 | 1 GPU, hours (PLM inference only; probes are tiny) | ProteinGym ✔ open | none |
| Q2 planted signal | 3 wks | 1 | 1 GPU, days | none beyond already-held panels | none |
| Q3 data qualification | 2 wks (parallel) | 1 | CPU | **Saifudeen licence unverified ⚠**; Duong-Ly CC BY-NC-ND ✔ (no derivative redistribution); PKIS2 CC BY 4.0 ✔ | licence |
| B1 Level A | 3 wks | 1–2 | 1 GPU, days | as above | none |
| B2 Level B | 3 wks | 1–2 | 1 GPU | + KLIFS (free academic; verify) | none |
| C cold-protein | 4 wks | 2 | 1–2 GPUs | + GPCRdb for replication | none |
| D zero/k-shot | 8 wks | 2 | 2–4 GPUs | + FS-Mol | none |

**To a defensible pipeline-qualified verdict: ~7 weeks, one researcher, one GPU used mostly for PLM inference.** To a defensible Level-A biological verdict: ~10 weeks.

Licence constraints to plan around: Duong-Ly's **ND** clause forbids redistributing derivative datasets, so any benchmark artefact must ship **code and manifests**, not repackaged values. PDBbind's user agreement similarly prohibits derivative distribution. PKIS2 and ProteinGym are the cleanest for anything intended to be shared.

---

## 17. The single best next action

**Build and externally validate the variant-coordinate layer (Stage Q0).**

Not the capability metric. Not the planted-signal harness. Not data acquisition.

**Why this and not something more scientific-sounding:**

1. **Five of the seven capability defects are downstream of its absence.** The midpoint-vs-mutation centring, the four residue mismatches, the historical BRAF numbering, the PDGFRα offsets, and the missing governance for multi-mutation, deletion, isoform, construct-offset and truncation cases are all one missing abstraction: a typed record in which the position index is a property of the **pair**, resolved against a **pinned reference**, with **provenance**, and which **fails loudly** rather than coercing.
2. **Two of the planted-signal defects are also downstream of it**, since the harness reuses the same representation path.
3. **Without it, every subsequent experiment silently re-inherits the same defects** — including, dangerously, experiments that would then *appear* to work.
4. **It is the only stage that can be validated against external ground truth.** Resolving a ProteinGym or MaveDB variant set and requiring ≥99.5% residue-identity agreement is a genuine external check. Everything else in the qualification sequence is, to some degree, self-certified.
5. **It is cheap and bounded:** roughly two weeks, one person, no GPU, no licence negotiation, no new data.
6. **It has a hard, unambiguous pass criterion**, which is exactly what a programme recovering from an ambiguous state needs next.

**Concretely, the first week:** implement `VariantRecord` (§11.2); pin a UniProt release; build the historical-alias table seeded with the BRAF V599E→V600E case and the mature-protein/signal-peptide offsets that most likely explain the PDGFRα entries; implement assertions I3, I4, I6, I7, I8, I10; implement content-hash seeding to retire `hash()` (I5). **Second week:** resolve an external ProteinGym variant set end-to-end and measure agreement; re-resolve all 76 Duong-Ly variants and classify every quarantine by cause; exercise the truncation policy on LRRK2 (G2019S at position 2019, far outside a typical ~1,022-token PLM window) and on the MTOR/ATM-class long kinases.

**Run in parallel, at low cost:** verify the Saifudeen licence and read the structure of Supplementary Table 4. If that panel is usable, it becomes the primary positive-control surface — <cite index="285-1">86 approved inhibitors against 758 kinases, 409 wild-type and 349 oncogenic variants, screened in duplicate at 1 µM with Km ATP, ~290,000 measurements</cite> — with wild types and variants measured inside a single study, which removes the cross-study batch term that currently contaminates every Duong-Ly WT→variant difference.

**Then, and only then:** Q1 capability (§5), Q2 planted signal (§6), Q3 data qualification, and the biological stages.

---

### Verification appendix

**Verified this session ✔:** Saifudeen et al. 2026 panel dimensions, duplicate design, 1 µM at Km ATP, ~290,000 measurements, 409 WT + 349 variants (doi:10.1038/s41587-026-03090-8); Duong-Ly 183 compounds × 76 mutants from 21 cognate wild-type kinases, >13,000 pairs, duplicate measurements with 33 discrepant pairs removed, **wild-type values imported from Anastassiadis 2011 rather than measured in the same screen**, CC BY-NC-ND 4.0; Anastassiadis 178 × 300; Davis 72 × 442 with blanks denoting Kd > 10 µM or non-detection at 10 µM; PKIS2 645 compounds, CC BY 4.0; KLIFS 85-position alignment; ProteinGym 217 DMS substitution assays / ~2.7M missense variants plus indel benchmark; Hewitt & Liang selectivity and control-task construction, including the finding that dropout does not improve selectivity; BRAF V599E→V600E as a documented reference-sequence renumbering; Landrum & Riniker cross-assay noise figures.

**Requires verification before use ⚠:** Saifudeen licence and Supplementary Table 4 structure, and whether the `kirhub.fredhutch.org` portal permits bulk download; per-kinase ATP concentration and basis for Anastassiadis and Duong-Ly (read from each study's methods table, not inferred); Duong-Ly and Anastassiadis compound concentrations; the exact compound overlap between the 183 and 178 sets by standardised structure; exact censored and out-of-responsive-window fractions for every panel; the mutant and phospho-variant composition of Davis's 442 assays; Hill-slope-1 adequacy for the logit transform, checked against the dose–response follow-ups both Reaction Biology studies report.

**Hypotheses, explicitly not verified:** that Km(ATP) shifts in activating mutants create a binding-mode-dependent confound in Km-ATP protocols; that mature-protein numbering explains the PDGFRα mismatches; that the ceiling-censoring rate in the responsive window will be the dominant power constraint on single-dose panels.
