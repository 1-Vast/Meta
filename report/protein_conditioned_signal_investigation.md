# Proving and Learning a Transferable Protein-Conditioned Interaction Signal for Cold-Target DTA

**A fresh investigation. No prior architecture, dataset, estimand or conclusion is assumed correct — including my own.**

Independent senior-researcher assessment · 17 August 2026

---

## Contents

1. [Executive judgment](#1-executive-judgment)
2. [Formal problem and estimand comparison](#2-formal-problem-and-estimand-comparison)
3. [Verified primary-literature map](#3-verified-primary-literature-map)
4. [Verified public-data-source table](#4-verified-public-data-source-table)
5. [False-positive and false-negative mechanisms](#5-false-positive-and-false-negative-mechanisms)
6. [Ranked portfolio of sixteen approaches](#6-ranked-portfolio-of-sixteen-approaches)
7. [Recommended positive control and acquisition plan](#7-recommended-positive-control-and-acquisition-plan)
8. [Recommended first low-cost experiment](#8-recommended-first-low-cost-experiment)
9. [Full control matrix](#9-full-control-matrix)
10. [Censoring and statistical analysis plan](#10-censoring-and-statistical-analysis-plan)
11. [Staged research programme](#11-staged-research-programme)
12. [Solved / unresolved / falsified-as-tested](#12-solved--unresolved--falsified-as-tested)
13. [Time, compute and data-access requirements](#13-time-compute-and-data-access-requirements)
14. [Final recommendation](#14-final-recommendation)

---

## 1. Executive judgment

### 1.1 The headline

**The prior programme's unresolved result is almost certainly explained by three specific, correctable defects rather than by absence of biological signal.** I can name them, and each is testable in days:

1. **The two informative experiments were not testing the same thing.** A global protein representation showed an unexplained positive; a mutation-position-specific model failed. These are not in tension, because *a global sequence embedding is nearly blind to a point mutation*. For a protein of length L, a single substitution perturbs a mean-pooled PLM embedding by O(1/L) of the between-protein scale — typically a cosine distance around 10⁻³ against a between-protein scale of order 1. A prediction head with any reasonable Lipschitz bound **cannot** produce a 1–2 log-unit affinity change from that perturbation without becoming wildly unstable on ordinary protein pairs. The mutation model's failure is therefore predicted by representation geometry alone and carries **no biological information whatsoever**. Meanwhile the global model was never tested on the mutation task, and its positive result on cross-protein data is exactly what a target-identity shortcut also produces. Neither experiment constrains the hypothesis.

2. **The dense kinase panels were used with the field's standard censoring convention, which manufactures fake structure and destroys real structure simultaneously.** Davis et al. state plainly in the supplementary description that blank cells were *tested* but had Kd > 10 µM or were not detected in the 10 µM primary screen ([Nat. Biotechnol. 2011](https://www.nature.com/articles/nbt.1990)). The universally redistributed ML version imputes these as Kd = 10,000 nM, i.e. pKd = 5.0. On a 72 × 442 = 31,824-cell matrix where secondary literature reports roughly 9,166 detected interactions, that means **~70% of the matrix becomes a single constant**. Every difference computed inside the constant region is exactly zero; every difference crossing its boundary is a censoring artifact of known sign. This does not merely add noise — it *creates* a low-rank pseudo-interaction structure that a model can fit, and it *annihilates* the genuine interaction structure in the weak-binding regime. A programme that found panels "promising but heavily censored" was seeing this.

3. **The estimand was over-constrained on the ligand side.** Strict core/context MMP matching requires the *same transformation* to recur across many proteins. It is rare, and the brief correctly reports it did not yield a supported evaluation surface. But the ligand-side matching constraint is unnecessary: in panel data you can hold the ligand **identical**, not merely similar, and vary the protein. Matching on ligand identity is strictly stronger than MMP matching and vastly more abundant. The MMP route was the hardest available path to the easiest available contrast.

### 1.2 Judgments, with stated confidence

| # | Judgment | Confidence |
|---|---|---|
| J1 | A protein-conditioned interaction term exists biologically and is large (selectivity spans >4 log units for fixed ligands across kinases) | 0.97 |
| J2 | The prior null on mutation-position models is uninformative for the reasons in §1.1(1) | 0.90 |
| J3 | Within-platform interaction signal is recoverable from KINOMEscan-class panels after correct censoring handling | 0.80 |
| J4 | It transfers to held-out kinase pocket-identity clusters at useful effect size | 0.55 |
| J5 | It transfers across superfamilies (kinase → GPCR → protease) from public data alone | 0.15 |
| J6 | Cross-platform disagreement, if observed, will be partly *biological/assay-physical* rather than evidence of absence (see Cheng–Prusoff argument, §10.4) | 0.75 |
| J7 | The correct near-term deliverable is within-superfamily zero-shot plus cross-superfamily few-shot, not universal zero-shot | 0.80 |

### 1.3 Three corrections to my own prior analysis

The brief instructs me not to assume prior conclusions are correct. Applying that to my own previous report:

- **Wrong kill switch.** I previously proposed cross-platform reproducibility of the interaction residual as a whole-programme terminator. That was a category error. KINOMEscan measures competition against an immobilised ligand; Reaction Biology measures catalytic inhibition at fixed ATP; kinobeads measure competition in cell lysate at endogenous ATP. For ATP-competitive inhibitors, the Cheng–Prusoff relation makes an activity-assay IC50 depend on [ATP]/K_m(ATP), and **K_m(ATP) varies per kinase**, so a fixed-ATP activity panel imposes a protein-specific distortion that a binding panel does not. Cross-platform disagreement is therefore *expected* and does not falsify within-platform biology. The kill switch must be **within-platform split-half reproducibility**.
- **A mathematically vacuous penalty.** I recommended an L1 sparsity penalty over attention mass. For softmax attention, Σᵢ|aᵢ| = Σᵢ aᵢ = 1 identically — the penalty is a constant and its gradient is zero. It does nothing. Correct instruments: entropy penalty, L1 on pre-softmax logits, sparsemax/α-entmax, top-k hard selection, or L0 gates via concrete relaxation. The brief names this defect explicitly and it was present in my own text.
- **A forbidden estimator.** I recommended kernel ridge for the identification stage. The brief rules out closed-form/pseudoinverse methods as deployable candidates. Resolution adopted here: classical estimators are confined to Stage 0 *audits*; Stage 1 uses low-capacity **normally trained** models (SGD, early stopping, weight decay) from the same estimator family that will later be scaled, so that a Stage-1 pass is evidence about the deployable family.

### 1.4 What I recommend in one sentence

Build a **matched-variant, single-platform protein contrast surface** from kinase panels that already contain wild-type *and* disease-relevant mutant constructs, define the estimand as the **centred selectivity contrast** (§2.4), instrument the pipeline with a **synthetic planted-signal control** and a **representation-capability pre-check**, and only then ask whether anything transfers.

---

## 2. Formal problem and estimand comparison

### 2.1 Decomposition

For a comparable measurement `y(l,p)` (one label system, one platform):

```
y(l,p) = μ + a(p) + b(l) + I(l,p) + ε(l,p,assay,doc,replicate)
```

- `a(p)` — protein main effect: construct, expression, activation state, assay offset, K_m(ATP), publication selection.
- `b(l)` — ligand main effect: promiscuity, lipophilicity, size, aggregation, purity, stock concentration.
- `I(l,p)` — **the object of Core Task 1.**
- `ε` — measurement and comparability noise.

A ligand-only model recovers `b`. A global protein embedding used as an identifier recovers `a`. Chemical similarity recovers a smoothed `b`. All three are nuisance under this decomposition, which is why "beating them" on raw `y` is weak evidence: `a(p)` alone predicts raw `y` extremely well on warm targets and degrades gracefully on cold ones through family similarity, reproducing the exact phenomenology usually reported as success.

### 2.2 The estimand comparison table

Notation: σ = within-platform per-measurement SD (log units). "Cancels" means exactly, under the additive model.

| # | Estimand | Measures | Cancels | Shortcuts still possible | Noise SD | Censoring sensitivity | Coverage needed | Zero-shot | Few-shot | Structural need |
|---|---|---|---|---|---|---|---|---|---|---|
| **E1** | `y(l,p)` raw | everything | nothing | target ID, ligand recall, family prior — all of them | σ | moderate (single value) | any | yes but uninformative | yes | none |
| **E2** | two-way residual `y − μ̂ − â(p) − b̂(l)` | interaction | `a`,`b` (estimated, not exact) | leakage through estimated `â` on cold proteins; `â` undefined for unseen p | ≈σ (+ estimation error) | high (main effects biased by censoring) | dense-ish matrix | **no** — `â(p)` unavailable for cold p | yes (k anchors give `â`) | none |
| **E3** | `Δ_L(l₁,l₂\|p)` ligand-pair diff within protein | local SAR | `a(p)`, protein-assay offset | ligand-only predicts the *mean* Δ; family prior | √2·σ | high (both terms) | ≥2 ligands/protein | yes | yes | none |
| **E4** | `Δ_P(l\|p,q)` selectivity diff for one ligand | protein discrimination | `b(l)` **exactly**, plus all compound-level artifacts (aggregation, purity, stock error) | `a(p)−a(q)` is a constant per pair → target-ID pair effect | √2·σ | high | same ligand on 2 proteins | yes | yes | none |
| **E5** | `ΔΔ` crossed double difference | pure interaction | `a`,`b`, additive assay offsets | very few; requires 4 measurements | **2σ** | very high (4 terms) | 2×2 block | yes | yes | none |
| **E6** | strict MMP transformation effect | transformation ↔ pocket | `a(p)`, much of `b` | transformation identity if τ appears in one family | √2·σ | high | same τ on ≥3 protein clusters | yes | yes | none |
| **E7** | high-similarity ligand-pair effect (Tanimoto ≥ t) | relaxed E6 | as E6, less completely | residual ligand main effect leaks in | √2·σ | high | abundant | yes | yes | none |
| **E8** | activity-cliff *direction* (sign of E3 on cliffs) | discrete SAR ordering | magnitude nuisances | 50% baseline; class imbalance | binary | **low** — sign often survives censoring | cliff pairs per protein | yes | yes | none |
| **E9** | WT→mutant response change | residue-level causation | `b(l)`, and nearly all of `a` (same construct family, same platform) | mutant-set imbalance (few kinases carry most mutants) | √2·σ | high | WT+mutant on same platform | **limited** (needs the WT) | yes | mutation position |
| **E10** | ortholog selectivity change | species-level residue effects | `b(l)`, most of `a` | ortholog pairs are few and correlated | √2·σ | high | ortholog panels | limited | yes | alignment |
| **E11** | aligned-pocket residue substitution effect | position-specific chemistry | `b(l)`; conditions on position | position confounded with kinase identity | √2·σ | high | many kinases sharing a position contrast | yes | yes | **aligned pocket** |
| **E12** | conformational-state-conditioned effect | DFG/αC-state dependence | `b(l)` | state annotation may be ligand-derived → circular | √2·σ | high | state-annotated structures | yes | yes | **structures** |

### 2.3 What the table implies

- **E5 (ΔΔ) is the purest but the noisiest.** Four measurements give 2σ. It is the right *conceptual* definition and the wrong *operational* one when a better-conditioned equivalent exists.
- **E2 cannot support zero-shot protein transfer at all**, because `â(p)` is undefined for an unseen protein. Any programme whose primary estimand is the two-way residual has a hidden dependency on the target being warm. This is worth checking against the prior work.
- **E4 is the highest-power protein contrast**, and it cancels the entire compound-level artifact class exactly — a much bigger practical win than usually recognised, because compound aggregation, purity and stock-concentration errors are among the largest real error sources and they are *perfectly* correlated across proteins for a fixed compound.
- **E6 (strict MMP) is the hardest surface to populate**, exactly as the prior work found. It should be demoted from primary estimand to a *confirmatory* stratum.
- **E8 (cliff direction) is the censoring-robust endpoint.** When one value is censored at >10 µM and the other is 30 nM, the *sign* is still known with certainty even though the magnitude is not. This is the endpoint that survives the Davis floor.

### 2.4 Recommended estimand: the Centred Selectivity Contrast (CSC)

Define, for ligand `l` and an ordered protein pair `(p,q)` over a reference ligand set `L` measured on both:

```
CSC(l | p,q) = [y(l,p) − y(l,q)] − (1/|L|) Σ_{l'∈L} [y(l',p) − y(l',q)]
```

**Properties:**

- Cancels `b(l)` **exactly** (it is a within-ligand difference).
- Cancels `a(p) − a(q)` **exactly** by construction of the centring term — including the whole per-protein-pair assay-offset block.
- Noise SD = σ·√(2 + 2/|L|) → **√2·σ** for large `|L|`, versus 2σ for pairwise ΔΔ. **It is the ΔΔ estimand with the reference ligand replaced by the panel mean, recovering a √2 power advantage.**
- A ligand-only model predicts CSC ≡ 0 identically.
- A target-ID model can only emit a constant per protein pair, which the centring removes → predicts CSC ≡ 0 identically.
- A family-prior model can only emit a constant per family pair → also removed to first order.
- Antisymmetric: `CSC(l|p,q) = −CSC(l|q,p)`. This is a free, cheap, decisive integrity check on any model.

**Why this is the right primary estimand:** it makes the two dominant confounds *structurally inexpressible* rather than merely disfavoured, while requiring only that the same compounds be measured on both proteins on one platform — which is precisely what a profiling panel is.

**Recommended hierarchy (primary → confirmatory):**

1. **CSC on single-platform panels** (primary; §8)
2. **E9 WT→mutant** as the positive control instance of CSC where `q` differs from `p` at one residue (§7)
3. **E8 cliff direction** as the censoring-robust endpoint reported alongside every CSC result
4. **E7 relaxed-similarity ligand pairs** for the deployment-facing SAR-ordering task
5. **E6 strict MMP** as a confirmatory stratum only, never as the evaluation surface
6. **E11 position-conditioned** only after E9 passes, and only for the *localisation* question (§7.4)

---

## 3. Verified primary-literature map

Links verified during this investigation. Items marked ⚠ are from secondary sources and are flagged for Stage-0 re-verification.

### 3.1 Negative results, shortcuts, and bias

| Work | Finding relevant here | Link |
|---|---|---|
| Volkov et al., *J. Med. Chem.* 65:7946 (2022) | Explicit non-covalent interaction descriptors gave **no advantage** over ligand-only or protein-only descriptors; nearest-neighbour lookup already strong ⇒ memorisation dominates | [10.1021/acs.jmedchem.2c00487](https://pubs.acs.org/doi/abs/10.1021/acs.jmedchem.2c00487) |
| Graber et al., *Nat. Mach. Intell.* (2025) | Characterises bias/leakage in protein–ligand datasets; clean splits change conclusions | [10.1038/s42256-025-01124-5](https://www.nature.com/articles/s42256-025-01124-5) |
| Landrum & Riniker, *JCIM* 64:1560 (2024) | Minimal curation: ~65% of same-target IC50 pairs differ >0.3 log, 27% >1 log. Maximal curation: 48% / 13%, Kendall τ ≈ 0.71. **Sets the detection floor.** | [10.1021/acs.jcim.4c00049](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049) · [code](https://github.com/rinikerlab/overlapping_assays) |
| Kramer & Gedeck, *JCIM* 50:1961 (2010) | Leave-cluster-out CV required for diverse protein sets | doi:10.1021/ci100264e |
| van Tilborg, Alenicheva & Grisoni, *JCIM* 62:5938 (2022) | MoleculeACE; descriptor models beat deep models on activity cliffs. **Two published corrections** — one for early-stopping description, one for a split-labelling bug — worth noting as a reminder that cliff labelling is error-prone | [10.1021/acs.jcim.2c01073](https://pubs.acs.org/jcisd8/article/62/23/5938/852680/Exposing-the-Limitations-of-Molecular-Machine) · [code](https://github.com/molML/MoleculeACE) |
| Wallach & Heifets, *JCIM* 58:916 (2018) | AVE bias — benchmarks reward memorisation | doi:10.1021/acs.jcim.7b00403 |
| Sieg, Flachsenberg & Rarey, *JCIM* 59:947 (2019) | Bias control in structure-based VS data | doi:10.1021/acs.jcim.8b00712 |
| Lapuschkin et al., *Nat. Commun.* 10:1096 (2019) | Clever-Hans predictors — the general diagnostic frame | doi:10.1038/s41467-019-08987-4 |
| ⚠ "Clever Hans in Chemistry" preprint (2025) | Claims chemist/lab-idiom signals confound public activity benchmarks — a *document-level* confound distinct from assay noise. Preprint; treat as motivating, not established | [arXiv:2512.20924](https://arxiv.org/pdf/2512.20924) |

### 3.2 Chemogenomics / proteochemometrics — direct prior art for this exact question

| Work | Relevance | Link |
|---|---|---|
| van Westen et al., *MedChemComm* 2:16 (2011) | PCM for selectivity design and **extrapolation to novel targets** — the intellectual ancestor | [10.1039/C0MD00165A](http://xlink.rsc.org/?DOI=c0md00165a) |
| van Westen et al., *J. Cheminform.* 5:41 & 5:42 (2013) | Benchmarking 13 amino-acid descriptor sets for PCM | [part 2](https://link.springer.com/article/10.1186/1758-2946-5-42) |
| van Westen et al., HIV antivirogram PCM | **Prior positive result on exactly the WT→mutant question**: PCM predicted phenotypic log-fold-change for novel HIV mutants, recovered known *and* novel resistance mutations, 84% correct resistance classification on the Stanford set | [PMC3578754](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3578754/) |
| Cortés-Ciriano et al., *J. Cheminform.* 6:35 (2014) | Bayesian PCM; extrapolation to new targets RMSE comparable to interpolation, but **non-uniform across target space** | [PMC4083135](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4083135/) |

**This matters a great deal.** The HIV antivirogram result is an existing, published demonstration that a residue-aware protein representation predicts ligand-response *changes* across protein variants, and generalises to unseen mutants. Any claim that residue-level protein conditioning is unlearnable must explain that result away.

### 3.3 Interaction architectures

- Bai et al., **DrugBAN**, *Nat. Mach. Intell.* 5:126 (2023) — bilinear attention + conditional domain adversarial adaptation: [10.1038/s42256-022-00605-1](https://www.nature.com/articles/s42256-022-00605-1) · [arXiv](https://arxiv.org/abs/2208.02194) · [code](https://github.com/peizhenbai/DrugBAN)
- **Independent reusability report**, Xu et al., *Nat. Mach. Intell.* (2024) — reproduces and re-examines the cross-domain claims: [10.1038/s42256-024-00822-w](https://www.nature.com/articles/s42256-024-00822-w). *Read this before adopting DrugBAN.*
- Singh, Sledzieski et al., **ConPLex**, *PNAS* 120:e2220778120 (2023) — protein-anchored contrastive co-embedding; prospective kinase validation including a 1.3 nM EPHB1 binder: [10.1073/pnas.2220778120](https://www.pnas.org/doi/10.1073/pnas.2220778120) · [code](https://github.com/samsledje/ConPLex)
- Chen et al., **TransformerCPI** — includes *label-reversal* experiments, a genuine shortcut control and rare in this literature.

### 3.4 Protein representation

- Rao et al., ICLR 2021 — PLM attention maps are unsupervised contact learners: [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.12.15.422761v1) · [OpenReview](https://openreview.net/forum?id=fylclEqgvgd) · [ESM](https://github.com/facebookresearch/esm)
- Kooistra et al., **KLIFS**, *NAR* 44:D365 (2016) — 85-residue kinase pocket alignment; reported superposition RMSD 0.8 ± 0.1 Å for superposing residues, 2.2 ± 0.2 Å for the full pocket: [10.1093/nar/gkv1082](https://academic.oup.com/nar/article/44/D1/D365/2502606) · [klifs.net](https://klifs.net/) · [original, *J. Med. Chem.* 57:249 (2014)](https://pubs.acs.org/doi/10.1021/jm400378w)
- GPCRdb — generic residue numbering for the independent-superfamily replication: [gpcrdb.org](https://gpcrdb.org)

### 3.5 Mutation, resistance and protein engineering — the positive-control literature

| Work | Why it matters | Link |
|---|---|---|
| Persky et al., *Nat. Struct. Mol. Biol.* (2020), "Defining the landscape of ATP-competitive inhibitor resistance residues in protein kinases" | DMS across multiple kinases finds **generalisable residues** mediating drug resistance across the kinome; resistance mutations engineered predictively into TBK1, CSNK2A1, BRAF; a generalisable activation site confirmed in BRAF, EGFR, HER2, MEK1. **This is published evidence for transferable, residue-level, protein-conditioned interaction structure.** | [10.1038/s41594-019-0358-z](https://www.nature.com/articles/s41594-019-0358-z) |
| Duong-Ly et al., *Cell Reports* 14:772 (2016) | **183 kinase inhibitors × 76 recombinant mutant kinases**, Reaction Biology platform — a purpose-built mutant-kinase profiling matrix. Open access **CC BY-NC-ND 4.0** (note: ND restricts derivative redistribution) | [10.1016/j.celrep.2015.12.080](https://www.cell.com/cell-reports/fulltext/S2211-1247(15)01536-3) |
| FGFR1–4 saturation mutagenesis (2025/26) | All 11,520 kinase-domain point mutations × 2 FGFR inhibitors; 474 activating, 738 resistance-mediating; captured 97% of clinical acquired-resistance mutations | [PMC12807871](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12807871/) |
| EGFR L858R DMS vs osimertinib / BLU-945, *npj Precis. Oncol.* (2025) | ~17,000-variant saturation library in Ba/F3; inhibitor-specific escape spectra | [10.1038/s41698-025-01086-2](https://www.nature.com/articles/s41698-025-01086-2) |
| Stanford HIVdb genotype–phenotype sets | 7 PIs, 6 NRTIs, 3 NNRTIs; response is **log fold-change vs wild type** — already a Δ estimand | [genotype-phenotype datasets](https://hivdb.stanford.edu/pages/genopheno.dataset.html) |

**Caveat on DMS data:** the readout is fitness/enrichment under drug pressure, not affinity. It is a different label system and must not be merged with Kd/IC50. It is nonetheless the best available evidence that residue-level drug-response rules generalise across kinases.

### 3.6 Few-shot, neural processes, meta-learning

- Stanley et al., **FS-Mol**, NeurIPS D&B 2021 — 5,120 protein-target tasks, 233,786 compounds, 4,938/40/157 split: [PDF](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/file/8d3bba7425e7c98c50f52ca1b52d3735-Paper-round2.pdf)
- Chen et al., **ADKF-IFT**, ICLR 2023: [arXiv:2205.02708](https://arxiv.org/pdf/2205.02708)
- Garnelo et al., Conditional Neural Processes, ICML 2018
- Schimunek et al., **MHNfs**, ICLR 2023
- **A Strong Baseline for Molecular Few-Shot Learning** (2024) — reports that simple fine-tuning/probing baselines beat meta-learning as tasks become imbalanced: [arXiv:2404.02314](https://arxiv.org/html/2404.02314v2). Sobering and should be a required baseline.

**Key structural observation:** FS-Mol tasks are protein-defined but the standard protocol conditions on a *support set*, not on protein features. That makes it an almost ideal ablation surface — support-only vs support + protein features, as a function of k, isolates exactly the incremental value of protein information, with k = 0 as the zero-shot point.

### 3.7 Pairwise / difference learning

- Fralish & Reker, *Front. Drug Discov.* 6 (2026) — review; notes paired potency differences partially normalise systematic inter-assay bias, citing Nelen et al. 2025: [10.3389/fddsv.2026.1859068](https://www.frontiersin.org/journals/drug-discovery/articles/10.3389/fddsv.2026.1859068/full)
- Tynes et al., **PADRE**, *JCIM* (2021) — pairwise difference regression with UQ: [10.1021/acs.jcim.1c00670](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00670)
- Tyrchan & Evertsson, MMP methods review, *J. Med. Chem.*: [10.1021/acs.jmedchem.2c01787](https://pubs.acs.org/doi/10.1021/acs.jmedchem.2c01787)

### 3.8 Censoring-aware learning (implementations that exist)

- **SparseChem** (MELLODDY) — censored regression with one-sided squared loss and an explicit censoring mask (+1 upper, −1 lower, 0 none): [arXiv:2203.04676](https://arxiv.org/pdf/2203.04676)
- **AMPL** (ATOM consortium) — maximum-likelihood mean estimation for partially censored pIC50: [arXiv:2002.12541](https://arxiv.org/pdf/2002.12541) · [code](https://github.com/ATOMconsortium/AMPL)

### 3.9 Adjacent fields — how others prove a modality is conditional rather than an identifier

1. **VQA / unimodal priors.** Agrawal et al., *VQA-CP* (CVPR 2018) rebuilt splits so answer priors are **inverted** between train and test. Transfer: a **prior-inverted cold-target split** in which the family-level selectivity direction in test is opposite to train. Prior-riding collapses; residue-level chemistry does not.
2. **Multimodal collapse.** Wang, Tran & Feiszli, "What Makes Training Multi-Modal Classification Networks Hard?" (CVPR 2020); Peng et al., **OGM-GE** (CVPR 2022). Transfer: log per-branch gradient norms and per-branch generalisation gaps. **A null with a collapsed protein branch is an optimisation failure, not a falsification.**
3. **Invariant / causal learning.** IRM (Arjovsky et al.), V-REx (Krueger et al.), GroupDRO (Sagawa et al.). Transfer: environment = assay/document. But first measure NMI(document, target) — if it approaches 1, the penalty is degenerate and the method is inapplicable to that data.
4. **System identification — persistent excitation.** A response operator is identifiable only under sufficient input excitation. Transfer: the ligand × protein matrix must contain enough *shared ligands across protein clusters* or `I` is unidentifiable **in principle**, no matter the model. This must be audited before training. It is the single most neglected idea in the DTA literature and it directly addresses "was the failure data scarcity?"
5. **Operator learning / hypernetworks.** Treating the protein as producing a *function over ligand space* rather than a vector reframes the identifiability question as operator recovery, and gives a natural few-shot mechanism.

---

## 4. Verified public-data-source table

**Verification status legend:** ✔ verified this session from primary/official source · ⚠ from secondary source, re-verify at Stage 0 · ✘ not verified.

**Capability legend:** ✅ primary evidence · ⚠️ replication/secondary only · ❌ cannot support a cold-protein interaction test.

### 4.1 Single-platform profiling panels (the core resource class)

| Dataset | Dimensions | Label semantics | Direction | Platform | Censoring | Access & licence | Cold-protein test |
|---|---|---|---|---|---|---|---|
| **Davis et al. 2011** | **72 inhibitors × 442 kinase assays** ✔ (= 31,824 cells; 442 comprises ~363 distinct kinase domains **plus mutant and phospho-state variants** ⚠) | **Kd in nM** ✔ | *smaller = stronger* | KINOMEscan competition binding, single lab | **Severe & explicit**: blanks were tested but Kd > 10 µM or undetected in the 10 µM primary screen ✔. Secondary sources report ~9,166 detected interactions ⚠ ⇒ ~70% censored | Supplementary XLS to [nbt.1990](https://www.nature.com/articles/nbt.1990); redistributed widely | ✅ **primary** — *only* with censoring handled |
| **Metz et al. 2011** | **≈3,858 compounds × 172 kinases** ⚠ | pKi | *larger = stronger* (already log, inverted) | Abbott, single lab | present; fraction unverified ✘ | *Nat. Chem. Biol.* 7:200, supplementary | ✅ **primary** — best ligand diversity |
| **Anastassiadis et al. 2011** | **178 inhibitors × 300 kinases** ✔ | **% inhibition at a single concentration** — activity, not affinity | *larger = stronger inhibition* (opposite orientation to Kd) | Reaction Biology, fixed ATP | floor/ceiling at 0%/100% | *Nat. Biotechnol.* 29:1039 supplementary | ⚠️ ordering only; **never merge with Kd** |
| **Duong-Ly et al. 2016** | **183 inhibitors × 76 mutant kinases** ✔ | % remaining activity | *smaller = stronger inhibition* | **Reaction Biology — same platform as Anastassiadis** | as above | [Cell Rep. 14:772](https://www.cell.com/cell-reports/fulltext/S2211-1247(15)01536-3); **CC BY-NC-ND 4.0** ✔ (ND limits derivative redistribution) | ✅ **the positive-control resource** |
| **Karaman et al. 2008** | 38 inhibitors × 317 kinases (287 kinases, 3 lipid kinases, **27 disease-relevant mutants**) ⚠ | Kd | smaller = stronger | KINOMEscan | yes | *Nat. Biotechnol.* 26:127 supplementary | ⚠️ small but mutant-rich |
| **PKIS (Elkins et al. 2016)** | **367 inhibitors × 224 kinases + 24 GPCRs** ✔; 196 unique kinases, 5 lipid kinases, **21 disease-relevant mutants** ⚠ | % inhibition @ 1 µM (Caliper) ⚠ | larger = stronger | Caliper | floor/ceiling | [nbt.3374](https://www.nature.com/articles/nbt.3374), open data | ⚠️ |
| **PKIS2 (Drewry et al. 2017)** | **645 inhibitors** ✔ × ~392–406 kinase assays ⚠ | % control (KINOMEscan single-conc.) | smaller %ctrl = stronger | DiscoverX KINOMEscan | yes | [PLOS ONE 12:e0181585](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0181585), **CC BY 4.0**, public domain release ✔ | ✅ **primary for ordering** |
| **Klaeger et al. 2017 (Kinobeads)** | 243 clinical drugs × ~253 proteins across 4 cell-line lysates ⚠ | **apparent Kd** from lysate competition chemoproteomics — different equilibrium semantics | smaller = stronger | Kinobeads/LC-MS | yes | [Science 358:eaan4368](https://www.science.org/doi/10.1126/science.aan4368); PRIDE PXD005336 | ✅ **independent replication cohort** |
| **KCGS (Wells et al. 2021)** | 187 inhibitors × 215 kinases ⚠ | % inhibition | — | — | yes | *Int. J. Mol. Sci.* 22:566, open | ⚠️ |

> **Critical structural finding.** The KINOMEscan-family panels *already contain* wild-type and disease-relevant mutant constructs measured with the same compounds on the same platform (Karaman: 27 mutants; PKIS: 21; Davis: variants included in the 442). Combined with Duong-Ly's dedicated 183 × 76 mutant matrix on the Reaction Biology platform — which shares its platform with Anastassiadis's 178 × 300 wild-type matrix — **a same-platform, ligand-identical, WT↔mutant contrast surface already exists in the public domain at a scale of dozens of mutants × hundreds of compounds.** A "small local point-mutant panel" that proved underpowered was not this surface.

> **Critical trap.** KINOMEscan lists phosphorylated and non-phosphorylated forms of the same kinase as separate assays. These have **identical amino-acid sequences** but genuinely different affinities. Any sequence-only protein representation assigns them identical features ⇒ irreducible label noise *and* an apparent violation of the function `y = f(sequence, ligand)`. They must be excluded or handled with an explicit activation-state covariate. Autoinhibited and truncated constructs pose the same problem.

### 4.2 Aggregated literature databases

| Dataset | Scale (verified) | Label semantics | Licence | Leakage risk | Cold-protein test |
|---|---|---|---|---|---|
| **ChEMBL 36** | 17,803 drug targets, 2.8M distinct compounds ✔ (released ~Sept 2025) | IC50/Ki/Kd/EC50 — **heterogeneous** | **CC BY-SA 3.0** | document, series, chemist-idiom; missing-not-at-random | ⚠️ only with within-document pairing |
| **BindingDB** | ~3.21M binding data, 11,414 proteins, >1.41M compounds ✔ | Ki/Kd/IC50 aggregated (literature + **US patents**) | ChEMBL-derived portion **CC BY-SA 3.0**; BindingDB-curated portion **CC BY 3.0** ✔ | re-aggregation ⇒ duplicate records across nominally independent sources | ⚠️ |
| **Papyrus** | ~60M points; ~1.24M flagged high-quality exact values ⚠ | standardised/normalised ChEMBL + ExCAPE + curated sets | paper CC BY 4.0; data DOI 10.4121/16896406 ⚠ | inherits ChEMBL confounds — **standardisation ≠ harmonisation** | ⚠️ curation infrastructure |
| **Drug Target Commons** | kinase-heavy, crowd-curated ✘ | annotated, variable | check at source | curation heterogeneity | ⚠️ |
| **KIBA** | 2,111 × 229 ⚠ | **merged Ki/Kd/IC50 into one synthetic score** | redistributed widely | merged label semantics; heavily overfit | ❌ **fails the brief's harmonisation requirement — exclude from evidence** |

### 4.3 Structural resources

| Dataset | Scale | Access / licence | Notes |
|---|---|---|---|
| **PDBbind** | v2020 ~19,443 entries ⚠; now **PDBbind+ with a subscription model** ✔ | **Registration required; the user agreement prohibits distribution of derivative datasets** ✔ (documented by the HiQBind authors, who could not publish their cleaned version) | ❌ for Q1; ⚠️ structural priors only. **The redistribution ban is a real programme constraint.** |
| **PLINDER** | **449,383 PLI systems**, >500 annotations each, similarity at protein/pocket/interaction/ligand levels, paired apo + predicted structures ✔ | [github.com/plinder-org/plinder](https://github.com/plinder-org/plinder) · [plinder.sh](https://www.plinder.sh/) | ✅ **the split-construction resource**, not an affinity resource |
| **Binding MOAD / BioLiP2** | large curated complexes ✘ | open; verify | ⚠️ structural priors |
| **KLIFS** | all human/mouse kinase catalytic-domain structures; 85-position pocket alignment ✔ | [klifs.net](https://klifs.net/), free academic; verify terms | ✅ **protein representation backbone** |
| **GPCRdb** | GPCR superfamily, generic numbering | [gpcrdb.org](https://gpcrdb.org) | ✅ independent-superfamily replication |
| **SAIR** | **5,244,285 structures across 1,048,857 protein–ligand systems**, co-folded with Boltz-1x from ChEMBL + BindingDB; ~97% PoseBusters-valid ✔ | **CC BY 4.0**, free for commercial and non-commercial use; on Hugging Face ✔ | ⚠️ structures are **synthetic**; affinity labels inherit full ChEMBL/BindingDB heterogeneity; co-folding introduces a **PDB-memorisation leakage vector**. Not a solution to label semantics. |

### 4.4 Variant / resistance resources

| Dataset | Content | Label | Access | Cold-protein test |
|---|---|---|---|---|
| **Stanford HIVdb genotype–phenotype** | 7 PIs, 6 NRTIs, 3 NNRTIs across thousands of protease/RT variants ✔ | **log fold-change vs WT** — phenotypic, cell-based, already a Δ | [hivdb.stanford.edu](https://hivdb.stanford.edu/pages/genopheno.dataset.html) | ✅ **highest-purity variant existence test**; ⚠️ for localisation (mutations co-occur in correlated patterns) |
| **Persky et al. NSMB 2020 kinase DMS** | multi-kinase resistance-residue landscape | enrichment/fitness | supplementary | ✅ **cross-kinase generalisation evidence**; different label system |
| **FGFR1–4 saturation scan** | 11,520 point mutations × 2 inhibitors ✔ | selection score | supplementary | ⚠️ two ligands only |
| **EGFR L858R DMS (BLU-945/osimertinib)** | ~17,000 variants ✔ | enrichment | supplementary | ⚠️ few ligands |
| **Ortholog panels** (e.g. human/rat adenosine receptors, van Westen) | small, single-lab | Ki | ChEMBL + primary papers | ✅ ideal similarity-matched control, small n |

### 4.5 Benchmarks / task suites

| Resource | Use | Limitation |
|---|---|---|
| **FS-Mol** (5,120 tasks / 233,786 compounds ✔) | the support-set-vs-protein-features ablation, k = 0 zero-shot point | task ≡ assay confound; binary by default |
| **MoleculeACE** (30 targets) | cliff-stratum definitions | single-target ⇒ no protein contrast ❌ |
| **TDC DTI** | standard splits incl. temporal | inherits ChEMBL/BindingDB semantics |
| **Polaris / ASAP antiviral challenge** | e.g. 965 Mpro complexes; 770 SARS-CoV-2 train, 98 + 97 test ⚠ | 2 targets ⇒ ❌ for protein-cold |
| **DEL sets (BELKA-class)** | enormous ligand coverage | ~3 proteins ⇒ ❌ for protein-cold |

---

## 5. False-positive and false-negative mechanisms

Each row: the defect, whether it produces a false positive (FP) or false negative (FN), and a **specific diagnostic** that detects it.

### 5.1 Label, unit and identifier defects

| # | Defect | FP/FN | Diagnostic test |
|---|---|---|---|
| D1 | Kd vs pKd / Ki vs pKi **orientation** flipped | FP *and* FN (can invert every sign) | **Named-anchor test.** Assert known facts before modelling: staurosporine is broadly sub-nanomolar; imatinib is potent on ABL1/KIT/PDGFR and weak on CDK2; lapatinib is EGFR/ERBB2-selective. If the label orientation disagrees with three named anchors, stop. Also check that the label's extreme tail contains recognisable drug–target pairs by name. |
| D2 | Unit conversion (nM vs µM vs M) | FP/FN | Distributional sanity: Davis Kd in nM should span roughly 10⁻¹–10⁴; a median near 10³ with a floor at exactly 10⁴ signals imputed censoring, not measurement. |
| D3 | **Floor values treated as exact labels** | **FP (severe)** | Recompute the exact fraction of the matrix equal to the floor constant. If a single value occupies >20% of cells, all variance decompositions are invalid. Repeat every analysis with (a) censored likelihood, (b) rank-only, (c) uncensored subset. **Divergence across the three ⇒ censoring-driven result.** |
| D4 | Inconsistent ligand-pair ordering | FP | **Antisymmetry test**: feed (l₁,l₂) and (l₂,l₁); require prediction sums to 0 within numerical tolerance. Canonicalise ordering by InChIKey before pairing. |
| D5 | Protein/ligand identifier mis-mapping | FP/FN | Round-trip mapping (gene symbol ↔ UniProt ↔ KLIFS ID) with a fixed-point check. Explicitly enumerate **non-sequence variants**: phospho/non-phospho forms, autoinhibited and truncated constructs, cyclin-partnered complexes. Any two "proteins" with identical sequence but different labels must be flagged. |
| D6 | Duplicate compounds under different IDs | FP | Standardise to InChIKey skeleton (strip salts, isotopes, stereo where the assay cannot resolve it) and count collisions across sources. |
| D7 | Label-system merging (Kd + Ki + IC50 + % inhibition) | FP | Refuse. Where merging already happened upstream (KIBA), exclude. Report per-label-system results separately and only then compare. |

### 5.2 Statistical and evaluation defects

| # | Defect | FP/FN | Diagnostic test |
|---|---|---|---|
| D8 | Bootstrap at row level rather than **protein-component level** | **FP (severe)** | Compare CI widths from row bootstrap vs protein-cluster bootstrap. Typical inflation is an order of magnitude in effective N. Pre-register the cluster bootstrap. Fit a crossed random-effects model (p-cluster, q-cluster, ligand-scaffold) and report the variance attributable to each. |
| D9 | Repeated rows treated as independent | FP | Count distinct (protein-cluster, scaffold-cluster) cells vs rows. Report both; if rows/cells > 5, per-row metrics are meaningless. |
| D10 | Metrics dominated by one protein pair or one family | FP | **Leave-one-family-out** and **leave-one-protein-pair-out** jackknife on the headline metric. Report the full per-pair distribution, not the pooled value. If dropping any single family moves the effect >50%, the result is not general. |
| D11 | Family imbalance (TK over-represented among inhibitors) | FP | Report per-family effect sizes and a family-balanced re-weighting. |
| D12 | Mutation-pair imbalance (most mutants sit on EGFR/ABL1/KIT/FLT3) | FP | Leave-one-parent-kinase-out. A mutant-panel result driven by EGFR alone is an EGFR result. |
| D13 | **Underpowered positive control read as biological falsification** | **FN (severe)** | Pre-specify a **minimum detectable effect** from the design (§10.5). If the design cannot detect τ = 0.5 log units at 80% power, a null is uninformative *by construction* and must be reported as "inconclusive," never as "absent." |

### 5.3 Leakage defects

| # | Defect | FP/FN | Diagnostic test |
|---|---|---|---|
| D14 | Train/test **pocket** similarity despite low full-sequence identity | **FP (severe)** | Cluster on the **85 aligned KLIFS pocket positions**, not global sequence. Report the joint distribution of (global identity, pocket identity) across the split boundary. Two kinases at 25% global identity can exceed 80% pocket identity. |
| D15 | Ligand recall (same compounds recur on every protein in a panel) | FP | Enforce ligand-identity exclusion, plus a strict arm where **both** ligands of a pair are novel. Also exclude near-duplicates at ECFP4 Tanimoto ≥ 0.9 and by MCS. |
| D16 | Assay/document identity inseparable from target identity | FP (and makes IRM-type methods degenerate) | Compute **NMI(document_id, target_id)**. If NMI > 0.9, no document control is possible on that dataset ⇒ must use single-platform panels instead. |
| D17 | Structural-prior leakage via co-folding features | FP | Stratify by PDB release date relative to the structure model's training cutoff; compare predicted-contact features against a sequence-only pocket baseline on post-cutoff entries only. |
| D18 | Series leakage (a congeneric series counter-screened across targets) | FP | Cluster by Murcko scaffold × document; treat the cluster, not the molecule, as the unit. |
| D19 | Ortholog straddling | FP | Force human/rat/mouse orthologs of the same target into the same split block. |

### 5.4 Model and objective defects

| # | Defect | FP/FN | Diagnostic test |
|---|---|---|---|
| D20 | **Global embeddings acting as target IDs** | FP | **ID-equivalence test.** Replace the protein representation with a *freely learned per-target embedding* of identical dimension, trained jointly. If performance is equal or better, the biological representation is functioning as an identifier and contributes nothing beyond identity. This is decisive and cheap. |
| D21 | **Protein-branch modality collapse** | **FN (severe)** | Log per-branch gradient norms per epoch; log per-branch generalisation gap; measure output sensitivity ∂ŷ/∂(protein input) empirically. If protein-branch gradient norm decays to <5% of the ligand branch within the first epochs, the null is an optimisation artifact. Remedy: OGM-GE-style gradient modulation, PCGrad, or separate learning rates / modality-specific supervision. |
| D22 | **Wrong-protein loss "succeeding" by degrading corrupted branches** | FP | Report absolute metrics for correct, control **and** ligand-only branches. Require: correct > ligand-only in absolute terms **and** controls ≈ chance. "Correct beats corrupted while both lose to ligand-only" is a failure, not a partial success. |
| D23 | **Attention sparsity penalties that are mathematically constant** | FN (silent no-op) | For softmax attention, Σᵢ|aᵢ| = Σᵢ aᵢ = 1 identically ⇒ an L1 penalty has zero gradient. **Unit-test every regulariser: assert that its gradient w.r.t. parameters is non-zero on random inputs.** Use entropy penalties, L1 on pre-softmax logits, sparsemax/α-entmax, top-k, or L0 gates instead. |
| D24 | **Permutations that preserve all usable information** | FN (vacuous control) | If the protein encoder is permutation-invariant (mean/sum pooling, bag-of-residues), residue-order shuffling provably cannot change the output. **Before running the control, assert the encoder is permutation-sensitive**: shuffle a random protein and confirm the embedding changes by more than numerical tolerance. A control that cannot fail is not a control. |
| D25 | **Representation cannot express the target effect** | **FN (severe, and the likely cause of the prior mutation-model null)** | **Representation-capability pre-check.** For every WT/mutant pair, compute ‖x(WT) − x(mut)‖ and compare against the distribution of ‖x(p) − x(q)‖ over random protein pairs. Report the ratio. If the mutant perturbation is <1% of the between-protein scale, a Lipschitz-bounded head cannot deliver a 1–2 log-unit response without instability elsewhere — **the experiment is impossible before it is run.** Remedy: residue-local features (aligned pocket one-hot, substitution encoding, mutation-position embedding) rather than mean-pooled global vectors. |
| D26 | Capacity mismatch between correct and control arms | FP | Match parameter count, rank, and input dimensionality exactly across arms; sweep capacity and report the whole curve, not one point. |
| D27 | Cross-platform biological difference misread as implementation error | FN | See §10.4 (Cheng–Prusoff). Before "fixing" a cross-platform disagreement, test whether it is predicted by assay physics. |

### 5.5 The specific diagnosis of the prior unresolved result

The reported pattern — global representation positive, mutation-position model negative — is jointly explained by **D25 + D20 + D14**, in that order of likelihood:

- **D25** makes the mutation-model null uninformative (a point mutation barely moves a global embedding).
- **D20** makes the global-model positive uninformative (it may be identity).
- **D14** would make the global-model positive *look* stronger than it is (pocket-level leakage across a split built on global identity).

None of these is a biological finding. All three are testable in under a week, without training a large model. **This is the single most actionable conclusion in this report.**

---

## 6. Ranked portfolio of sixteen approaches

Ranking weights: P(valid evidence) × biological relevance × shortcut resistance ÷ (data risk × compute). Each entry answers the brief's twelve required questions in order.

---

### **A1 — Centred Selectivity Contrast regression on single-platform panels** ★ Rank 1

1. **Hypothesis.** For a fixed ligand, the *deviation* of its selectivity between two proteins from the panel-average selectivity of that protein pair is determined by the chemistry of the two pockets interacting with that ligand's specific features.
2. **Input/dataset.** Davis (Kd) primary; PKIS2 and Metz as separate label systems; KLIFS 85-position aligned pocket features; ECFP4 + physicochemical ligand features.
3. **Representation.** `f_θ(z(l), x(p), x(q)) → ĈSC`, with hard antisymmetry imposed architecturally, e.g. `f = g(z, x_p, x_q) − g(z, x_q, x_p)`.
4. **Objective.** Interval-censored (Tobit-style) loss on CSC where determinate; ranking loss (Somers' D surrogate) on partially determinate pairs; both reported.
5. **Cannot reduce to target ID.** Target ID emits a constant per protein pair; the centring removes it exactly. Ligand-only emits 0. Both are structurally incapable, not merely disfavoured.
6. **Controls.** Full matrix, §9. Positive: WT/mutant pairs. Negative: shuffled, family-preserving shuffled, similarity-matched wrong, random capacity-matched, residue-damaged.
7. **Falsification.** §8 — days, CPU-class.
8. **Cost.** Very low. Weeks of one researcher; small GPU or none.
9. **Failure modes.** Panel censoring; ligand recall; kinase-only scope; phospho-variant contamination.
10. **Zero-shot relevance.** High — CSC is directly the cold-target SAR-ordering quantity.
11. **k = 1,2,3,5 few-shot.** Very high: k anchors identify `a(q)`; CSC supplies the rest. **CSC is exactly the part few-shot anchors cannot supply.**
12. **Novelty.** The centring construction and its √2 power advantage over pairwise ΔΔ appear to be **novel as a stated estimand**; selectivity modelling itself is established prior art.

---

### **A2 — Mutation-conditioned response model** ★ Rank 2

1. **Hypothesis.** A single residue substitution at a pocket position produces a ligand-specific, transferable shift in affinity.
2. **Input.** Duong-Ly 183 × 76 mutants (Reaction Biology) + Anastassiadis 178 × 300 WT (same platform); KINOMEscan in-panel variants; Stanford HIVdb as an independent label system.
3. **Representation.** `Δ̂(l, p→p') = h(z(l), pos, aa_wt, aa_mut, local pocket context)`.
4. **Objective.** Censored regression on Δ plus sign classification; per-mutation random effect.
5. **Cannot reduce to target ID.** WT and mutant share >99% sequence; distinguishing them *is* the task. A target ID must be given the answer to succeed.
6. **Controls.** Randomise which position is mutated; substitute a chemically matched **distal surface residue** at matched BLOSUM distance — the ideal similarity-matched control because it is matched at the *residue*, not protein, level.
7. **Falsification.** Two-way ANOVA on the drug × variant matrix; is the interaction term significant after main effects? One afternoon.
8. **Cost.** Very low.
9. **Failure modes.** D12 (mutant imbalance), D25 (representation blindness), correlated mutation patterns in HIVdb.
10. **Zero-shot.** Moderate — needs the WT reference; but this is realistic for resistance-prediction deployment.
11. **Few-shot.** High.
12. **Novelty.** Established prior art (van Westen HIV PCM; Persky DMS). **Its value here is as the positive control, not as a novel method.**

---

### **A3 — Aligned-pocket residue vocabulary + residue–fragment cross-attention** ★ Rank 3

1. **Hypothesis.** Interaction information is localised at specific *aligned positions* (gatekeeper, hinge, back pocket), transferable across all members of an alignment.
2. **Input.** KLIFS 85-position matrix (one-hot AA + Z-scales + optional pocket-restricted ESM-2); ligand fragments/pharmacophores.
3. **Representation.** `[85 × d]` protein tokens × `[F × d]` ligand fragment tokens → cross-attention → low-rank bilinear head.
4. **Objective.** CSC loss + **entropy or α-entmax** sparsity over positions (never L1 on softmax — see D23) + antisymmetry.
5. **Cannot reduce to target ID.** Single-position input edits must produce coherent output changes matching real mutational data — an ID cannot respond to an edit.
6. **Controls.** All of A1's, plus **column-shuffled alignment** (permute the 85 positions consistently across all proteins: preserves everything except positional semantics) — the sharpest control available for this design.
7. **Falsification.** Fit an L1-regularised *linear* model on (position × AA-class) × ligand-cluster interaction terms first. No sparse informative position set ⇒ cross-attention will not find one.
8. **Cost.** 1–2 GPU-weeks.
9. **Failure modes.** Alignment quality outside the superfamily; family classification learned from the pocket string; diffuse uninterpretable attention.
10. **Zero-shot.** High within an aligned superfamily; undefined outside it.
11. **Few-shot.** High.
12. **Novelty.** Incremental architecturally; **the column-shuffle control and the entropy-sparsity correction are the contributions.**

---

### **A4 — Sparse low-rank bilinear operator with explicit main-effect decomposition** ★ Rank 4

1. **Hypothesis.** Isolating the interaction term as a *named component* prevents it absorbing main effects and makes its magnitude directly reportable.
2. **Input.** Any panel.
3. **Representation.** `ŷ = μ + a_θ(p) + b_φ(l) + ⟨U x(p), V z(l)⟩` with `U, V` low-rank and group-sparse.
4. **Objective.** Censored likelihood + orthogonality penalty forcing the bilinear output decorrelated from both main effects + conditional-information penalty (maximise `I(y; x(p) | l)`, penalise `I(x(p); target_id)`) + gradient-conflict control.
5. **Cannot reduce to target ID.** `a_θ(p)` is *given* the identity job explicitly; whatever the bilinear term adds is above and beyond it. **Reportable quantity: variance explained by the interaction component on cold targets** — a direct numerical answer to Core Task 1.
6. **Controls.** Rank sweep with capacity-matched random embeddings at every rank; branch gradient logging.
7. **Falsification.** Fit with real vs random `x(p)` at several ranks; identical cold-target interaction variance ⇒ the objective isolates nothing.
8. **Cost.** Low–medium; 3–5 GPU-days.
9. **Failure modes.** Orthogonality penalty mis-tuned; unstable adversarial CMI estimators.
10. **Zero-shot.** Very good.
11. **Few-shot.** Very good — `a(p)` is precisely what k anchors calibrate.
12. **Novelty.** The decomposition is classical (PCM/ANOVA); **the combination with conditional-information penalties and cold-target interaction-variance reporting is a genuine training innovation.**

---

### **A5 — Selectivity-profile learning (protein as query over a fixed ligand basis)** ★ Rank 5

1. **Hypothesis.** A protein's *entire* response profile over a reference ligand set is predictable from its pocket; the profile shape (not its mean) is the interaction.
2. **Input.** Dense panels with a common ligand basis.
3. **Representation.** `x(p) → ŷ(·|p) ∈ ℝ^{|L|}`, evaluated after per-protein centring (which is CSC by another route).
4. **Objective.** Profile-shape loss (cosine/Spearman after centring) + censored likelihood per entry.
5. **Cannot reduce to target ID.** Centring removes the per-protein constant; the shape is what remains.
6. **Controls.** Predict the profile of a held-out protein cluster; controls as A1.
7. **Falsification.** Nearest-pocket-neighbour profile transfer as the baseline: does a learned model beat "copy the profile of the most pocket-similar training kinase"? **If not, the model has learned similarity, not chemistry.** This baseline is essential and usually omitted.
8. **Cost.** Low.
9. **Failure modes.** Requires a common ligand basis ⇒ panel-only; profile shape dominated by promiscuity.
10. **Zero-shot.** High.
11. **Few-shot.** Very high — anchors directly fill known profile entries.
12. **Novelty.** Incremental; the nearest-pocket-neighbour baseline is the important part.

---

### **A6 — Protein-conditioned transformation learning on *relaxed-similarity* ligand pairs** ★ Rank 6

1. **Hypothesis.** A chemical change's effect depends on the pocket it contacts.
2. **Input.** ChEMBL/Papyrus with **within-document** pairing; relaxed similarity (Tanimoto ≥ 0.7 or MCS-anchored) rather than strict core/context MMP.
3. **Representation.** `(transformation encoding, pocket(p)) → Δ residual`.
4. **Objective.** Residual after subtracting the transformation's global mean effect over training proteins.
5. **Cannot reduce to target ID.** The target has expectation 0 over proteins.
6. **Controls.** τ held out; document held out; A1's negative battery.
7. **Falsification.** Count transformations recurring on ≥3 protein clusters at ≥50 instances. **The prior programme's finding that strict MMP is too sparse should be re-run at relaxed thresholds before the route is closed.**
8. **Cost.** Low.
9. **Failure modes.** Document confounding (D16); relaxation reintroduces ligand main effect.
10. **Zero-shot.** High and deployment-relevant.
11. **Few-shot.** High.
12. **Novelty.** Incremental; **the correction is relaxing the ligand-matching constraint that made the surface too sparse.**

---

### **A7 — Weakly supervised pocket discovery via multiple-instance learning** ★ Rank 7 — *high-risk / high-reward*

1. **Hypothesis.** Affinity supervision alone can identify which region governs the response, without complex coordinates.
2. **Input.** Sequences + candidate pockets (fpocket / P2Rank / AF2-derived); affinity labels only.
3. **Representation.** MIL: bag = candidate pockets; attention or noisy-OR aggregation.
4. **Objective.** Censored likelihood + MI penalty between pocket-selection distribution and target identity.
5. **Cannot reduce to target ID.** The selection variable is *within*-protein.
6. **Controls.** Random pocket; largest-pocket heuristic; **plus a pre-registered external validation**: pocket-recovery AUC against held-out annotated sites (KLIFS/BioLiP2) never used in training, threshold declared in advance (AUC ≥ 0.75).
7. **Falsification.** **Oracle-first ordering, mandatory**: on proteins with known sites, does *being given* the true pocket beat whole-sequence input? If not, discovering it cannot help.
8. **Cost.** 2–4 GPU-weeks.
9. **Failure modes.** Degenerate selection; selection tracking pocket size/hydrophobicity; **attention read as mechanism without external validation** — explicitly forbidden by the brief.
10. **Zero-shot.** Strong if it works, including on novel folds where alignments do not exist.
11. **Few-shot.** Moderate.
12. **Novelty.** Genuinely novel in this application; MIL itself is established.

---

### **A8 — MSA / co-evolution and PLM contact priors as local features** ★ Rank 8

1. **Hypothesis.** Co-evolutionary couplings and PLM attention identify pocket-lining and allosterically coupled positions, providing an unsupervised local prior.
2. **Input.** MSAs; ESM-2 attention maps (Rao et al.).
3. **Representation.** Per-position coupling/attention statistics appended to the aligned-pocket vocabulary.
4. **Objective.** As A3.
5. **Cannot reduce to target ID.** Features are positional, shared across proteins.
6. **Controls.** Shuffled MSA; PLM attention from a scrambled sequence.
7. **Falsification.** Do the prior's top positions overlap KLIFS pocket positions above chance?
8. **Cost.** Low–medium (MSA generation dominates).
9. **Failure modes.** Contacts ≠ ligand contacts; MSA depth varies by family, creating a family-correlated feature.
10. **Zero-shot.** Moderate.
11. **Few-shot.** Moderate.
12. **Novelty.** Established prior art.

---

### **A9 — Predicted structure / co-folding features under leakage discipline** ★ Rank 9

1. **Hypothesis.** Predicted complexes supply interaction hypotheses usable as features.
2. **Input.** Boltz-2/AF3-class predictions; SAIR as a pre-computed corpus (CC BY 4.0).
3. **Representation.** Predicted contact maps / interface descriptors as *inputs* to A3.
4. **Objective.** As A3. **No claim of Cartesian recognition is made or permitted.**
5. **Cannot reduce to target ID.** Features are pair-specific — *but see failure modes*.
6. **Controls.** Scrambled poses; features from a different ligand co-folded with the same protein; PDB-release-date stratification.
7. **Falsification.** On post-cutoff complexes only, do predicted-contact features beat the aligned pocket string?
8. **Cost.** **High** (co-folding inference at panel scale).
9. **Failure modes.** Serious. Co-folding models train on the PDB, so predicted features can **re-import training-set memorisation as a leakage vector disguised as physics**. An independent 2026 evaluation reports Boltz-2 affinity outputs appear largely **pose-insensitive** ([arXiv:2603.05532](https://arxiv.org/html/2603.05532v1)) — if correct, the affinity head is a learned chemogenomic prior, not a structural one.
10. **Zero-shot.** Potentially strong for novel folds; currently unvalidated for this purpose.
11. **Few-shot.** Moderate.
12. **Novelty.** Incremental; the leakage discipline is the contribution.

---

### **A10 — Contrastive learning with *measured selectivity* negatives** ★ Rank 10

1. **Hypothesis.** A metric space where "close" means "binds this pocket" encodes transferable interaction structure.
2. **Input.** Panels supply exact selectivity pairs: same ligand, strong on p, weak on q — **measured**, not assumed.
3. **Representation.** Protein-anchored embedding (ConPLex lineage) with an affinity-signed margin.
4. **Objective.** InfoNCE/triplet with margin proportional to measured |Δ|.
5. **Cannot reduce to target ID.** Negatives are matched on ligand identity; only pocket differences can satisfy the loss.
6. **Controls.** Random vs property-matched vs selectivity-mined negatives, reported separately; AVE bias score for every split.
7. **Falsification.** Linear probe on pretrained pocket embeddings — are measured selectivity pairs separable at all?
8. **Cost.** Low–medium.
9. **Failure modes.** Decoy bias; binary framing discards magnitude; missing-not-at-random false negatives.
10. **Zero-shot.** Good for retrieval; weaker for regression-grade CSC.
11. **Few-shot.** Good.
12. **Novelty.** Incremental over ConPLex; **measured rather than assumed negatives is the improvement.**

---

### **A11 — Conditional neural processes / function-space meta-learning** ★ Rank 11

1. **Hypothesis.** Protein features act as an informative prior over the per-target function, reducing the anchors needed.
2. **Input.** FS-Mol (5,120 tasks) + panels.
3. **Representation.** Support encoder → task vector `r`; decoder over `(r, x(p), z(l))`.
4. **Objective.** CNP likelihood, censoring-aware.
5. **Cannot reduce to target ID.** **The decisive ablation is `r` alone vs `r` + `x(p)` as a function of k.** A target ID is redundant with `r` once k ≳ 5, so ID-like features give converging curves; genuinely conditional information gives a *persistent* gap and non-zero skill at **k = 0**.
6. **Controls.** Shuffled protein features across tasks; assay-matched-but-different-protein tasks; k = 0 with full negative battery.
7. **Falsification.** Run the k-curve with a small deep-kernel model. If the k = 0 correct-vs-shuffled gap has CI covering 0, the representation carries no conditional information.
8. **Cost.** ~1 GPU-week (infrastructure exists).
9. **Failure modes.** Modality collapse; task ≡ assay confound; simple fine-tuning baselines may beat meta-learning ([arXiv:2404.02314](https://arxiv.org/html/2404.02314v2)) — include them.
10. **Zero-shot.** This *is* the zero-shot formulation.
11. **Few-shot.** Highest direct product relevance of any entry — k = 1,2,3,5 is native.
12. **Novelty.** Established; **the k-curve gap ablation as an identifiability test is the contribution.**

---

### **A12 — Task-conditioned operator / hypernetwork** ★ Rank 12

1. **Hypothesis.** A protein should generate a *function over ligand space*, not a vector to concatenate.
2. **Input.** As A11.
3. **Representation.** `x(p) → θ_p` (low-rank weight modulation, FiLM, or LoRA-style) parameterising a ligand scorer.
4. **Objective.** Censored likelihood + low-rank penalty on `θ_p` variation.
5. **Cannot reduce to target ID.** Constrain `θ_p` to a low-dimensional subspace with dimension ≪ number of targets, making per-target memorisation impossible by counting.
6. **Controls.** Subspace-dimension sweep; random `x(p)` at matched dimension.
7. **Falsification.** At subspace dimension 1–4, does anything survive?
8. **Cost.** Medium.
9. **Failure modes.** Optimisation instability; the subspace collapses to a promiscuity axis.
10. **Zero-shot.** High.
11. **Few-shot.** High — adapt `θ_p` from k anchors.
12. **Novelty.** Genuinely novel in DTA; established in operator learning / hypernetworks.

---

### **A13 — Invariant / causal learning across assay–document environments** ★ Rank 13

1. **Hypothesis.** Real interaction effects are invariant across assays; artifacts are not.
2. **Input.** ChEMBL with assay_id/doc_id.
3. **Representation.** Any; the contribution is the penalty.
4. **Objective.** IRM / V-REx / GroupDRO with environment = document.
5. **Cannot reduce to target ID.** Target correlates with document, so the invariance penalty actively suppresses ID solutions.
6. **Controls.** ERM baseline; random environments.
7. **Falsification.** **Measure NMI(doc_id, target_id) first.** Near 1 ⇒ degenerate ⇒ method inapplicable.
8. **Cost.** Low.
9. **Failure modes.** IRM under-performs ERM under many realistic shifts; too few environments per target.
10. **Zero-shot.** Indirect (robustness, not information).
11. **Few-shot.** Indirect.
12. **Novelty.** Established elsewhere; incremental here.

---

### **A14 — Censoring-aware, uncertainty-aware learning** ★ Rank 14 *(infrastructural, not optional)*

1. **Hypothesis.** Correct treatment of detection floors changes the conclusion.
2. **Input.** All panels.
3. **Representation.** Any.
4. **Objective.** Interval-censored likelihood (SparseChem-style masks; AMPL-style ML estimation) + heteroscedastic or deep-ensemble UQ at the **protein-cluster** level.
5. **Cannot reduce to target ID.** N/A — this is infrastructure every other entry depends on.
6. **Controls.** Censored vs rank-only vs uncensored-subset, all three reported.
7. **Falsification.** If conclusions flip across the three, the result is censoring-driven and must be reported as such.
8. **Cost.** Very low; reference implementations exist.
9. **Failure modes.** Mis-specified censoring limits; per-assay limits differ.
10–11. **Relevance.** Enables everything.
12. **Novelty.** Established prior art — **and apparently not applied in the prior programme, which is why it is listed.**

---

### **A15 — Active acquisition / prospective panel design** ★ Rank 15

1. **Hypothesis.** If public data lack excitation, the efficient move is to buy the missing cells.
2. **Input.** Commercial panel screening (KINOMEscan scanMAX = 468 kinases including mutant forms; scanEDGE = 97).
3. **Representation.** D-optimal or information-gain design over (ligand cluster × pocket cluster) cells.
4. **Objective.** Maximise expected information about `I` per assay dollar.
5. **Cannot reduce to target ID.** Design targets interaction contrasts directly.
6. **Controls.** Pre-register predictions before ordering the screen — a true prospective test.
7. **Falsification.** Cost model: a scanEDGE-scale screen of ~20 designed compounds is a modest, defined budget.
8. **Cost.** Money, not compute.
9. **Failure modes.** Cost; turnaround; needs synthesised compounds.
10–11. **Relevance.** High — a prospective result is the strongest possible evidence.
12. **Novelty.** Novel as a *designed identifiability experiment*.

---

### **A16 — A new cross-target interaction benchmark (KIN-MUT / CTIB)** ★ Rank 16 *(highest long-term leverage)*

1. **Hypothesis.** No existing benchmark can answer this question; build the one that can.
2. **Input.** Davis + Metz + PKIS2 + Duong-Ly + Anastassiadis + KLIFS + PLINDER similarity metrics.
3. **Deliverable.** Records of the form `(ligand, protein p, protein q, CSC value, censoring interval, platform, pocket-identity cluster ids, mutation annotation, scaffold cluster)`, shipped **with** the control splits (shuffled / family-preserving-shuffled / similarity-matched-wrong / random / residue-damaged) so they cannot be skipped, and with a published per-stratum noise ceiling.
4. **Objective.** N/A — infrastructure.
5–7. As above.
8. **Cost.** CPU; engineering-heavy; ~6–10 weeks.
9. **Failure modes.** Licence constraints on redistribution — **Duong-Ly is CC BY-NC-ND (no derivatives) and PDBbind forbids derivative distribution**, so the benchmark must ship *code + manifests* rather than repackaged data for those sources.
10–11. Enables everything downstream.
12. **Novelty.** Publishable; clearly distinct from PLINDER (structure/pose) and MoleculeACE (single-target cliffs).

### 6.1 Required selections

| Slot | Choice |
|---|---|
| Highest probability of valid evidence | **A1** — CSC on single-platform panels |
| Strongest biological positive control | **A2** on the Duong-Ly × Anastassiadis same-platform WT↔mutant surface, preceded by the synthetic planted-signal control (§7.1) |
| Best data-centric | **A16** — KIN-MUT / CTIB construction |
| Best representation-centric | **A3** — aligned-pocket vocabulary + cross-attention with corrected sparsity |
| Best training innovation | **A4 + A14** — censoring-aware sparse bilinear decomposition with orthogonality and conditional-information penalties |
| High-risk / high-reward | **A7** — weakly supervised pocket discovery with pre-registered external validation |
| Minimal diagnostic closing only one route | **Within-platform split-half reproducibility of CSC on Davis** (§8.4) — closes the Davis-CSC route only |
| Publishable new benchmark | **A16** |

---

## 7. Recommended positive control and acquisition plan

### 7.1 P0 — Synthetic planted-signal control (run before any real data)

**This must be step zero and it is not optional.** Generate a synthetic ligand × protein matrix with:
- realistic main effects `a(p)`, `b(l)` drawn from the empirical Davis marginals;
- a **planted** interaction term of known variance τ² generated from a known sparse function of a small set of "pocket positions";
- noise at the measured within-platform σ;
- the **same censoring pattern as Davis**, applied at the same threshold.

Then run the entire pipeline end to end and require it to recover the planted signal at the design power, across τ ∈ {0.2, 0.4, 0.8, 1.6} log units.

**Why this dominates every biological control:** it separates "the pipeline cannot detect an effect of size τ" from "biology has no effect of size τ." Without it, every null is ambiguous — which is exactly the state the prior programme is in. It also directly produces the **minimum detectable effect** needed by D13, and it exercises D3 (censoring), D8 (bootstrap level), D21 (collapse) and D23 (dead regularisers) with ground truth available.

### 7.2 Positive-control panel comparison

| Control | Ligand matching | Protein contrast | Platform consistency | Effective independent units | Answers Q-A? | Answers Q-B (localisation)? | Verdict |
|---|---|---|---|---|---|---|---|
| **Duong-Ly mutants × Anastassiadis WT** | **identical compounds** where sets overlap | single/double residue substitutions, 76 mutants | **same platform (Reaction Biology)** ✔ | mutants across many parent kinases | **Yes — strongly** | Partially (positions vary across parents) | ★ **Primary** |
| **KINOMEscan in-panel variants** (Davis / Karaman / PKIS) | identical | ~21–27 disease-relevant mutants + parents | same platform, binding assay (no Cheng–Prusoff distortion) | fewer mutants, more compounds | **Yes** | Partially | ★ **Co-primary** (binding-assay complement) |
| **Gatekeeper substitutions across kinases** (ABL1 T315I, EGFR T790M, KIT T670I, ALK L1196M) | identical | *same positional change in different kinases* | mixed | ~4–8 kinases | Yes | **Yes — this is the sharpest localisation test available** | ★ **The transferability test**: train on some gatekeeper kinases, predict the gatekeeper effect on a held-out kinase |
| **Stanford HIVdb** | fixed drug set | thousands of variants | consistent phenotype platform | high row count, **low effective n** (correlated mutation patterns) | **Yes** | **No** — co-occurrence makes single-position attribution unidentifiable | ★ Independent-label-system replication |
| **Ortholog pairs** (human/rat adenosine receptors) | identical | multi-residue, natural | single lab | small | Yes | Weak | Supporting |
| **Kinase DMS** (Persky; FGFR; EGFR L858R) | 1–2 ligands | saturation | fitness readout | very high positions, very low ligands | Yes, for *residue* effects | **Yes** | Supporting, **different label system — do not merge** |
| **Prospective synthetic (A15)** | designed | designed | one platform | designed | Yes | Yes | Best evidence, highest cost |

### 7.3 Acquisition plan

| Step | Action | Effort | Licence note |
|---|---|---|---|
| 1 | Download Davis supplementary XLS from [nbt.1990](https://www.nature.com/articles/nbt.1990); reconstruct the raw matrix **with blanks preserved as censored, not imputed** | 1 day | supplementary data |
| 2 | Download Anastassiadis (nbt.2017) and Duong-Ly (Cell Rep. 14:772) supplementary tables | 1 day | **Duong-Ly is CC BY-NC-ND 4.0 — no derivative redistribution; use in place, ship code not data** |
| 3 | Download PKIS2 from [PLOS ONE 12:e0181585](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0181585) | 1 day | CC BY 4.0, public domain release ✔ |
| 4 | Download Metz supplementary (Nat. Chem. Biol. 7:200) | 1 day | verify terms |
| 5 | Pull KLIFS 85-position alignments via API / OpenCADD-KLIFS for every kinase in the panels | 2 days | free academic; verify |
| 6 | Map all kinase assay names → UniProt → KLIFS, **explicitly annotating** mutant, phospho-state, construct and complex variants | 3–5 days | — |
| 7 | Pull Stanford HIVdb genotype–phenotype sets | 1 day | [link](https://hivdb.stanford.edu/pages/genopheno.dataset.html) |
| 8 | Optional: Klaeger kinobeads (PRIDE PXD005336) as a third platform | 2 days | — |
| 9 | **Do not** acquire KIBA. **Do not** merge label systems. | — | — |

**Total: ~2 weeks of one person.** No cost, no compute, no access negotiation.

### 7.4 Separating question A from question B

- **Question A — can any valid protein representation predict protein-to-protein ligand-response differences?** Tested by A2/A1 on the surfaces above. This is the load-bearing question for Core Task 1.
- **Question B — can the model localise the effect to the correct residue/region/mechanism?** Tested by (i) the gatekeeper cross-kinase transfer test, (ii) input-edit sensitivity, (iii) for A7, pre-registered pocket-recovery AUC against held-out annotations.

**Ruling, stated in advance:** *B failing while A passes is an expected and acceptable outcome and must not be reported as failure of protein conditioning.* Localisation is a strictly harder, and sometimes formally unidentifiable, problem — when mutations co-occur (HIVdb) or when a family has only one variant position sampled, no estimator can attribute the effect. Conversely, **A passing via a global representation is not evidence unless it also survives the cold-pair, ID-equivalence (D20), matched-control and pocket-leakage (D14) tests.**

### 7.5 Power requirements

Power must be computed on **effective independent units**, defined as:

```
N_eff = (# independent protein-pocket clusters at ≤50% pocket identity)
        × (# independent ligand scaffold clusters)
```

counted **after** removing censored-censored pairs, not as raw rows. Minimum requirements for a Stage-1 claim:

| Quantity | Minimum |
|---|---|
| Independent protein-pocket clusters | ≥ 30 |
| Independent ligand scaffold clusters | ≥ 20 |
| Determinate (non-censored-censored) CSC observations per cold cluster | ≥ 50 |
| Independent parent kinases carrying mutants (for A2) | ≥ 8 |
| Power to detect cluster-level ρ = 0.20 | ≥ 80% |
| Minimum detectable effect, declared in advance from P0 | τ ≤ 0.5 log units |

If the audit shows fewer, **the honest output is "underpowered," and Stage 1 must not be run as a falsification test.**

---

## 8. Recommended first low-cost experiment

### 8.1 Name and one-line statement

**CSC-Cold**: on the Davis binding panel with censoring handled correctly, test whether the centred selectivity contrast is predictable for protein pairs involving held-out pocket-identity clusters, using a low-capacity normally trained model, with the full control matrix.

### 8.2 Order of operations (each gate blocks the next)

```
G0  Synthetic planted-signal control (P0) ................ pipeline can detect τ ≥ 0.5?
G1  Representation-capability pre-check (D25) ............ can x(·) express the contrasts?
G2  Label/identifier audit (D1–D7) ....................... orientation, units, variants
G3  Censoring audit (D3) ................................. exact floor fraction
G4  Excitation audit (system-ID) ......................... is I identifiable in principle?
G5  Leakage audit (D14–D19) .............................. pocket clusters, NMI(doc,target)
G6  Biological positive control (§7.2) ................... WT↔mutant recovered?
G7  CSC-Cold with full control matrix .................... the actual test
```

**Steps G0–G5 involve no model training and take about three weeks.** They are where the prior programme's ambiguity gets resolved.

### 8.3 CSC-Cold specification

- **Data.** Davis (Kd, KINOMEscan). Replications on Metz (pKi) and PKIS2 (% control) as *separate* label systems.
- **Exclusions.** Phospho-state and construct variants removed or given an explicit state covariate (D5). Mutant variants **retained and labelled** — they are the positive control.
- **Protein features.** KLIFS 85-position aligned pocket: one-hot AA + Z-scale physicochemical descriptors. Second arm: pocket-restricted ESM-2 per-residue embeddings (85 positions only, never mean-pooled over the full chain).
- **Ligand features.** ECFP4 (2048, radius 2) + RDKit descriptors.
- **Estimand.** CSC per §2.4, with the reference set `L` restricted to compounds measured on both proteins.
- **Model.** Low-capacity, **normally trained**: a rank-8 factorised bilinear model or a 2-layer MLP (≤64 hidden units), SGD/Adam, weight decay, early stopping on a held-out cluster. Antisymmetry imposed architecturally. **No closed-form fits, no pseudoinverses, no test-time gradients** — so a Stage-1 pass is evidence about the deployable family.
- **Loss.** Interval-censored likelihood; plus a rank-only arm.
- **Splits.** Leave-one-Manning-group-out **and** a stricter constraint: no test kinase >50% identity over the 85 KLIFS positions to any training kinase. Two regimes: *semi-cold* (one protein of the pair cold) and *fully cold* (both cold). 5 folds × 5 seeds.
- **Ligand exclusion.** Test CSC values must involve compounds absent from training; strict arm requires this for the reference set too.

### 8.4 Primary endpoints

1. **Sign accuracy of CSC** on determinate observations with |CSC| > 1.0 log unit, correct vs each control.
2. **Spearman(ĈSC, CSC)** on cold folds, reported as a fraction of the empirical noise ceiling `r_max = τ/√(τ² + 2σ²)`.
3. **Paired correct-minus-control effect**, bootstrapped at the **protein-cluster** level.
4. **Absolute cold-target DTA improvement** over ligand-only + protein-main-effect (guards D22).
5. **ID-equivalence gap** (D20): biological representation minus free per-target embedding.

### 8.5 The minimal route-specific diagnostic

**Within-platform split-half reproducibility of CSC.** Split Davis's compound set in half; compute CSC for each protein pair independently from each half; correlate. This estimates how much of the interaction structure is reproducible *within one platform, one lab, one assay format* — removing every cross-platform confound.

- **ρ ≥ 0.4** → proceed with full confidence.
- **0.2 ≤ ρ < 0.4** → proceed with reduced scope and recalibrated effect sizes.
- **ρ < 0.2 with upper CI < 0.3** → **close the Davis-CSC route only.** This does *not* falsify biology; it says this matrix, after censoring, lacks reproducible interaction structure. Move to Metz (more ligands) and the mutant surface.

Cost: one day, no GPU. Note the contrast with my earlier proposal: this is **within-platform**, precisely because cross-platform disagreement is confounded with assay physics (§10.4).

---

## 9. Full control matrix

Every row is run, reported, and pre-registered. "Expected if signal is real" assumes the CSC estimand.

### 9.1 Baseline / capability controls

| Control | Construction | Expected | Detects |
|---|---|---|---|
| Ligand-only | drop protein features | **exactly 0 by construction** — non-zero proves a pipeline bug | pipeline integrity |
| Protein-main-effect only | per-protein constant | **exactly 0 on CSC**; large on raw y | demonstrates the confound |
| Free per-target embedding (**ID-equivalence**) | learned embedding, matched dim | correct representation must **beat** it | D20 identity shortcut |
| Nearest-pocket-neighbour profile transfer | copy most pocket-similar training kinase | model must beat it | similarity ≠ chemistry |
| Ligand-only + protein-main-effect on raw y | additive model | correct model must beat in absolute terms | D22 |

### 9.2 Protein corruption controls

| Control | Construction | Expected | Detects |
|---|---|---|---|
| Shuffled protein | permute protein→feature assignment | ≈ 0 | generic shortcut |
| **Family-preserving shuffled** | permute only *within* Manning group / family | ≈ 0 | family-prior riding |
| Similarity-matched wrong protein | substitute a paralog matched on pocket identity ±5% | substantially degraded, not 0 | fine- vs coarse-grained use |
| Random capacity-matched (A) | iid Gaussian, same dim | ≈ 0 | ID memorisation |
| Random capacity-matched (B) | real embedding matrix, rows permuted within family | ≈ 0 | taxonomy riding |
| Residue-order permutation | shuffle residue order | degraded — **only valid after asserting encoder permutation-sensitivity (D24)** | positional structure |
| Composition scramble | keep positions, randomise identities | degraded | residue chemistry |
| **Column-shuffled alignment** | permute the 85 positions consistently across all proteins | degraded | positional semantics (sharpest for A3) |
| Distal-residue substitution | mutate a matched surface residue instead of a pocket residue | ≈ 0 | residue-level specificity (A2) |

### 9.3 Data-structure controls

| Control | Construction | Expected | Detects |
|---|---|---|---|
| Ligand-identity exclusion | test ligands absent from train | effect persists | ligand recall |
| Strict exclusion | reference set also novel | effect persists, attenuated | deep recall |
| Scaffold-novel subset | Murcko scaffold held out | effect persists | scaffold recall |
| Protein-family-cold | leave-one-Manning-group-out | effect persists | family transfer |
| Pocket-identity-cold | ≤50%, ≤40%, ≤30% thresholds | monotone decay expected | graded transfer |
| Document held-out | (ChEMBL arms only) | effect persists | assay/document artifact |
| **Prior-inverted split** | invert family-level selectivity direction between train and test | effect persists, attenuated | prior riding (VQA-CP transfer) |
| Activity-cliff strata | |CSC| > 1, > 2 log units | effect **larger** on cliffs | genuine interaction vs smoothing |
| Censoring strata | fully determinate / partially / rank-only | conclusions agree | D3 |
| Temporal arm | publication-date split where available | effect persists | realism |

### 9.4 Integrity assertions (must pass or the run is void)

| Assertion | Test |
|---|---|
| Antisymmetry | `ĈSC(l\|p,q) + ĈSC(l\|q,p) = 0` within tolerance |
| Regulariser is live | gradient of every penalty w.r.t. parameters is non-zero on random inputs (D23) |
| Encoder is permutation-sensitive | shuffling residues changes the embedding (D24) |
| Representation can express the contrast | ‖x(WT) − x(mut)‖ / ‖x(p) − x(q)‖ reported (D25) |
| No branch collapse | per-branch gradient norms logged every epoch (D21) |
| Capacity matched | parameter counts identical across all arms (D26) |
| Bootstrap level correct | cluster-level, not row-level (D8) |

---

## 10. Censoring and statistical analysis plan

### 10.1 The censoring problem, precisely

Davis blanks mean "tested; Kd > 10 µM or undetected at 10 µM." These are **left-censored in pKd** (equivalently right-censored in Kd). They are *not* missing, and they are *not* pKd = 5.

Propagation to differences:

| Case | Result |
|---|---|
| both values exact | CSC term exact |
| one exact, one censored | **one-sided censored** — the difference is bounded, and its **sign is often known with certainty** |
| both censored | **completely uninformative** — must be dropped, never imputed as 0 |

The last row is the reason imputation is so destructive: imputing pKd = 5 turns every uninformative censored–censored pair into an exact zero, which is then fitted as if it were a measurement. On a matrix that is ~70% censored, this fabricates the majority of the apparent structure.

### 10.2 Four analyses, all reported

1. **Interval-censored likelihood (primary).** Tobit-style one-sided losses with an explicit mask, per SparseChem's +1/−1/0 convention ([arXiv:2203.04676](https://arxiv.org/pdf/2203.04676)); or ML mean estimation for partially censored values per AMPL ([arXiv:2002.12541](https://arxiv.org/pdf/2002.12541)).
2. **Rank-only endpoints (co-primary).** Restrict to pairs whose order is *determinate* — a >10 µM value is unambiguously weaker than a 30 nM value even though the magnitude is unknown. Report Somers' D / concordance on determinate pairs. **This is the endpoint that survives the Davis floor and it should carry equal weight to the censored likelihood.**
3. **Uncensored-subset (sensitivity).** Complete cases only. Note and report the induced bias: this selects for potent pairs and therefore for promiscuous compounds and druggable kinases.
4. **Imputed baseline (diagnostic only, never evidence).** Run the field-standard pKd = 5 imputation *purely to quantify how much it changes the answer*. The delta between (4) and (1) is a publishable number in its own right.

**Decision rule:** if (1), (2) and (3) disagree in direction, the finding is censoring-driven and must be reported as such.

### 10.3 Where ΔΔ / CSC cancels — and where it does not

**Cancels exactly:**
- additive per-protein constants (construct, expression, offset)
- additive per-ligand constants (promiscuity, aggregation, purity, stock error)
- additive per-protein-pair assay offsets (via the centring term)

**Does not cancel:**
| Residual | Mechanism |
|---|---|
| **Multiplicative scale differences** | a platform that compresses dynamic range scales *all* differences, so CSC magnitudes are not comparable across platforms even when ranks are |
| **Ligand × platform interaction** | solubility, aggregation and DMSO tolerance differ by assay format; a compound may behave differently in lysate vs purified enzyme |
| **Protein × platform interaction (the important one)** | see §10.4 |
| **Censoring** | as §10.1 — a bounded value stays bounded through the difference |
| **Non-additive activation state** | phospho vs non-phospho forms with identical sequence |

### 10.4 Cheng–Prusoff: why cross-platform disagreement is not falsification

For an ATP-competitive inhibitor in an **activity** assay, the observed IC50 depends on the ATP concentration relative to that kinase's K_m(ATP). Since K_m(ATP) varies substantially across kinases, a panel run at *fixed* ATP imposes a **protein-specific, mechanism-dependent distortion** on the measured potency. A **binding** assay (KINOMEscan competition) does not have this dependence.

Consequences that must be stated before any cross-platform comparison:
- Davis (binding) and Anastassiadis/Duong-Ly (activity at fixed ATP) are **expected** to disagree on selectivity contrasts in a systematic, kinase-specific way.
- Klaeger (lysate, endogenous ATP, cellular context) differs again.
- Therefore **a cross-platform null is weak evidence about biology and strong evidence about assay physics.**
- Conversely, cross-platform *agreement* is unusually strong evidence, because it must survive this distortion.

### 10.5 The four-level distinction (mandated by the brief)

| Level | Claim | Test | What a null means |
|---|---|---|---|
| **L1** | within-platform existence | split-half reproducibility of CSC on one platform | this matrix lacks reproducible interaction structure after censoring — closes one route |
| **L2** | cross-platform reproducibility | CSC correlation on the shared compound × kinase overlap | platforms measure different observables — expected, see §10.4 |
| **L3** | cross-platform model transfer | train on platform 1, test on platform 2 | representation is platform-entangled — a calibration problem |
| **L4** | universal biological absence | **not establishable from this evidence base under any outcome** | — |

**No experiment in this programme can support L4.** Any negative conclusion must be phrased as "falsified as tested on {dataset, estimand, representation, platform}."

### 10.6 Statistical model

Fit crossed random effects:

```
CSC_ijk = f_θ(z_i, x_j, x_k) + u_j + u_k + v_i + e_ijk
```

with `u` over protein-pocket clusters, `v` over ligand scaffold clusters. Report:

- variance components with CIs;
- **cluster bootstrap** resampling protein clusters (taking all pairs among sampled clusters) — this, not the row bootstrap, defines significance;
- leave-one-family-out and leave-one-parent-kinase-out jackknife on every headline number;
- per-protein-pair distributions, never only pooled values;
- 5 seeds for confirmation, per-seed values shown;
- all results as a **fraction of the noise ceiling**, with the ceiling estimated empirically from split-half reproducibility rather than assumed.

---

## 11. Staged research programme

```
════════════════════════════════════════════════════════════════════
STAGE 0 — DATA, SEMANTICS, CENSORING, EXCITATION, POWER   (3 weeks, no model)
════════════════════════════════════════════════════════════════════
 0.1  Acquire panels (§7.3), blanks preserved as censored
 0.2  Label orientation + unit audit (D1, D2) via named anchors
 0.3  Identifier audit (D5): enumerate mutant / phospho / construct variants
 0.4  Censoring audit (D3): exact floor fraction per dataset
 0.5  Excitation audit: shared ligands across pocket clusters
 0.6  Leakage audit (D14–D19): pocket clustering; NMI(doc, target)
 0.7  Noise ceiling: within-platform split-half reproducibility
 0.8  Power: N_eff and minimum detectable effect
 ── GATE 0 ── ≥30 pocket clusters · ≥20 scaffold clusters · ≥50 determinate
    CSC obs per cold cluster · split-half ρ ≥ 0.2 · power ≥80% for ρ=0.20
    FAIL → the tested dataset is closed; move to the next panel or to A15
           (acquisition). NOT a biological conclusion.

════════════════════════════════════════════════════════════════════
STAGE 0-P — BIOLOGICAL POSITIVE CONTROL                    (2 weeks)
════════════════════════════════════════════════════════════════════
 0-P.0  SYNTHETIC PLANTED-SIGNAL control (§7.1) at τ ∈ {0.2,0.4,0.8,1.6}
 0-P.1  REPRESENTATION-CAPABILITY pre-check (D25): ‖x(WT)−x(mut)‖ ratios
 0-P.2  WT↔mutant CSC on Duong-Ly × Anastassiadis (same platform)
 0-P.3  WT↔mutant CSC on KINOMEscan in-panel variants (binding assay)
 0-P.4  Gatekeeper cross-kinase transfer (train on some, predict held-out)
 0-P.5  Stanford HIVdb replication (independent label system)
 ── GATE 0-P ── planted signal recovered at design power · representation
    demonstrably able to express the WT/mutant contrast · known
    protein-conditioned effects recovered with controls flat
    FAIL at 0-P.0/0-P.1 → PIPELINE or REPRESENTATION defect. Fix. Do not
    interpret any downstream null.
    FAIL at 0-P.2–0-P.5 with 0-P.0/0-P.1 passing → informative negative on
    residue-level conditioning for the tested representation only.

════════════════════════════════════════════════════════════════════
STAGE 1 — TRANSFERABLE-SIGNAL IDENTIFICATION   (4 weeks, low-capacity,
                                                normally trained)
════════════════════════════════════════════════════════════════════
 1.1  CSC-Cold (§8) — semi-cold and fully cold regimes
 1.2  Complete control matrix (§9), including ID-equivalence and
      nearest-pocket-neighbour baselines
 1.3  Censoring triple analysis (§10.2)
 1.4  Replication on a second label system (Metz or PKIS2)
 ── GATE 1 ── ALL of:
    (a) cold-fold Spearman ≥ 0.25, and ≥ 0.15 for every one of 5 seeds
    (b) sign accuracy ≥ 60% on |CSC| > 1 log unit
    (c) every corruption control within [−0.05, +0.05], CI covering 0
    (d) correct − control ≥ 0.20, cluster-bootstrap CI excluding 0
    (e) correct model beats ligand-only + main-effect in ABSOLUTE terms
    (f) ID-equivalence gap > 0 (biology beats a free per-target embedding)
    (g) beats nearest-pocket-neighbour transfer
    (h) not carried by ≤2 families (drop-one-family jackknife)
    (i) replicates in a second label system at ≥60% of primary effect
    (j) all §9.4 integrity assertions pass
    ── ONLY GATE 1 AUTHORISES STAGE 2 ──

════════════════════════════════════════════════════════════════════
STAGE 2 — LOCAL INTERACTION REPRESENTATION                (6 weeks)
════════════════════════════════════════════════════════════════════
 2.1  A3 aligned-pocket cross-attention (corrected sparsity)
 2.2  Column-shuffled-alignment control
 2.3  Positional sensitivity vs known determinants (external validation)
 2.4  A7 oracle-pocket test before any discovery attempt
 ── GATE 2 ── a sparse, chemically interpretable position set beats the
    dense baseline AND positional sensitivity agrees with independent
    mutational data. FAIL → signal exists but is not localisable with this
    representation; proceed to Stage 3 with the dense model, do not scale.

════════════════════════════════════════════════════════════════════
STAGE 3 — ONE CONCENTRATED TRAINING INNOVATION            (6 weeks)
════════════════════════════════════════════════════════════════════
 3.1  A4 + A14: censoring-aware sparse bilinear decomposition with
      orthogonality and conditional-information penalties
 3.2  Attribution: ablate each mechanism separately; the gain must be
      attributable to ONE named mechanism
 3.3  Gradient-conflict / collapse monitoring throughout
 ── GATE 3 ── cold-target RMSE improvement ≥ 0.10 log units over the
    Stage-1 model, cluster-bootstrap CI excluding 0, 5 seeds; improvement
    LARGER on the activity-cliff stratum than off it; controls flat.

════════════════════════════════════════════════════════════════════
STAGE 4 — ZERO-SHOT AND FEW-SHOT DTA                      (8 weeks)
════════════════════════════════════════════════════════════════════
 4.1  A11 CNP / A12 operator with k = 0,1,2,3,5,16
 4.2  Support-only vs support + protein features at every k
 4.3  Strong fine-tuning/probing baselines included (they often win)
 4.4  Prospective or temporal holdout
 ── GATE 4 ── non-zero k=0 skill with controls flat, AND a persistent
    (non-converging) support-only vs support+protein gap out to k = 5.
```

**Stage-closure rule (mandated).** A failed stage closes only the tested {dataset × estimand × representation × platform} combination. It never closes biology. Every negative report must name all four coordinates in its title.

---

## 12. Solved / unresolved / falsified-as-tested

### 12.1 SOLVED

All of the following, simultaneously:

1. **Existence.** Within-platform split-half CSC reproducibility ρ ≥ 0.4 on ≥1 panel, with ≥50 determinate observations per cold cluster.
2. **Positive control.** Synthetic planted signal recovered at design power; representation-capability pre-check passed; known WT→mutant effects recovered with controls flat.
3. **Transfer.** Gate 1 (a)–(j) passed on pocket-identity-cold splits at ≤50% threshold.
4. **Signed alignment.** Sign accuracy ≥ 60% on |CSC| > 1 log unit. *A large prediction change under protein perturbation that does not align with signed differences counts as zero evidence.*
5. **Self-improvement.** The correct-protein model beats ligand-only + protein-main-effect in **absolute** terms. "Correct beats corrupted while both lose to ligand-only" is a failure.
6. **Not identity.** ID-equivalence gap > 0 and nearest-pocket-neighbour transfer beaten.
7. **Breadth.** Survives drop-one-family; replicates in a second label system; ideally a second superfamily (kinase + GPCR via GPCRdb).
8. **Robustness.** Persists under censored, rank-only and uncensored analyses; persists on the prior-inverted split; 5 seeds.
9. **Utility.** Gate 3 met — it actually improves cold-target DTA, with cliff-stratum gains exceeding non-cliff gains.

### 12.2 UNRESOLVED (the most probable honest outcome)

Any of:
- Gate 0 or Gate 0-P.0/0-P.1 fails → **pipeline or representation defect**; nothing about biology has been learned. *This is where the prior programme currently sits.*
- Gate 1 passes but Gate 2 fails → signal exists, is **not localisable** with current representations. Response: representation research. **Explicitly not: a larger model.**
- Effect confined to one superfamily → restate as **family-conditional transferability** and re-scope the product.
- Effect present but below the noise ceiling for practical use → a **data-generation** problem (A15), not a modelling problem.
- Any null obtained with a failed positive control, insufficient excitation, mis-handled censoring, a collapsed protein branch, a dead regulariser, a vacuous permutation control, or an expressively incapable representation → **not a falsification**; the experiment did not test the hypothesis.

### 12.3 FALSIFIED-AS-TESTED

Requires **all** of:
1. Synthetic planted-signal control passed at τ ≤ 0.5 (the pipeline demonstrably detects effects of the relevant size);
2. Representation-capability pre-check passed (the representation can express the contrast);
3. Gate 0 excitation and power requirements met;
4. Censoring handled by all three analyses, agreeing;
5. No branch collapse; all §9.4 integrity assertions passed;
6. Label orientation and identifier mapping verified against named anchors;
7. Biological positive control **passed** (proving the pipeline recovers known protein-conditioned biology);
8. And *then* no cold-fold CSC signal across ≥2 datasets, ≥2 protein representations, ≥5 seeds.

**Even then the claim is bounded:** *"no transferable protein-conditioned interaction signal was recoverable for {estimand} on {datasets} using {representations} at {platform} under {power}."* That is a statement about a tested combination. It is not a statement about biology, and per §10.5 no experiment in this programme can license one.

---

## 13. Time, compute and data-access requirements

| Stage | Duration | People | Compute | Data access | Cost risk |
|---|---|---|---|---|---|
| 0 | 3 weeks | 1 | CPU only | all open; supplementary downloads | none |
| 0-P | 2 weeks | 1 | CPU, optional 1 GPU | Duong-Ly (**CC BY-NC-ND — no derivative redistribution**); HIVdb open | licence care |
| 1 | 4 weeks | 1–2 | 1 GPU, days | as above | none |
| 2 | 6 weeks | 2 | 1–2 GPUs, 1–2 weeks | + KLIFS, GPCRdb | none |
| 3 | 6 weeks | 2 | 2–4 GPUs | as above | none |
| 4 | 8 weeks | 2 | 2–4 GPUs | + FS-Mol | none |
| A16 benchmark | 6–10 weeks, parallel | 1 | CPU | ship **code + manifests**, not repackaged data, for NC-ND and PDBbind-derived content | licence |
| A15 prospective (optional) | 8–12 weeks | 1 + chemistry | — | commercial panel screening | **budget** |

**To a defensible Gate-1 decision: ~9 weeks, one researcher, essentially no compute cost.**

Access constraints to plan around:
- **PDBbind**: registration required; user agreement **prohibits distributing derivative datasets**; now on a subscription model (PDBbind+). Avoid making it load-bearing.
- **Duong-Ly (Cell Reports)**: CC BY-NC-ND 4.0 — use in place; do not redistribute derivatives.
- **ChEMBL**: CC BY-SA 3.0 — share-alike propagates to derived datasets.
- **BindingDB**: CC BY-SA 3.0 for ChEMBL-derived records, CC BY 3.0 for BindingDB-curated records — mixed licensing within one file.
- **PKIS2** (PLOS ONE) and **SAIR** (CC BY 4.0): cleanest terms for a redistributable benchmark.

---

## 14. Final recommendation

**The next action should be pipeline correction, executed as a hard two-to-three-week prerequisite, immediately followed by data acquisition of the matched-variant panel surface. It should not be representation research, and it should certainly not be model training.**

Reasoning, in order:

1. **Model training is unauthorised** because no gate has been passed and, more importantly, because the two informative prior experiments are individually uninterpretable (§1.1, §5.5). Training a larger model on an uninstrumented pipeline converts an ambiguous null into an uninterpretable one.
2. **Representation research is premature** because the leading candidate explanation for the prior mutation-model failure — D25, a global embedding being geometrically incapable of expressing a point mutation — is not a research question. It is a one-hour measurement. Run the measurement before commissioning research.
3. **Pipeline correction is the binding constraint.** Six specific instruments are missing and each is cheap: the synthetic planted-signal control, the representation-capability pre-check, the ID-equivalence test, correct censoring handling, cluster-level bootstrapping, and live-regulariser/permutation-sensitivity assertions. Without them no result — positive or negative — carries information. The prior programme's "unexplained positive" and "unexplained failure" are both symptoms of their absence.
4. **Data acquisition follows immediately and is nearly free.** The matched-variant, same-platform surface the programme needs — Duong-Ly's 183 × 76 mutant matrix on the Reaction Biology platform, paired with Anastassiadis's 178 × 300 wild-type matrix on the same platform, plus the mutant constructs already embedded in the KINOMEscan panels — exists in the public domain and takes about two weeks to assemble. A programme that concluded its point-mutant control was "too weak" was not using this surface.
5. **The estimand should change on day one.** Move from strict MMP transformations, which the prior work correctly found too sparse, to the **centred selectivity contrast**, which holds the ligand *identical* rather than merely similar, cancels both dominant confounds structurally, and carries a √2 power advantage over pairwise ΔΔ.

**Ordered plan:** pipeline correction (2–3 weeks) → data acquisition (2 weeks, overlapping) → Stage 0 audit → Stage 0-P positive control → Stage 1. Model training is authorised only at Gate 1, roughly nine weeks out.

**One prediction, recorded so it can be wrong.** I expect Stage 0-P to *pass* on the KINOMEscan and Duong-Ly mutant surfaces once censoring and representation capability are handled — i.e. I expect the prior mutation-model null to reverse. If it does not reverse *after* the planted-signal and capability pre-checks both pass, that is a genuinely informative negative about residue-level conditioning, and it would be the first one this programme has produced.

---

### Verification appendix

**Verified this session from primary or official sources:** Davis 72 × 442 and the explicit statement that blanks correspond to Kd > 10 µM or non-detection at 10 µM; Anastassiadis 178 × 300; Duong-Ly 183 inhibitors × 76 mutant kinases and its CC BY-NC-ND 4.0 licence; PKIS 367 inhibitors × 224 kinases + 24 GPCRs; PKIS2 645 inhibitors and its CC BY 4.0 public release; KINOMEscan scanMAX = 468 kinases including mutant forms; KLIFS 85-position alignment with 0.8 ± 0.1 Å / 2.2 ± 0.2 Å superposition RMSDs; ChEMBL 36 at 17,803 targets and 2.8M compounds; BindingDB at ~3.21M data / 11,414 proteins / >1.41M compounds with mixed CC BY-SA 3.0 and CC BY 3.0 licensing; PDBbind's registration requirement, derivative-distribution prohibition and move to a subscription model; PLINDER at 449,383 systems; SAIR at 5,244,285 structures across 1,048,857 systems under CC BY 4.0; FS-Mol at 5,120 tasks / 233,786 compounds; Landrum & Riniker's curation noise figures; Stanford HIVdb genotype–phenotype coverage.

**Requires Stage-0 re-verification from supplementary files:** Metz ≈3,858 × 172; the exact Davis censored fraction (the ~9,166 detected-interaction figure is from secondary literature); Klaeger 243 × ~253; Karaman 38 × 317 with 27 mutants; PKIS 21 disease-relevant mutants; KIBA 2,111 × 229; the precise mutant and phospho-variant composition of Davis's 442 assays. Several of these came from a single review table whose row/citation alignment was ambiguous; they are reported here as approximate and must be read off the primary supplementary files before any power calculation depends on them.
