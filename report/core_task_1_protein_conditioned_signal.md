# Core Task 1 — Does a Transferable, Affinity-Relevant, Protein-Conditioned Interaction Signal Exist and Can It Be Learned for Cold-Target DTA?

**Role:** independent senior researcher (computational biology / medicinal chemistry / ML / scientific evaluation)
**Status:** pre-registration-grade research program. Hypothesis generation is broad; evidence standards are deliberately severe.
**Date:** 17 August 2026

---

## 0. Executive scientific judgment

**Short version.** The signal almost certainly *exists* in nature. Whether it is *recoverable at scale from public data* and *transferable to genuinely cold protein clusters* is unproven, and nearly every published claim to the contrary is confounded. The bottleneck is not architecture. It is the **choice of estimand**. I recommend reformulating Core Task 1 so that the prediction target is an interaction contrast in which the ligand-only baseline and the target-ID baseline are *structurally forced to zero*, rather than merely "outperformed."

**Six judgments, with confidence:**

1. **A protein-conditioned interaction term exists biologically (confidence: very high, >0.95).** Selectivity is the existence proof. A single kinase inhibitor spans >4 log units of Kd across the kinome; that dispersion is by definition non-additive in a ligand main effect + protein main effect model. The question was never "does I(ligand, protein) exist" but "is it estimable and transferable."

2. **Most reported cold-target DTA gains are not evidence for it (confidence: high, ~0.85).** Volkov et al. showed that explicit protein–ligand interaction descriptors gave no advantage over ligand-only or protein-only descriptors, and that nearest-neighbour lookups already perform well — i.e. memorization dominates ([J. Med. Chem. 2022](https://pubs.acs.org/doi/abs/10.1021/acs.jmedchem.2c00487)). Graber et al. reach the same conclusion from the data-bias side ([Nat. Mach. Intell. 2025](https://www.nature.com/articles/s42256-025-01124-5)). Leave-cluster-out validation has been known to be necessary since Kramer & Gedeck (2010).

3. **The dominant confound is that a global protein embedding is a near-perfect target identifier.** With ~10³–10⁴ targets and a per-target affinity mean that varies by several log units, a model that learns only `a(protein)` scores extremely well on warm targets and degrades gracefully on cold ones via family similarity — producing exactly the phenomenology usually reported as "protein information helps."

4. **The second confound is assay heterogeneity, and it is large enough to swamp the effect being sought.** Landrum & Riniker quantified it: with minimal curation, ~65% of same-compound/same-target IC50 pairs from different literature assays differ by >0.3 log units and 27% by >1 log unit; maximal metadata curation improves this to 48% / 13% with Kendall's τ ≈ 0.71, at the cost of dataset size ([JCIM 2024](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049)). Any claimed interaction effect smaller than this noise floor is not interpretable.

5. **The right estimand is a double difference (confidence: high).** For two ligands (l₁, l₂) and two proteins (p, q):
   `ΔΔ = [f(l₁,p) − f(l₂,p)] − [f(l₁,q) − f(l₂,q)]`
   This cancels **every** ligand main effect, **every** protein main effect, and **every additive assay offset that is constant within an assay**. A ligand-only model predicts ΔΔ = 0 identically. A target-ID model predicts ΔΔ = 0 identically. Therefore *any* predictive skill above zero is protein-conditioned by construction. This is not a metric preference; it is a change in what is being estimated. Empirical support for the noise-cancellation property of paired data comes from Nelen et al. and is reviewed in [Fralish & Reker, *Front. Drug Discov.* 2026](https://www.frontiersin.org/journals/drug-discovery/articles/10.3389/fddsv.2026.1859068/full).

6. **My prior on the outcome (calibrated, stated in advance so it can be wrong):**
   - Interaction signal reproducible *within* a superfamily across independent assay platforms: **p ≈ 0.75**
   - Interaction signal transferable to held-out protein clusters *within* a superfamily (e.g. cold kinase groups): **p ≈ 0.5**
   - Transferable across superfamilies (kinase → GPCR → protease) at useful effect size from public data alone: **p ≈ 0.15**
   - Therefore the realistic deliverable is **within-superfamily zero-shot + cross-superfamily few-shot**, and the program should be scoped that way from day one rather than discovering it at month nine.

**The single most valuable thing in this document** is §7 (T0-KILL): a one-day, zero-GPU experiment that can terminate the entire program by showing that the interaction residual is not even reproducible *between independent experiments*, in which case no model can learn it and Core Task 1 is unresolvable on public affinity data.

---

## 1. Formal decomposition: the three questions, separated

Let `f(l,p)` be a comparable affinity measurement (fixed label system, e.g. pKd). Write the cross-classified decomposition:

```
f(l,p) = μ + a(p) + b(l) + I(l,p) + e(l,p,assay,doc)
```

- `a(p)`: protein main effect — target-level "bindability", assay platform offset, dynamic range, publication selection. **This is what a target-ID / global protein embedding recovers.**
- `b(l)`: ligand main effect — promiscuity, lipophilicity, molecular size, reactivity. **This is what a ligand-only model recovers.**
- `I(l,p)`: the interaction term. **This is Core Task 1.**
- `e`: measurement + assay-comparability noise, quantified above.

The three questions the brief insists on separating map cleanly onto this:

| Question | Formal statement | How it is answered | Requires a neural model? |
|---|---|---|---|
| **Q1. Is transferable protein-conditioned information present?** | Is `Var[I]` large relative to `Var[e]`, and is `I` reproducible across independent assay platforms? | Variance-component estimation + cross-panel reproducibility of the interaction residual | **No.** Mixed-effects models, no learning. |
| **Q2. Can a representation recover it?** | Does there exist `g(z(l), x(p))` with `corr(g, I) > 0` on **held-out protein clusters**? | Kernel/GP two-way model on the interaction residual with cold-protein splits | **No.** Ridge/kernel methods suffice and are preferable. |
| **Q3. Does it improve cold-target DTA?** | Does `μ̂ + â(p) + b̂(l) + ĝ` beat `μ̂ + â(p) + b̂(l)` on cold-target RMSE/CI/Spearman? | Full DTA evaluation with all controls | Only after Q1 and Q2 pass. |

**Conflation is the field's core methodological error.** Papers routinely answer Q3 affirmatively while having answered neither Q1 nor Q2, because Q3 is satisfiable by `â(p)` alone.

**Structural consequence (important).** Because ligand-only and target-ID models both predict ΔΔ ≡ 0, the standard control battery changes character. The controls stop being "baselines to beat" and become **falsifiers of the pipeline**: if shuffled-protein or random-embedding branches achieve non-zero ΔΔ skill, that is proof of leakage in the split, not of a weaker-but-real effect.

**Noise ceiling, explicitly.** If within-assay per-measurement SD is σ, then ΔΔ (four measurements) has noise SD 2σ. If the true interaction contrast has SD τ, the maximum achievable correlation is

```
r_max = τ / sqrt(τ² + 4σ²)
```

For σ = 0.3 log units: τ = 0.8 → r_max ≈ 0.80; τ = 0.4 → r_max ≈ 0.55; τ = 0.2 → r_max ≈ 0.32. **All reported performance must be expressed as a fraction of r_max**, not as raw correlation, or the results are uninterpretable across datasets.

---

## 2. Literature map (primary sources, direct links)

### 2.1 Negative results and shortcut evidence — read these first

| Work | Why it matters here | Link |
|---|---|---|
| Volkov et al., *J. Med. Chem.* 2022 — "On the Frustration to Predict Binding Affinities…" | Explicit non-covalent interaction descriptors gave no advantage over ligand-only/protein-only; nearest-neighbour lookup already strong. The canonical negative control for this program. | [doi:10.1021/acs.jmedchem.2c00487](https://pubs.acs.org/doi/abs/10.1021/acs.jmedchem.2c00487) |
| Graber et al., *Nat. Mach. Intell.* 2025 — "Resolving data bias improves generalization in binding affinity prediction" | Characterizes leakage/bias in protein–ligand datasets; clean splits change conclusions. | [doi:10.1038/s42256-025-01124-5](https://www.nature.com/articles/s42256-025-01124-5) |
| Landrum & Riniker, *JCIM* 2024 — "Combining IC50 or Ki Values from Different Sources Is a Source of Significant Noise" | The quantitative assay-noise floor. Sets the detection limit for the entire program. | [doi:10.1021/acs.jcim.4c00049](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049) |
| van Tilborg, Alenicheva & Grisoni, *JCIM* 2022 — MoleculeACE | Activity-cliff stratified evaluation; descriptor models beat deep models on cliffs. Supplies the cliff strata required by the brief. | [doi:10.1021/acs.jcim.2c01073](https://pubs.acs.org/jcisd8/article/62/23/5938/852680/Exposing-the-Limitations-of-Molecular-Machine) · [code](https://github.com/molML/MoleculeACE) |
| Mastropietro, Pasculli & Bajorath, *Nat. Mach. Intell.* 2023 — "Learning characteristics of GNNs predicting protein–ligand affinities" | What GNNs actually key on. | doi:10.1038/s42256-023-00756-9 |
| Chatterjee et al., *Nat. Commun.* 2023 — AI-Bind | Degree/annotation shortcuts in DTI networks; generalization to unseen nodes. | doi:10.1038/s41467-023-37572-z |
| Kramer & Gedeck, *JCIM* 2010 — leave-cluster-out CV | The original statement that random splits are invalid for scoring functions on diverse protein sets. | doi:10.1021/ci100264e |

### 2.2 Chemogenomics / proteochemometrics — the intellectual ancestor of Core Task 1

- van Westen et al., *MedChemComm* 2011 — PCM as a tool for selective compound design and **extrapolation to novel targets**: [doi:10.1039/C0MD00165A](http://xlink.rsc.org/?DOI=c0md00165a)
- van Westen et al., *J. Cheminform.* 2013 — benchmarking 13 amino-acid descriptor sets, parts 1 & 2: [part 2](https://link.springer.com/article/10.1186/1758-2946-5-42)
- Cortés-Ciriano et al. — PCM in a Bayesian framework: [PMC4083135](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4083135/)
- **Relevance:** PCM already established that concatenated ligand+protein descriptors can extrapolate in target space. It also established the failure mode: with concatenation and a flexible learner, the protein block collapses to an identifier. Modern deep DTA largely rediscovered PCM without inheriting its controls.

### 2.3 Interaction-modelling architectures

- Bai et al., *Nat. Mach. Intell.* 2023 — **DrugBAN**, bilinear attention + conditional domain adversarial adaptation for cross-domain DTI: [doi:10.1038/s42256-022-00605-1](https://www.nature.com/articles/s42256-022-00605-1) · [arXiv](https://arxiv.org/abs/2208.02194) · [code](https://github.com/peizhenbai/DrugBAN)
- Reusability report on DrugBAN (Xu et al., *Nat. Mach. Intell.* 2024) — independent reproduction of the cross-domain claims: [doi:10.1038/s42256-024-00822-w](https://www.nature.com/articles/s42256-024-00822-w). **Read the reusability report before adopting the method.**
- Singh, Sledzieski et al., *PNAS* 2023 — **ConPLex**: protein-anchored contrastive co-embedding on PLM features, with prospective kinase validation: [doi:10.1073/pnas.2220778120](https://www.pnas.org/doi/10.1073/pnas.2220778120) · [code](https://github.com/samsledje/ConPLex)
- Chen et al. — TransformerCPI, including **label-reversal experiments** (a genuine shortcut control, rare in this literature).

### 2.4 Protein representation at residue level

- Rao et al., ICLR 2021 — transformer PLM attention maps are unsupervised contact learners: [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.12.15.422761v1) · [OpenReview](https://openreview.net/forum?id=fylclEqgvgd) · [ESM code](https://github.com/facebookresearch/esm)
- Lin et al., *Science* 2023 — ESM-2 / ESMFold (per-residue embeddings at scale).
- Kooistra et al., *NAR* 2016 — **KLIFS**: a consistent 85-residue kinase pocket alignment with sub-Å superposition RMSD for aligned residues; the single most useful protein representation resource for this program: [doi:10.1093/nar/gkv1082](https://academic.oup.com/nar/article/44/D1/D365/2502606) · [klifs.net](https://klifs.net/) · [original JMC 2014](https://pubs.acs.org/doi/10.1021/jm400378w)
- GPCRdb (generic residue numbering) — the GPCR analogue, enabling the same aligned-pocket construction in an independent superfamily: [gpcrdb.org](https://gpcrdb.org)

### 2.5 Few-shot / meta-learning / neural processes

- Stanley et al., NeurIPS D&B 2021 — **FS-Mol**: 5,120 protein-target tasks, 233,786 compounds, 4,938/40/157 train/val/test task split: [paper PDF](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/file/8d3bba7425e7c98c50f52ca1b52d3735-Paper-round2.pdf)
- Chen et al., ICLR 2023 — **ADKF-IFT** (adaptive deep kernel GPs): [arXiv:2205.02708](https://arxiv.org/pdf/2205.02708)
- Garnelo et al., 2018 — Conditional Neural Processes (the amortized-conditioning formalism used in Direction D4)
- Schimunek et al., ICLR 2023 — MHNfs (context molecules for few-shot)
- **Critical caveat:** FS-Mol tasks are *protein-defined*, but the standard protocol conditions on a support set, **not** on protein features. This makes FS-Mol an almost perfect ablation surface: support-set-only vs support-set + protein features isolates exactly the incremental value of protein information.

### 2.6 Pairwise / matched-molecular-pair learning (the machinery for the recommended estimand)

- Fralish & Reker, *Front. Drug Discov.* 2026 — review of pairwise learning; notes that paired potency differences partially normalize systematic inter-assay bias (citing Nelen et al. 2025): [doi:10.3389/fddsv.2026.1859068](https://www.frontiersin.org/journals/drug-discovery/articles/10.3389/fddsv.2026.1859068/full)
- Tynes et al., *JCIM* 2021 — **PADRE**, pairwise difference regression with UQ: [doi:10.1021/acs.jcim.1c00670](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00670)
- Tyrchan & Evertsson — MMP analysis methods and applications review, *J. Med. Chem.*: [doi:10.1021/acs.jmedchem.2c01787](https://pubs.acs.org/doi/10.1021/acs.jmedchem.2c01787)
- Hussain & Rea, *JCIM* 2010 — MMP fragmentation algorithm; `mmpdb` (Dalke et al.) is the practical open implementation.

### 2.7 Structure and co-folding (use with discipline)

- Passaro et al., 2025 — **Boltz-2**, joint structure + affinity co-folding: [bioRxiv 2025.06.14.659707](https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1.full.pdf) · [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12262699/)
- **Counter-evidence to weigh:** an independent 2026 evaluation reports that Boltz-2 affinity outputs appear to derive largely from features **independent of the final ligand pose** ([arXiv:2603.05532](https://arxiv.org/html/2603.05532v1)). If true, co-folding affinity heads are *not* evidence of Cartesian recognition and must be treated as learned chemogenomic priors, not physics.
- Durairaj et al., 2024 — **PLINDER**: 449,383 PLI systems, similarity metrics at protein/pocket/interaction/ligand levels, leakage-minimizing splits: [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.17.603955v3) · [code](https://github.com/plinder-org/plinder) · [plinder.sh](https://www.plinder.sh/)

### 2.8 Adjacent fields — how they prove "modality X contributes conditional information, not an ID"

This is where the brief's request for cross-disciplinary transfer pays off. Four transferable protocols:

1. **VQA / multimodal NLP — the unimodal-prior problem.** Agrawal et al., *VQA-CP* (CVPR 2018) reconstructed test splits so that answer priors are *inverted* relative to train. The transfer: build a **prior-inverted cold-target split** where the family-level average selectivity direction in test is the opposite of train. A model riding on family priors collapses; a model using residue-level chemistry does not.
2. **Multimodal collapse / greedy learners.** Wang, Tran & Feiszli, "What Makes Training Multi-Modal Classification Networks Hard?" (CVPR 2020); Peng et al., OGM-GE (CVPR 2022). Transfer: measure **per-modality gradient contribution** and per-modality generalization gap; if the protein branch's gradient norm collapses early in training, modality collapse — not absence of signal — explains a null result. This distinction is essential to avoid falsely falsifying Core Task 1.
3. **Causal / invariant learning.** Arjovsky et al., IRM; Krueger et al., REx; Sagawa et al., GroupDRO. Transfer: treat **assay/document as the environment variable**. A protein-conditioned effect that is not invariant across documents is an assay artifact.
4. **System identification / operator learning.** The classic identifiability question — can a system's response operator be identified from input–output pairs without persistent excitation? Transfer: the ligand×protein matrix must have sufficient **excitation** (shared ligands across protein clusters) or `I` is *unidentifiable in principle*, regardless of model. Stage 0 must verify excitation before anything is trained. This is the most underused idea in the DTA literature.

### 2.9 Shortcut-learning meta-literature

- Lapuschkin et al., *Nat. Commun.* 2019 — "Unmasking Clever Hans predictors"
- Geirhos et al., *Nat. Mach. Intell.* 2020 — "Shortcut learning in deep neural networks"
- Wallach & Heifets, *JCIM* 2018 — AVE bias; benchmarks that reward memorization
- Sieg, Flachsenberg & Rarey, *JCIM* 2019 — "In Need of Bias Control…" (structure-based VS)
- Chen et al., *PLoS ONE* 2019 — hidden bias in DUD-E
- Tran-Nguyen et al., *JCIM* 2020 — LIT-PCBA (unbiased experimental VS sets)
- A 2025 preprint reports that *chemist-style* signals (author/lab-specific structural idiom) confound activity prediction on public benchmarks — a document-level confound distinct from assay noise: [arXiv:2512.20924](https://arxiv.org/pdf/2512.20924). Treat as motivating a document control, not as established fact.

---

## 3. Public data source comparison

Legend for **Cold-protein test capability**: ✅ suitable as primary evidence · ⚠️ usable as secondary/replication only · ❌ cannot support the test.

| Resource | Label semantics | Assay comparability | Protein × ligand coverage | Structural coverage | Access | Principal leakage risk | Cold-protein test |
|---|---|---|---|---|---|---|---|
| **Davis et al. 2011** (KINOMEscan) | Kd, competition binding, single platform | **Excellent** (one platform, one lab) | ≈72 inhibitors × ≈442 kinases, near-dense | Indirect (via KLIFS) | Supplementary to *Nat. Biotechnol.* 29:1046 | Same ligands recur across all proteins → **ligand-identity leakage is guaranteed unless explicitly excluded** | ✅ **primary** |
| **Anastassiadis et al. 2011** | % inhibition @ single conc. — **activity, not affinity** | Excellent internally | ≈178 × ≈300 | Indirect | *Nat. Biotechnol.* 29:1039 supplementary | Same as above; also censoring at 0/100% | ⚠️ ordering only, never merge with Kd |
| **Metz et al. 2011** | pKi, kinase panel | Good internally | ≈3,800 compounds × ≈172 kinases | Indirect | *Nat. Chem. Biol.* 7:200 | Series clustering | ✅ **primary** (best ligand diversity of the panels) |
| **Klaeger et al. 2017** (Kinobeads) | apparent Kd from **lysate competition chemoproteomics** — different equilibrium semantics | Excellent internally, **not** interchangeable with biochemical Kd | 243 clinical drugs × several hundred proteins | Indirect | [*Science* 358:eaan4368](https://www.science.org/doi/10.1126/science.aan4368); PRIDE PXD005336 | Clinical drugs are heavily represented in ChEMBL → cross-source contamination | ✅ **independent replication cohort** |
| **PKIS / PKIS2** | % inhibition panels | Good internally | ~hundreds × hundreds | Indirect | Open (GSK/SGC) | Congeneric series | ⚠️ replication |
| **ChEMBL** (v35+) | IC50/Ki/Kd/EC50, **heterogeneous** | **Poor across assays** unless maximally curated | ~2.4M compounds, ~15k targets, extremely sparse | Via mapping | CC BY-SA 3.0, [ebi.ac.uk/chembl](https://www.ebi.ac.uk/chembl/) | Document/series/chemist-style; missing-not-at-random negatives | ⚠️ only with within-document pairing |
| **BindingDB** | Ki/Kd/IC50 aggregated (ChEMBL, PubChem, patents) | Poor across sources | Large | Links to PDB | [bindingdb.org](https://www.bindingdb.org) — verify current terms | Re-aggregation → duplicate records across "independent" sources | ⚠️ |
| **Papyrus** | Standardized/normalized ChEMBL + ExCAPE + curated sets; ~60M points, ~1.24M flagged high-quality exact values | Improved but still cross-assay | Very large | No | [doi:10.1186/s13321-022-00672-x](https://link.springer.com/article/10.1186/s13321-022-00672-x); data DOI 10.4121/16896406 | Inherits ChEMBL confounds; standardization ≠ harmonization | ⚠️ good for curation infrastructure |
| **KIBA** | **Merged** Ki/Kd/IC50 into a single synthetic score | **Fails the brief's harmonization requirement** | 2,111 × 229 | No | Widely redistributed | Merged label semantics; heavily overfit benchmark | ❌ **do not use as evidence** |
| **Drug Target Commons** | Crowd-curated with assay annotation | Variable, but annotated | Kinase-heavy | No | [drugtargetcommons.fimm.fi](https://drugtargetcommons.fimm.fi) | Curation heterogeneity | ⚠️ |
| **PDBbind** | Kd/Ki/IC50 attached to crystal complexes | Poor (aggregated literature) | ~20k complexes, ~5k proteins | ✅ full | Registration; non-commercial terms — **verify** | Documented train/test leakage; small | ❌ for Q1; ⚠️ for structural priors |
| **Binding MOAD / BioLiP2** | Curated complexes + some affinity | Mixed | Large structural | ✅ full | Open (verify) | Functional-annotation focus; no ML splits | ⚠️ structural priors only |
| **PLINDER** | Structural systems + annotations (+ affinity annotations in progress) | n/a | 449,383 systems | ✅ full + apo + predicted | [github.com/plinder-org/plinder](https://github.com/plinder-org/plinder) | **This is the leakage-control resource, not the affinity resource** | ✅ **for split construction** |
| **KLIFS** | Not affinity — 85-residue aligned kinase pocket + interaction fingerprints | n/a | All human/mouse kinase domain structures | ✅ | [klifs.net](https://klifs.net/) (free academic; verify) | n/a | ✅ **protein representation backbone** |
| **GPCRdb** | Generic residue numbering, structures, mutations | n/a | GPCR superfamily | ✅ | [gpcrdb.org](https://gpcrdb.org) | n/a | ✅ **independent superfamily replication** |
| **FS-Mol** | Binary activity (regression variant available), ChEMBL-derived | Task-internal | 5,120 tasks / 233,786 compounds | No | [NeurIPS D&B 2021](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/file/8d3bba7425e7c98c50f52ca1b52d3735-Paper-round2.pdf) | Tasks defined by assay → task = assay confound | ✅ **for the support-set ablation** |
| **MoleculeACE** | pKi/pIC50, 30 single targets, cliff-annotated | Per-target curated | 30 targets | No | [github.com/molML/MoleculeACE](https://github.com/molML/MoleculeACE) | Single-target → no protein contrast | ❌ for Q1; ✅ for cliff strata definitions |
| **Stanford HIVdb** | Fold-change resistance, fixed drugs × protein variants | Good within phenotype platform | ~10 drugs × thousands of variants | ✅ (protease/RT) | [hivdb.stanford.edu](https://hivdb.stanford.edu) | Variant correlation structure | ✅ **highest-purity existence test** (§ D5) |
| **Ortholog panels** (e.g. human/rat adenosine receptors, van Westen) | Ki, same lab | Excellent | Small | Some | ChEMBL + primary papers | Small n | ✅ ideal similarity-matched control |
| **BELKA / DEL** | DEL enrichment, not affinity | Internally consistent | Billions of molecules × **3 proteins** | No | Kaggle/Leash | 3 proteins | ❌ cannot test protein-cold |

**Harmonization ruling.** Do not merge label systems. The one defensible harmonization is the one the brief implicitly asks for: **differences within an assay are comparable even when absolute values across assays are not.** This is exactly why the ΔΔ estimand is the correct scientific object here, and it is supported empirically (Nelen et al. via Fralish & Reker 2026; Landrum & Riniker 2024). Every direction below therefore trains **per label system**, with cross-system agreement used only as external replication.

---

## 4. Ranked solution portfolio

Ranking criterion: `P(produces confirmatory evidence for Core Task 1 within 6 months) × (value of that evidence) ÷ cost`.

---

### **D1 — Protein-conditioned transformation transfer (matched-pair ΔΔ learning)** ★ Rank 1 — *selected: highest-probability direction*

1. **Biological hypothesis.** A chemical transformation τ (e.g. Cl→CF₃ at a specific attachment vector) produces a potency change whose *sign and magnitude depend on the local protein environment it contacts* (gatekeeper size, hinge donor/acceptor pattern, back-pocket polarity). SAR rules are pocket-conditional, not universal.
2. **Why the protein must alter the ligand effect.** A transformation that adds bulk is tolerated by a small-gatekeeper pocket and rejected by a large one; a transformation adding an H-bond donor pays off only where a complementary acceptor exists. These are first-order medicinal chemistry facts.
3. **Minimum information / structure.** Same transformation τ observed on ≥2 protein clusters, each observation being a within-assay pair. Minimum viable: ~5,000 (τ, protein) observations spanning ≥30 protein clusters at <50% pocket identity, with ≥3 clusters held out.
4. **Representation and objective.** Ligand side: MMP transformation encoding (fragment pair + environment radius) from `mmpdb`. Protein side: aligned pocket residue vector (KLIFS/GPCRdb) or pocket-restricted PLM embeddings. Objective: regress the **interaction residual** `Δ(τ|p) − E_p[Δ(τ|·)]` — i.e. subtract the transformation's global average effect first, so the model *cannot* score by learning universal SAR. Loss: Huber on residual + antisymmetry regularizer.
5. **Why it is not a target ID.** The prediction target is a residual whose expectation over proteins is zero. A target ID carries no information about it. Formally: `E_p[target] = 0` ⟹ any constant-per-protein feature has zero predictive power.
6. **Failure modes / shortcuts.** (a) Transformation identity leaking through to the residual if the transformation appears on only one protein family — mitigate by requiring τ to appear on ≥3 clusters; (b) series/document leakage — same paper reports both proteins; (c) apparent effect driven entirely by a handful of promiscuous "magic methyl" transformations.
7. **Exact negative controls.** Ligand-only (predicts 0 by construction — verify empirically as a pipeline check); shuffled protein; similarity-matched wrong protein (paralog matched on pocket identity ±5%); capacity-matched random embedding in two flavours — iid Gaussian and *family-structure-preserving scramble*; pocket residue permutation (order shuffled, composition preserved); composition scramble (positions preserved, identities randomized); document-held-out; τ-held-out; ligand-identity exclusion.
8. **Low-cost falsification.** §7's T1 experiment. Kernel ridge on ~10⁴ residuals, CPU only, <1 day.
9. **Compute / complexity.** Stage 1: CPU-hours. Stage 2 (cross-attention version): 1 GPU-week. Low.
10. **Zero-/few-shot support.** Excellent. Zero-shot: rank analogues on a novel target. Few-shot: k anchor compounds identify `a(p)`; the transformation model supplies the ordering. This maps directly onto real lead-optimization use.

---

### **D2 — Alignment-anchored residue vocabulary + fragment cross-attention** ★ Rank 2 — *selected: representation-centric direction*

1. **Hypothesis.** Interaction information lives at specific *aligned positions*, not in a pooled global embedding. A model given a fixed-length aligned pocket (KLIFS 85 positions; GPCRdb generic numbering) can learn position-specific chemistry that transfers across proteins sharing the alignment.
2. **Why the protein alters the effect.** Position 45 (gatekeeper) determines back-pocket access; the hinge triplet determines donor/acceptor complementarity. These are positionally identifiable across the whole superfamily.
3. **Minimum information.** An alignment covering the target family + ≥50 proteins with ≥30 shared ligands each. KLIFS provides the kinase case with sub-Å superposition for aligned residues.
4. **Representation / objective.** Protein: `[85 × d]` matrix (one-hot AA + Z-scale physicochemical + optional pocket-restricted ESM-2 embeddings). Ligand: fragment/pharmacophore tokens. Objective: cross-attention → low-rank bilinear interaction head, trained on the interaction residual, with a position-sparsity penalty (L1 over attention mass across the 85 positions).
5. **Why not a target ID.** Perform **single-position substitution probes**: mutate one aligned residue *in the input only* and require the prediction change to correlate with real mutational data (gatekeeper mutants, resistance panels). An ID-based model produces no coherent response to single-position edits.
6. **Failure modes.** Alignment quality degrades outside the superfamily; the model may learn family classification from the pocket string; attention may be diffuse and uninterpretable.
7. **Negative controls.** All of D1's, plus **position-shuffled alignment** (columns permuted consistently across all proteins — preserves everything except positional semantics; this is the sharpest single control available for this direction).
8. **Low-cost falsification.** Before any deep model: fit an L1-regularized linear model on `(position, amino-acid class) × ligand-cluster` interaction terms. If no sparse set of positions carries signal, cross-attention will not find one.
9. **Compute.** 1–2 GPU-weeks. Medium.
10. **Zero-/few-shot.** Zero-shot within the aligned superfamily is natural. Cross-superfamily requires a learned alignment-free variant — a known hard problem.

---

### **D3 — Weakly supervised latent binding-region discovery via multiple-instance learning** ★ Rank 3 — *selected: high-risk / high-reward direction*

1. **Hypothesis.** Even without legal complex coordinates, affinity supervision across many ligands can identify *which* protein region governs the response — recoverable as a latent variable.
2. **Why the protein alters the effect.** Only pocket residues contact the ligand; the rest of the sequence is nuisance. If a model must select a region, that selection is the interaction hypothesis.
3. **Minimum information.** Sequence + candidate pocket proposals (fpocket / P2Rank / AF2-derived), affinity labels only. No complex structures needed.
4. **Representation / objective.** MIL: bag = candidate pockets; instance score = pocket-conditioned affinity prediction; aggregation = attention or noisy-OR. Add a mutual-information penalty between the pocket-selection distribution and target identity.
5. **Why not a target ID.** The selection variable is *within*-protein; a target ID cannot pick a pocket.
6. **Failure modes.** Degenerate selection (always pocket #1); selection correlated with pocket size/hydrophobicity rather than the true site; **attention treated as explanation without validation** — the exact error the brief forbids.
7. **Negative controls.** Random pocket assignment; largest-pocket heuristic; and critically an **external, pre-registered validation**: pocket-recovery AUC against held-out annotated sites (KLIFS/BioLiP2) that were *never* used in training. This converts "the attention looks biological" into a falsifiable prediction with a stated threshold (pre-register AUC ≥ 0.75).
8. **Low-cost falsification.** On a small set of proteins with known sites, check whether an oracle that is *given* the true pocket outperforms whole-sequence input. If knowing the true pocket does not help, MIL discovery of it cannot help either. **This ordering — oracle first, discovery second — should be mandatory.**
9. **Compute.** 2–4 GPU-weeks. Medium-high.
10. **Zero-/few-shot.** Strong if it works; the discovered region generalizes to novel folds where alignments do not exist. This is why it is the high-reward option.

---

### **D4 — Conditional neural process / task-conditioned operator over targets** ★ Rank 4

1. **Hypothesis.** A protein is a *task*. Protein features should reduce the number of support examples needed — measurable as a shift in the few-shot learning curve.
2. **Why the protein alters the effect.** If protein features carry conditional information, they act as an informative prior over the task function.
3. **Minimum information.** Many tasks with variable support sizes: FS-Mol (5,120 tasks) is purpose-built.
4. **Representation / objective.** CNP/ANP: encoder over support (ligand, affinity) pairs → task representation `r`; concatenate protein features `x(p)`; decode query affinity. Train with the standard CNP likelihood.
5. **Why not a target ID.** The decisive ablation is **`r` alone vs `r` + `x(p)`, as a function of support size k**. A target ID is redundant with `r` once k ≳ 5, so an ID-like protein feature yields curves that converge; genuinely conditional information yields a *persistent* gap and a leftward shift, and non-zero skill at **k = 0**.
6. **Failure modes.** Protein branch ignored (modality collapse); tasks defined by assay, so "task" conflates protein and assay.
7. **Negative controls.** Shuffled protein features across tasks; random capacity-matched features; **assay-matched-but-different-protein** tasks; k=0 evaluation with all controls.
8. **Low-cost falsification.** Run the k = 0 vs k = 16 curve with a simple deep-kernel GP. If the k=0 gap between correct and shuffled protein features is within CI of zero, the representation carries no conditional information.
9. **Compute.** 1 GPU-week (FS-Mol infrastructure exists). Low-medium.
10. **Zero-/few-shot.** This *is* the zero-/few-shot formulation. Highest direct product relevance.

---

### **D5 — Protein-variant conditioning: orthologs, point mutants, resistance panels** ★ Rank 5 — *selected: data-centric direction*

1. **Hypothesis.** If the *same* ligand set is measured against proteins differing by a handful of residues, in the same assay, any systematic affinity change is unambiguously protein-conditioned.
2. **Why the protein alters the effect.** Gatekeeper mutations (ABL T315I, EGFR T790M/C797S), HIV protease resistance positions, ortholog substitutions in the pocket.
3. **Minimum information.** Fixed ligand set × protein variants, single assay. Stanford HIVdb; kinase resistance series; human/rat ortholog panels.
4. **Representation / objective.** Predict `Δf(l, wild-type → mutant)` from (ligand, mutated position, substitution). Objective: signed regression + sign classification.
5. **Why not a target ID.** Variants of the same protein share ~99% sequence; a global embedding barely distinguishes them. **Distinguishing them is the task.**
6. **Failure modes.** Small n; strong correlation among resistance variants; fold-change is a ratio with its own error model; publication bias toward large effects.
7. **Negative controls.** Randomize which position mutated; substitute a chemically matched but structurally distal position (surface residue at matched BLOSUM distance) — this is the ideal "similarity-matched wrong protein" control because it is matched at the *residue* level, not the protein level.
8. **Low-cost falsification.** Two-way ANOVA on the drug × variant matrix: is the interaction term significant after protein and drug main effects? One afternoon.
9. **Compute.** Trivial. CPU.
10. **Zero-/few-shot.** Limited direct transfer to novel families, but it is the **cleanest identifiability proof available** and therefore the correct sanity gate: a method that fails here should not be trusted on cold families.

---

### **D6 — Explicit additive/interaction decomposition with orthogonality + conditional-information penalty** ★ Rank 6 — *selected: training-objective innovation*

1. **Hypothesis.** Making the interaction term a *named, isolated model component* prevents it from absorbing main effects and makes its contribution directly measurable.
2. **Why the protein alters the effect.** Same as D1/D2; this is about estimation, not biology.
3. **Minimum information.** Any ligand×protein matrix with adequate excitation.
4. **Representation / objective.**
   `f̂ = μ + a_θ(p) + b_φ(l) + ⟨U x(p), V z(l)⟩`
   with (i) an **orthogonality penalty** forcing the bilinear term's output to be uncorrelated with both main effects; (ii) a **conditional-information penalty** — maximize `I(y ; x(p) | l)` while penalizing `I(x(p) ; target_id)` via a CLUB/InfoNCE-style bound; (iii) **gradient-conflict control** (PCGrad or OGM-GE-style modulation) so the protein branch is not starved by the faster-learning ligand branch.
5. **Why not a target ID.** `a_θ(p)` is *given* the target-ID job explicitly. Anything the bilinear term contributes is, by construction, above and beyond it. The reportable quantity becomes **variance explained by the interaction component on cold targets** — a direct answer to Core Task 1.
6. **Failure modes.** Orthogonality penalty too weak (leakage into the bilinear term) or too strong (kills real signal); adversarial CMI estimators are unstable; rank of `U,V` acts as capacity and must be matched in controls.
7. **Negative controls.** Rank sweep with all controls at each rank; **capacity-matched random embedding at identical rank** (essential — otherwise the comparison is confounded by capacity); gradient-norm logging per branch to detect collapse.
8. **Low-cost falsification.** Fit the decomposition with `x(p)` = real vs random at several ranks. If the interaction variance on cold targets is identical, the objective is not isolating anything.
9. **Compute.** Low-medium; 3–5 GPU-days.
10. **Zero-/few-shot.** Very good: `a(p)` is exactly the term few-shot anchors calibrate, and the interaction term is exactly the transferable part.

---

### **D7 — Assay-environment invariance / domain generalization**  ★ Rank 7

1. **Hypothesis.** Genuine interaction effects are invariant across assays and documents; artifacts are not.
2. **Why the protein alters the effect.** Not a new mechanism — a filter on candidate mechanisms.
3. **Minimum information.** Environment labels (ChEMBL assay_id, doc_id) on every datapoint. Non-negotiable metadata requirement.
4. **Representation / objective.** IRM / V-REx / GroupDRO with environment = document (or assay platform); report worst-environment performance.
5. **Why not a target ID.** Target ID correlates with document (many targets appear in few documents), so the invariance penalty actively suppresses ID-like solutions.
6. **Failure modes.** IRM is known to underperform ERM under many realistic shifts; too few environments per target; environment ≈ target, making the penalty degenerate.
7. **Negative controls.** ERM baseline; random environment assignment; environment-shuffled-within-target.
8. **Low-cost falsification.** Measure environment/target confounding first (mutual information between doc_id and target_id). If it is near-maximal, the method is inapplicable and this must be known before implementation.
9. **Compute.** Low.
10. **Zero-/few-shot.** Indirect; improves robustness rather than adding conditional information.

---

### **D8 — Protein-conditioned metric learning with selectivity-mined hard negatives**  ★ Rank 8

1. **Hypothesis.** A metric space in which "close" means "binds this pocket" carries transferable interaction structure (ConPLex line of work).
2. **Why the protein alters the effect.** The anchor is the protein; ligand embeddings are compared conditionally.
3. **Minimum information.** Positives/negatives per protein; ideally *selectivity pairs* — same ligand, active on A and inactive on B.
4. **Representation / objective.** Protein-anchored triplet/InfoNCE loss (per ConPLex, [PNAS 2023](https://www.pnas.org/doi/10.1073/pnas.2220778120)) but with negatives mined for **measured selectivity**, not random or chemically-similar decoys, and with an affinity-signed margin rather than binary contrast.
5. **Why not a target ID.** Hard negatives are matched on ligand identity, so the only way to satisfy the loss is to encode pocket differences.
6. **Failure modes.** Decoy bias (the classic virtual-screening failure — see AVE bias literature); binary framing discards magnitude; false negatives from missing-not-at-random data.
7. **Negative controls.** Random negatives vs property-matched vs selectivity-mined, reported separately; decoy-property distribution audit; AVE bias score reported for every split.
8. **Low-cost falsification.** Check whether measured selectivity pairs are separable by a linear probe on pretrained PLM pocket embeddings.
9. **Compute.** Low-medium.
10. **Zero-/few-shot.** Good for retrieval/screening; weaker for regression-quality ΔΔ.

---

### **D9 — Co-folding / structural priors as *features*, never as ground truth**  ★ Rank 9

1. **Hypothesis.** Predicted complexes (Boltz-2, AF3-class) supply interaction hypotheses that improve cold-target prediction.
2. **Why the protein alters the effect.** Predicted contacts are protein-specific.
3. **Minimum information.** Compute budget for co-folding at scale; a pocket definition.
4. **Representation / objective.** Use predicted contact maps / interface features as *inputs* to D1/D2; do **not** claim atomic recognition.
5. **Why not a target ID.** Contact features are pair-specific.
6. **Failure modes.** Serious. Co-folding models are trained on the PDB, so predicted features may re-import training-set memorization — a *leakage vector disguised as physics*. Independent evaluation suggests Boltz-2 affinity outputs may be substantially pose-insensitive ([arXiv:2603.05532](https://arxiv.org/html/2603.05532v1)), which would mean the affinity head is a chemogenomic prior, not a structural one.
7. **Negative controls.** Scrambled-pose features; features from a *different* ligand co-folded with the same protein; PDB-release-date-stratified evaluation to expose memorization; comparison against a sequence-only pocket baseline.
8. **Low-cost falsification.** On a held-out set with known complexes and *no* PDB entry before the model cutoff, test whether predicted-contact features add anything over the aligned pocket string.
9. **Compute.** **High** (co-folding inference at panel scale). This is the main reason it ranks below cheaper directions with the same or better evidential value.
10. **Zero-/few-shot.** Potentially strong for novel folds; currently unvalidated for this purpose.

---

### **D10 — Curated cross-target matched-pair evaluation surface (new dataset artifact)**  ★ Rank 10

1. **Hypothesis.** No existing benchmark can answer Core Task 1; therefore build one.
2. **Why the protein alters the effect.** n/a — infrastructure.
3. **Minimum information.** ChEMBL/Papyrus + panels + `mmpdb` + KLIFS/GPCRdb alignments + PLINDER similarity metrics for split construction.
4. **Deliverable.** **CTMP** (Cross-Target Matched Pairs): every record = (τ, ligand pair, protein, within-assay Δ, document id, censoring flags, pocket-identity cluster id). Ships with pre-computed cold splits at multiple pocket-identity thresholds and a published noise ceiling per stratum.
5. **Why not a target ID.** n/a.
6. **Failure modes.** Curation errors become field-wide; scope creep.
7. **Negative controls.** Ship the control splits (shuffled/random/similarity-matched) *inside* the release so they cannot be skipped.
8. **Low-cost falsification.** Pilot on kinases only before generalizing.
9. **Compute.** CPU; engineering-heavy.
10. **Zero-/few-shot.** Enables everything else. Highest long-term leverage, lowest immediate scientific novelty.

---

### **D11 — Prior-inverted and counterfactual protein splits (evaluation innovation)**  ★ Rank 11

1. **Hypothesis.** Borrowing from VQA-CP: construct cold-target splits where family-level SAR priors are *inverted* between train and test.
2. **Why the protein alters the effect.** n/a — diagnostic.
3. **Minimum information.** Family annotations + enough pairs to invert priors while preserving marginals.
4. **Objective.** Report performance on standard cold split *and* prior-inverted split. The gap quantifies prior-riding.
5. **Why not a target ID.** Directly measures the ID/prior component.
6. **Failure modes.** Inverted splits may be so small or so extreme that all methods fail, yielding no discrimination.
7. **Negative controls.** Marginal-matched random inversion.
8. **Low-cost falsification.** Construct on Davis panel in one day.
9. **Compute.** Trivial.
10. **Zero-/few-shot.** Diagnostic only — but it is the single most decisive test of the shortcut hypothesis I know of, and it does not appear to be in use in this field.

---

### Selections required by the brief

| Slot | Direction |
|---|---|
| Highest-probability | **D1** — protein-conditioned MMP transformation transfer (ΔΔ) |
| High-risk / high-reward | **D3** — weakly supervised latent binding-region discovery (MIL) with pre-registered external pocket validation |
| Data-centric | **D5** (+ **D10**) — variant/ortholog panels, then the CTMP evaluation surface |
| Representation-centric | **D2** — alignment-anchored residue vocabulary + fragment cross-attention |
| Training-objective innovation | **D6** — explicit additive/interaction decomposition with orthogonality + conditional-information + gradient-conflict control |
| Minimal program-terminating diagnostic | **T0-KILL** (§7) — cross-panel reproducibility of the interaction residual |

---

## 5. Recommended first experiment — the Cold-Target ΔΔ Sign Test (CT-ΔΔ)

**One sentence.** On single-platform kinase panels, test whether the *interaction residual* of ligand-pair potency differences is predictable on kinase groups held out at the pocket-identity level, using a model too small to memorize.

### Design

- **Datasets (kept strictly separate, never merged):** Davis 2011 (Kd) as primary; Metz 2011 (pKi) as independent replication; Klaeger 2017 kinobeads (apparent Kd) as a third, mechanistically distinct replication. Anastassiadis 2011 (% inhibition) used for **ordering only**.
- **Protein representation:** KLIFS 85-position aligned pocket → one-hot AA + Z-scales; optional pocket-restricted ESM-2 embeddings as a second arm.
- **Ligand representation:** ECFP4 (2048 bits) + MMP transformation encoding where available.
- **Estimand:** for every ligand pair (l₁,l₂) and protein p, compute `Δ(l₁,l₂|p)`, then the residual `R = Δ(l₁,l₂|p) − mean_p' Δ(l₁,l₂|p')` where the mean is taken **over training proteins only**.
- **Model:** kernel ridge regression with product kernel `k_pair(ligand pair) ⊗ k_pocket(protein)`. Deliberately low-capacity. **No neural network at this stage** — a null result from a big model is uninterpretable (it could be optimization failure), whereas a null from kernel ridge on a well-conditioned problem is informative.
- **Splits:** leave-one-kinase-group-out (Manning groups: TK, TKL, STE, CK1, AGC, CAMK, CMGC, atypical), *plus* the stricter constraint that no test kinase shares >50% identity over the 85 KLIFS pocket positions with any training kinase. 5 outer folds × 3 seeds.
- **Ligand-identity exclusion:** because panels reuse the same compounds across all proteins, ligand identity leaks by construction. Test pairs must involve at least one ligand absent from training — and a second, stricter arm where **both** are absent.
- **Censoring:** panel values reported as bounded (">10 µM", "<1 nM") make Δ a bound, not a value. Use a censored (Tobit) likelihood; **never impute**. Report the censored fraction per stratum.

### Primary endpoints

1. **Signed-difference accuracy** on `|R| > 1.0` log units (real cliffs), correct protein vs each control.
2. **Spearman(R̂, R)** on cold folds, reported as a fraction of the empirical noise ceiling `r_max`.
3. **Paired correct-minus-control effect**, bootstrapped **at the kinase-group level** (not the datapoint level).

### Control matrix (all run; all reported; none optional)

| Control | Construction | Expected if signal is real |
|---|---|---|
| Ligand-only | Drop protein features | Skill ≈ 0 **by construction** — a non-zero result proves a pipeline bug |
| Target-ID | One-hot protein | ≈ 0 on residual; large on raw affinity (demonstrates the confound) |
| Shuffled protein | Permute protein labels | ≈ 0 |
| Similarity-matched wrong protein | Substitute a paralog matched on pocket identity ±5% | Substantially degraded |
| Capacity-matched random (A) | iid Gaussian embeddings, same dim | ≈ 0 |
| Capacity-matched random (B) | Real embedding matrix, rows permuted **within family** | ≈ 0 — isolates taxonomy-riding |
| Residue-order permutation | Shuffle the 85 positions | Degraded |
| Composition scramble | Keep positions, randomize identities | Degraded |
| Column-shuffled alignment | Permute alignment columns consistently across all proteins | Degraded (tests positional semantics) |
| Document control | Held-out documents; same-doc vs cross-doc strata | Effect persists cross-document |
| Prior-inverted split (D11) | Invert family-level SAR priors | Effect persists, attenuated |
| Antisymmetry check | Require `R̂(l₁,l₂) = −R̂(l₂,l₁)` | Violation ⇒ shortcut |

### Mandatory positive control (run first)

Before interpreting any cold-target result: on *warm* kinases, the same pipeline must recover **known** selectivity determinants — e.g. the gatekeeper effect on type-II inhibitors, or hinge-region substitutions. A pipeline that cannot recover established biology in the easy setting cannot falsify Core Task 1 in the hard one. **A null result is only interpretable after the positive control passes.**

### Cost

~3 person-weeks; CPU only; no GPU. Deliverable is a decision, not a model.

---

## 6. Leakage and noise controls (complete checklist)

### Leakage

1. **Ligand identity** — exclude across splits; report the strict variant (both ligands novel).
2. **Near-duplicate ligands** — ECFP4 Tanimoto ≥ 0.9, MCS overlap, and InChIKey-skeleton matching to catch salts, tautomers, stereoisomers, isotopologues.
3. **Series leakage** — cluster by Murcko scaffold × document; a congeneric series counter-screened across targets is a single unit of evidence, not hundreds.
4. **Pocket-level protein leakage** — ⚠️ **the most commonly missed control.** Global sequence identity can be low while the binding pocket is nearly identical. Split on **pocket identity** (KLIFS 85 positions / GPCRdb generic positions / PLINDER pocket similarity), not full-sequence identity.
5. **Family and fold splits** — Manning groups, GPCR classes, ECOD/CATH for cross-family claims.
6. **Ortholog blocking** — human/rat/mouse orthologs of the same target must never straddle a split.
7. **Document/assay leakage** — no document in both train and test for paired contrasts; report same-doc vs cross-doc strata separately.
8. **Structural-prior leakage** — if using co-folding features, stratify by PDB release date relative to the structure model's training cutoff.
9. **Temporal realism** — a time-split arm as a secondary check (SIMPD-style simulated time splits are a reasonable proxy).
10. **Benchmark contamination** — Klaeger's clinical drugs are extensively represented in ChEMBL; treat "independent replication" claims with suspicion unless compound overlap is quantified and reported.

### Noise

11. **Per-stratum noise ceiling** from replicate and cross-panel agreement; report all results as a fraction of ceiling.
12. **Within-assay-only differences** for every Δ; cross-assay differences are inadmissible as primary evidence.
13. **Censoring** handled with a censored likelihood; censored fraction reported.
14. **Missing-not-at-random** — panels are near-dense (good); ChEMBL is not (absence ≠ inactivity). Never construct implicit negatives from absence in the primary analysis.
15. **Unit of analysis** — cluster bootstrap at the protein-cluster level, with random effects for protein cluster, scaffold, and document. Per-datapoint p-values are inflated by orders of magnitude here and are the most common statistical error in this literature.
16. **Multiplicity** — pre-register the primary endpoint; everything else is exploratory and labelled as such.
17. **Seeds** — 3 seeds for screening, 5+ for confirmation; report per-seed values, not just the mean.
18. **Modality-collapse monitoring** — log per-branch gradient norms and per-branch generalization gaps. **A null result with a collapsed protein branch is an optimization failure, not a falsification of Core Task 1.** This distinction must be enforced before any negative conclusion is drawn.

---

## 7. T0-KILL — the minimal diagnostic that can terminate the program

**Question.** Is the interaction residual reproducible *between independent experiments*?

**Procedure.**
1. Identify (compound, kinase) pairs measured in ≥2 independent panels (Davis / Metz / Klaeger / PKIS overlap).
2. Within each panel independently, fit `f = μ + a(p) + b(l)` and extract residual `Î`.
3. Correlate `Î_panelA` against `Î_panelB` over the shared overlap.

**Interpretation.**
- **Spearman ≥ 0.4** with n ≥ 200 → the interaction term is real and reproducible. Proceed. This is Core Task 1's Q1 answered affirmatively **without training anything**.
- **0.2 – 0.4** → weak but present; proceed with reduced scope and recalibrated effect-size expectations.
- **< 0.2 with upper CI < 0.3** → **STOP.** The interaction residual is not reproducible between experiments performed on the same molecules and the same proteins. No model can learn what two laboratories cannot agree on. Core Task 1 is **unresolvable on public affinity data**, and the correct response is to redirect resources toward data generation (a dedicated cross-target matched-pair panel) rather than toward modelling.

**Cost:** ~1 day, one analyst, no GPU. **Run this before anything else in the program.**

---

## 8. Staged research program and decision tree

```
STAGE 0 — DATA AUDIT & POWER  (3–4 weeks, CPU only)
├─ 0.1 T0-KILL: cross-panel reproducibility of interaction residual
│     └─ FAIL → TERMINATE modelling track; pivot to data generation
├─ 0.2 Variance components: Var[a], Var[b], Var[I], Var[e] per dataset
├─ 0.3 Excitation audit: # protein clusters (≤50% pocket id) × # shared ligands
│     └─ FAIL (identifiability impossible) → expand data before proceeding
├─ 0.4 Leakage audit: pocket-identity clustering, document/target MI, AVE bias
├─ 0.5 Noise ceiling per stratum; censoring rates
└─ 0.6 Power: ≥80% power for r=0.2 at protein-CLUSTER level
       GATE 0 → all of: Var[I] ≥ 25% of explainable variance (cluster-bootstrap
                95% CI excluding 15%); T0-KILL Spearman ≥ 0.2; ≥30 independent
                protein clusters with ≥20 shared ligands each; power ≥ 80%.

STAGE 1 — IDENTIFIABILITY  (4–6 weeks, CPU only, NO large model)
├─ 1.0 Positive control: recover known selectivity determinants on warm targets
│     └─ FAIL → pipeline is broken; fix before interpreting anything
├─ 1.1 CT-ΔΔ sign test on cold kinase groups (kernel ridge)
├─ 1.2 Full control matrix (§5)
├─ 1.3 Replication on a second, independent label system (Metz or Klaeger)
└─ 1.4 D5 variant/ortholog ANOVA as an orthogonal existence proof
       GATE 1 → ALL of:
         (a) cold-fold residual Spearman ≥ 0.25 as point estimate, and ≥ 0.15
             for EVERY seed (3 seeds);
         (b) signed-difference accuracy ≥ 60% on |ΔΔ| > 1 log unit;
         (c) every control within [−0.05, +0.05] Spearman, CI covering 0;
         (d) correct-minus-control paired difference ≥ 0.20 with cluster
             bootstrap 95% CI excluding 0;
         (e) effect persists cross-document and in the strict ligand-exclusion arm;
         (f) effect NOT carried by ≤2 protein families (drop-one-family analysis);
         (g) replication in ≥1 independent label system at ≥60% of primary effect;
         (h) antisymmetry satisfied.
       ── ONLY GATE 1 AUTHORIZES MODEL-SCALE TRAINING ──

STAGE 1.5 — REPRESENTATION RECOVERY  (3–4 weeks)
├─ D2 sparse positional probe → which aligned positions carry signal?
├─ D3 ORACLE test: does knowing the true pocket help at all?
└─ D6 decomposition at matched rank, real vs random embeddings
       GATE 1.5 → a sparse, chemically interpretable position set is identified,
                  AND the oracle-pocket arm beats whole-sequence input.
                  FAIL → the signal exists (Gate 1) but current representations
                  cannot recover it → representation research, NOT bigger models.

STAGE 2 — MODEL-SCALE TRAINING  (6–10 weeks, GPU)
├─ D2 cross-attention + D6 objective, per label system
├─ Full cold-target DTA metrics: RMSE, Pearson, Spearman, CI
├─ Activity-cliff strata (MoleculeACE-style definitions, per target)
├─ Prior-inverted split (D11)
└─ Modality-collapse and gradient-conflict monitoring throughout
       GATE 2 → cold-target RMSE improvement ≥ 0.10 log units over
                (ligand-only + protein-main-effect) with cluster-bootstrap CI
                excluding 0, 5 seeds; CI improvement ≥ 0.02; controls flat;
                AND the improvement is larger on the activity-cliff stratum than
                on the non-cliff stratum (the signature of genuine interaction
                learning rather than smoothing).

STAGE 3 — FEW-SHOT / META-LEARNING  (6–8 weeks)
├─ D4 CNP with k = 0, 1, 4, 16, 64 support anchors
├─ Ablation: support-only vs support + protein features at each k
└─ Prospective-style temporal holdout
       GATE 3 → non-zero k=0 skill with controls flat, AND a persistent
                (non-converging) gap between support-only and support+protein
                out to k = 16.
```

---

## 9. Conditions under which Core Task 1 is SOLVED, UNRESOLVED, or FALSIFIED

### SOLVED

All of the following, simultaneously:

1. **Existence.** T0-KILL Spearman ≥ 0.4, or Var[I] ≥ 25% of explainable variance with cluster-bootstrap CI excluding 15%, in ≥2 independent datasets.
2. **Recovery on cold targets.** Interaction-residual Spearman ≥ 0.25 (5 seeds, all ≥ 0.15) on protein clusters held out at ≤50% pocket identity, expressed also as ≥ 40% of the empirical noise ceiling.
3. **Signed alignment, not magnitude change.** Signed-difference accuracy ≥ 60% on |ΔΔ| > 1 log unit. A large prediction shift under protein perturbation that does not align with signed affinity differences counts as **zero evidence**.
4. **Controls flat.** Shuffled, capacity-matched random (both flavours), residue-permuted, composition-scrambled, and column-shuffled branches all within [−0.05, +0.05] with CIs covering 0. Similarity-matched wrong protein substantially degraded.
5. **Self-improvement, not control-degradation.** The correct-protein model must improve on **absolute** cold-target metrics against ligand-only + protein-main-effect (RMSE ≥ 0.10 log units, CI excluding 0). "Correct beats corrupted" while both lose to ligand-only is a **failure**, not a partial success.
6. **Breadth.** Not carried by ≤2 protein families; survives drop-one-family analysis; replicates in ≥2 label systems and ideally ≥2 superfamilies (kinase + GPCR).
7. **Assay-independence.** Persists cross-document and under the prior-inverted split.
8. **Utility.** Gate 2 passed — i.e. it actually improves cold-target DTA, with the cliff-stratum improvement exceeding the non-cliff improvement.

### UNRESOLVED (the most likely honest outcome)

Any of:
- Existence established (Gate 0/1 pass) but Gate 1.5 fails → **signal is present, representations cannot yet recover it.** Correct response: representation research and better data. **Explicitly not**: a larger neural model. Scaling in response to a failed identifiability test converts an interpretable null into an uninterpretable one.
- Effect present but confined to one superfamily → **restate the claim as family-conditional transferability** and re-scope the product to within-superfamily zero-shot + cross-superfamily few-shot.
- Effect present but below the noise ceiling for practical use → data-generation problem, not a modelling problem.
- Any Stage-1 null obtained with a collapsed protein branch, a failed positive control, or insufficient excitation → **not a falsification**; the experiment did not test the hypothesis.

### FALSIFIED (on these data)

All of:
1. T0-KILL fails: interaction residuals do not reproduce across independent panels (Spearman < 0.2, upper CI < 0.3, n ≥ 200); **and**
2. Stage 1 shows no cold-target residual signal (Spearman CI covering 0) across ≥2 datasets, ≥2 protein representations (aligned pocket + PLM), and ≥3 seeds; **and**
3. The positive control passed (known determinants recovered on warm targets) — proving the pipeline works; **and**
4. Protein-branch gradient norms confirm no modality collapse; **and**
5. The Stage-0 excitation audit confirms the design was identifiable in principle.

**Scope of falsification.** Even then, the correct claim is narrow and must be stated as such: *"transferable protein-conditioned interaction signal is not recoverable from currently available public affinity data at the effect sizes and noise levels measured here."* That is a statement about public data, not about biology. Selectivity exists; if we cannot measure it through this lens, the deficiency is in the data or the lens.

---

## 10. What I would refuse to do, and why

- **Scale up after a failed identifiability test.** A bigger model that "works" after a Stage-1 null has almost certainly found a shortcut the small model was too weak to exploit. Capacity is not evidence.
- **Report attention maps as biological explanation.** D3 pre-registers an external pocket-recovery AUC with a stated threshold precisely so that interpretability claims are falsifiable rather than illustrative.
- **Claim atomic or Cartesian recognition.** KLIFS alignment supports *positional* (topological) claims about aligned residues. Predicted complexes supply hypotheses. Neither licenses a claim about Cartesian recognition without legal common-frame complex coordinates.
- **Use KIBA as primary evidence,** or merge IC50/Ki/Kd into a single regression target. The harmonization is not scientifically defensible; the defensible one is within-assay differencing.
- **Report datapoint-level significance.** The unit of analysis is the protein cluster. With ~30 clusters, effective N is ~30, not ~10⁵.
- **Present an exploratory screen as confirmatory.** Every table in the final report carries an explicit label: pre-registered primary, pre-registered secondary, or exploratory.

---

## 11. Immediate next actions (first 30 days)

| Day | Action | Output |
|---|---|---|
| 1–2 | T0-KILL cross-panel reproducibility | Go / no-go on the entire program |
| 3–7 | Variance-component + excitation audit on Davis, Metz, Klaeger | Var[I] estimates with CIs; identifiability verdict |
| 8–12 | Build pocket-identity clustering from KLIFS; construct cold splits + all control splits | Reusable split artifact (CTMP v0) |
| 13–18 | Warm-target positive control (gatekeeper/hinge determinants) | Pipeline validation |
| 19–28 | CT-ΔΔ sign test + full control matrix, 3 seeds | Gate 1 decision |
| 29–30 | Write-up with pre-registered Stage 2 protocol | Authorization document (or termination memo) |

**Total to a defensible Gate-1 decision: ~6 weeks, CPU only, one researcher.** No model-scale training is authorized before that decision.

---

### Verification notes

Panel dimensions quoted in §3 (Davis ≈72×442, Anastassiadis ≈178×300, Metz ≈3,800×172, Klaeger 243 drugs) are from recall and must be re-confirmed against the primary supplementary files during Stage 0.1 — the exact counts affect the power calculation. Licensing terms for BindingDB, PDBbind, Papyrus, KLIFS and PLINDER should likewise be re-verified at the source before any redistribution of derived artifacts; the table records access route and known constraints, not legal advice.
