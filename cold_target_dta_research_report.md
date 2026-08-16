# Cold-Target Few-Shot DTA: Method Research Report

**Scope & method.** First-party sources only (papers, official code repos, official docs), gathered via targeted web searches. Where a 2025–26 paper is paywalled (AdaMBind, CrossLinker, LigUnity details), mechanism statements are flagged as abstract-level. **Constraints assumed throughout:** BindingDB pK regression; protein = ESM-2 pooled + per-residue tokens (128 slots); ligand = molecular graph (GINE); no 3D complex coordinates for any deployment pair; strict cold-target episodic few-shot (k=0, 1, 2, 3, 5; disjoint support/query); single-stage episodic training, one optimizer update, no MAML inner loop, no ridge/closed form, no test-time gradients; frozen precomputed features allowed if reported as external data. Failure modes to fix: (1) zero-shot within-target constant collapse (CI ≈ 0.52–0.58, worst on activity cliffs); (2) few-shot collapse to a scalar level shift; (3) CI↔MSE tension (MSE gains = shrinkage to target mean).

---

## 1. ActFound — pairwise meta-learning bioactivity foundation model

**Paper/code:** [bioRxiv 2023.10.30.564861](https://www.biorxiv.org/content/biorxiv/early/2023/11/02/2023.10.30.564861.full.pdf), [GitHub BFeng14/ActFound](https://github.com/BFeng14/ActFound), independent [reusability report (Nat. Mach. Intell. 2026)](https://link.springer.com/article/10.1038/s42256-026-01187-y).

**Core mechanism.** ActFound does *not* use MAML, and does not use anchor-based contrastive mining across assays. A task is one assay (a set of compounds measured in one assay); training is a **within-assay pairwise comparison**: given two compounds from the same assay, the model predicts which is more active. The molecule encoder is pretrained at scale (~1.6M bioactivity records; press report [Sina/Beijing Univ.](https://finance.sina.cn/2024-08-16/detail-inciuytn5935615.d.html)) and the pairwise head is meta-trained across thousands of assays. Unseen assays are handled in zero-shot by scoring all query compounds against each other pairwise; with a few labeled compounds the same pairwise machinery becomes an episodic few-shot task.

**Why it works.** Pairwise/relative formulation is invariant to per-assay label calibration (each assay has its own distribution, scale, and noise). The model never commits to a global absolute scale; it learns the molecular features that *order* compounds within a target/assay — exactly the quantity that decides CI. Reported results show strong zero-shot and few-shot gains on bioactivity benchmarks, with few-shot gains scaling with support size; the reusability report re-evaluates generalization on antibacterial natural products (headline numbers should be read with care, but the formulation is methodologically decisive).

**Data/training.** ChEMBL-scale multi-assay corpus; episodic task sampling per assay; pairwise loss; no 3D structures. Single forward pass per episode — no inner loop.

**Transferable elements.**
- Replace (or co-train) absolute pK regression with **within-target pairwise ranking** as the primary objective — additive per-target constants are irrelevant to pair order, directly attacking failure (1).
- Episode = target's support+query ligands; pairwise loss over query pairs and support–query pairs; single-stage, one update — satisfies our constraints.
- The "assay invariance" insight is the same disease as our k=1 scalar level shift: a pairwise loss has no gradient toward a constant shift (a constant output gives chance pairwise accuracy, so the optimizer is forced to produce ligand-dependent variation).

## 2. FS-CAP, AdaMBind, MetaDTA — per-target adaptation family

**FS-CAP** — [arXiv 2311.16328](https://arxiv.org/abs/2311.16328), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10705577/).
*Core:* "target-free" few-shot compound activity: an assay is a task; **no protein features at all**. Episodic meta-training over assays; a molecular GNN (GINE-class) encoder embeds compounds; a query's activity is predicted by **attention-weighted aggregation over the support set** (soft addressing over support embeddings and their labels, iteratively refined, Matching-Networks-style). Because weights depend on each query's own embedding, the correction is a per-query function of the support — never a scalar.
*Transferable:* support-conditional attention as the few-shot mechanism; assay-as-task episodic sampling; disjoint support/query.

**AdaMBind** — [Nat. Commun. 2026](https://www.nature.com/articles/s41467-026-70554-5) (abstract/figures level; paywalled). *Core:* meta-learning plus a **task-adaptive module** that modulates the DTA predictor per target/assay; few-shot adaptation demonstrated on BindingDB-type benchmarks. *Transferable:* per-target modulation without inner loops — one forward pass computing target-conditioned parameters, the same role as FiLM/adaLN (Section 8).

**MetaDTA** — [KAIST record](https://koasas.kaist.ac.kr/handle/10203/299261?mode=full), [Semantic Scholar](https://www.semanticscholar.org/paper/MetaDTA%3A-Meta-learning-based-drug-target-binding-Lee-Yoo/3861669e900f88cf5127850d9031da610216d84a). MAML-based DTA meta-learning with ESM protein features on leave-target-out episodic splits; showed few-shot gains on unseen targets. The MAML inner loop is excluded by our constraints, but its **episodic cold-target task construction** (train tasks = targets held out entirely; support/query disjoint) is the evaluation protocol we should adopt.

**Transferable elements.** Episodic leave-one-target-out task sampling (MetaDTA); attention-over-support for per-query correction (FS-CAP); task-conditioning modules instead of MAML (AdaMBind); a target-free ablation is a useful control to isolate whether protein conditioning is doing ranking work or just adding a constant.

## 3. PSICHIC — sequence-only protein–ligand interaction

**Paper/code:** [Nat. Mach. Intell. 2024](https://www.nature.com/articles/s42256-024-00847-1), [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.17.558145v1.full), [GitHub huankoh/PSICHIC](https://github.com/huankoh/PSICHIC).

**Core mechanism.** Inputs are protein sequence + ligand SMILES — no 3D pose at inference. Protein residues are embedded with ESM-2 and arranged into a learnable **"physicochemical" graph** (node/edge features carrying learnable physicochemical descriptors); the ligand is a molecular graph; ligand atoms and protein residues interact through **attention**, producing a pair-specific interaction fingerprint and an attention-based estimate of which residues contact the ligand.

**Why it works.** (i) residue-level (not just pooled) protein tokens let the model localize interacting residues; (ii) ligand↔protein cross-attention builds a pair-specific fingerprint rather than a bag-of-features interaction; (iii) contrastive self-supervised pretraining on a large corpus, then supervised affinity fine-tuning. It reported state-of-the-art on several DTA benchmarks (including cold/unseen targets) and correct residue-level interaction prediction without 3D.

**Data caveat.** Its *pretraining* exploited protein–ligand complex data (3D-derived) — we have no 3D for deployment pairs, so only the architecture is directly reusable; the pretraining corpus would count as external data if we wanted to replicate the contrastive stage.

**Transferable elements.** Per-residue ESM-2 tokens + GINE ligand + **cross-attention interaction fingerprint** as our encoder skeleton; attention-based residue weighting (cheap, uses our 128 residue slots); contrastive pretraining objective over affinity pairs; proof that sequence+2D can beat pose-based methods on ranking tasks.

## 4. ConPLex — contrastive DTA

**Paper/code:** [PNAS 2023](https://www.pnas.org/doi/10.1073/pnas.2220778120), [GitHub samsledje/ConPLex](https://github.com/samsledje/ConPLex), [MLSB 2022 (decoys)](https://www.mlsb.io/papers_2022/Contrasting_drugs_from_decoys.pdf).

**Core mechanism.** Two encoders — frozen ESM-2 protein embeddings (pocket/full sequence) and a trainable SMILES encoder — projected to a shared space, trained with an **InfoNCE-style contrastive loss** where each drug is the *anchor*, its target is the positive, and other targets in the batch are negatives (a "drug-anchor" formulation for stability; variants add decoys as hard negatives). Zero-shot inference = cosine similarity between drug and protein embeddings. No 3D, no docking.

**Why it works.** Contrastive alignment makes the embedding space encode interaction propensity rather than absolute affinity, generalizes across targets, and the anchor formulation + hard negatives calibrate the score. Reported strong DTI/virtual-screening performance including generalization to unseen protein targets.

**Data/training.** BindingDB-scale interaction data; single-stage contrastive training (no inner loop) — compatible with our constraints; frozen ESM-2 is external data.

**Transferable elements.** Frozen ESM-2 + ligand encoder aligned in one contrastive space; zero-shot interaction scoring by dot product (a natural k=0 output before a regression head); in-batch negatives + decoys; note that ConPLex is binary (interact/no), so affinity *ranking* needs its ranking machinery layered on top — but the representation it induces is cliff-relevant (nearest neighbors in this space are the "similar ligands").

## 5. kNN-DTA / NN-DTA / Ada-kNN-DTA — retrieval-augmented DTA

**Paper:** [arXiv 2407.15202](https://arxiv.org/abs/2407.15202) / [KDD 2024](https://dl.acm.org/doi/abs/10.1145/3627673.3679704).

**Core mechanism.** Non-parametric: frozen pre-trained encoders (ESM-2 protein, molecular embeddings) define an embedding space; for a query (drug, target), retrieve the k nearest neighbors from the training corpus and predict affinity as a **similarity-weighted combination of neighbor labels** (NN-DTA); Ada-kNN-DTA adaptively selects k per query. Reported strong results on BindingDB/DAVIS **cold-target splits** and robustness when labeled data are scarce.

**Why it works.** Label pooling over neighbors is inherently per-query (weights differ across queries), so predictions cannot collapse to a constant; retrieval exploits the large weakly-labeled corpus; in cold-target regimes where parametric models shrink to the mean, a non-parametric estimate anchored on real measured neighbors degrades gracefully.

**Data/training.** The kNN core needs no gradient training — only the embedding space (precomputed, external). This matches our frozen-feature constraint exactly.

**Transferable elements.** Similarity-weighted neighbor-label aggregation as the few-shot correction (fixes failure (2): weights are functions of each query's embedding); for k=0, retrieve **same-target ligands from the training corpus** (BindingDB has many ligands per target) — a non-collapsing zero-shot baseline; if the corpus lacks the cold target, pool ligands from *similar targets* by pooled-ESM-2 distance (CrossLinker-style relational prior, Section 10); adaptive k; report frozen ESM-2/GINE features as external data.

## 6. PBCNet / PBCNet2.0 — Siamese relative binding affinity

**Paper:** [Nat. Comput. Sci. 2023 (PBCNet)](https://www.nature.com/articles/s43588-023-00529-9); [bioRxiv 2025.06.04.657800 (PBCNet2.0)](https://www.biorxiv.org/content/10.1101/2025.06.04.657800), [Nat. Chem. Biol. 2025 (PBCNet2.0)](https://pubmed.ncbi.nlm.nih.gov/42286270/).

**Core mechanism.** PBCNet is a **Siamese network over pairs of protein–ligand complexes (with 3D poses)** that predicts *relative* binding affinity between two ligands of the same protein, trained with a pairwise ranking/relative-regression objective; designed for lead-optimization ranking of congeneric series. PBCNet2.0 upgrades the geometry encoder to Cartesian tensors and targets pharmacologic probe discovery with better ranking accuracy.

**Constraint.** Both require binding poses / atomistic 3D per pair → **not usable in our setting (0/17,717 pairs)**.

**Transferable elements.** The **relative-recognition idea** — Siamese pair comparison + pairwise ranking loss over same-target ligand pairs — transfers wholesale to sequence+graph inputs; evaluate relative predictions by ranking accuracy (not RMSE); the "compare two ligands of the same protein" task formulation is exactly our within-target CI problem.

## 7. Activity cliffs (ACs) and molecular matched pairs (MMPs)

**Key refs:** [ICML 2024 — Activity Cliff-Aware Molecular Representations](https://icml.cc/virtual/2024/38085) (+ [JKU record](https://research.jku.at/de/publications/towards-learning-activity-cliff-aware-molecular-representations/)); [Activity Cliff-Informed Contrastive Learning (2023)](https://sciety.org/articles/activity/10.21203/rs.3.rs-2988283/v3) (+ [PMC PDF](https://pmc.ncbi.nlm.nih.gov/articles/PMC11643338/pdf/nihpp-rs2988283v2.pdf)); [Activity-cliff awareness for robust graph learning (ChemRxiv 2023)](https://chemrxiv.org/doi/full/10.26434/chemrxiv-2023-5cz7s/v5); [Property Cliffs Reveal Hidden Errors (arXiv 2026)](https://ar5iv.labs.arxiv.org/html/2605.17265); [MolFeSCue (few-shot + contrastive, 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10984949/).

**Core ideas.** An AC is a pair of similar ligands (Tanimoto-similar, often an MMP differing by one substituent) with a large activity gap (≥1 log unit in our case). Models systematically underperform on cliffs. Methods: (i) **AC-aware contrastive representations** — cliff partners are treated as hard-negative pairs so the embedding separates them; (ii) **AC-aware training weighting** in graph models; (iii) **cliff-stratified evaluation** — reporting performance separately on cliff-rich subsets reveals model blind spots that pooled metrics hide.

**Transferable elements.** Mine within-target MMP/similar pairs with |ΔpK| ≥ 1.0 (support–support, support–query, query–query) and **oversample them in the pairwise ranking loss** — the direct fix for failure (1) on cliffs; add a margin/contrastive term pulling cliff partners apart in the (GINE) ligand embedding space; stratify all CI/Spearman reporting by cliff membership so "barely above chance" is measured, not averaged away.

## 8. Few-shot regression / meta-learning toolbox

**Key refs:** [Matching Networks (NIPS 2016)](https://arxiv.org/abs/1606.04080); [Prototypical Networks (NIPS 2017)](https://arxiv.org/abs/1703.05175); [Set Transformer (ICML 2019)](https://arxiv.org/abs/1810.00825); [FEAT (CVPR 2020)](https://arxiv.org/abs/1812.03664); [CNP/NP (ICML 2018)](https://arxiv.org/abs/1807.01613), [NP](https://arxiv.org/abs/1807.01622), [ANP](https://arxiv.org/abs/1808.03856); [Hypernetworks](https://arxiv.org/abs/1609.09106); [FiLM (ICLR 2018)](https://arxiv.org/abs/1709.07871); [DiT/adaLN (ICLR 2023)](https://arxiv.org/abs/2212.09748); [Context-enriched molecule representations for few-shot drug discovery (2023)](https://arxiv.org/abs/2305.09481).

**Per method, in one line each.**
- *Matching Networks:* soft-attention over support labels (per-query weights) — the canonical support-conditional predictor; exactly the mechanism that prevents a scalar shift.
- *Prototypical Networks:* support mean as prototype; in regression the prototype degenerates to the **mean label** — a cautionary tale: mean-based conditioning reproduces failure (2).
- *Set Transformer:* permutation-invariant set encoding (ISAB + pooling-by-multihead-attention) — the right primitive for encoding a support set into a fixed-size context vector.
- *FEAT:* a set-to-set transformer transforms support prototypes *and* the query embedding with the same task-adapted function — task-adaptive **embeddings**, not scalar shifts.
- *CNP/NP/ANP:* encode {(x_i, y_i)} support pairs into a context vector; decode per-query mean/variance conditioned on x_query — purpose-built few-shot regression with continuous labels; NP adds a latent z, ANP adds attention (the best fit for our setting).
- *Hypernetworks / FiLM / adaLN:* context → weights or (γ, β) feature-wise modulation (adaLN = scale/shift from a conditioning embedding, as in diffusion transformers) — cheap per-target modulation with no inner loop.
- *Context-enriched molecular representations (Schimunek et al.):* applies attention/set-based context enrichment to *molecular* embeddings for few-shot activity prediction — a molecular-domain FEAT; directly on our problem.

**Transferable elements.** Attention-over-support (Matching/ANP/FS-CAP) for per-query correction; Set Transformer-style encoders for support-set context; CNP-style **(x, y) support encoding with labels as input** — labels enter the context, so the model can exploit them without test-time gradients; FiLM/adaLN for target conditioning; avoid pure label-means (ProtoNet regression) — it is the collapse mode we already observe.

## 9. Learning-to-rank for regression targets

**Key refs:** [RankNet (ICML 2005)](https://www.microsoft.com/en-us/research/publication/learning-to-rank-using-gradient-descent/); [LambdaRank (NIPS 2006)](https://papers.neurips.cc/paper_files/paper/2006/hash/af44c4c56f385c43f2529f9b1b018f6a-Abstract.html); [ListNet](https://arxiv.org/abs/0712.0241); [LambdaLoss (CIKM 2018)](https://arxiv.org/abs/1802.02265); [allRank implementation](https://github.com/vishalbelsare/allRank); [listwise loss + transfer for DTI top-k ranking (NTHU thesis)](https://etd.lib.nycu.edu.tw/cgi-bin/gs32/hugsweb.cgi?o=dnthucdr&s=id=%22G021070644250%22).

**Core ideas.** *RankNet:* pairwise logistic loss on (i, j) pairs — smooth, but weights all pairs equally. *LambdaRank:* rescales each pair's gradient by the NDCG delta of swapping the pair — optimizes the ranking metric directly (the natural choice when CI is the metric). *ListNet:* listwise KL divergence over "top-one" probabilities; *LambdaLoss:* a unified framework showing LambdaRank is the metric-optimizing gradient and deriving principled listwise variants. **Crucial property for regression targets:** all ranking losses are invariant to per-task additive constants — the level shift that currently destroys our few-shot ranking has *zero gradient* under a ranking objective. Property/DTA applications: LigUnity (Section 10) trains with a within-assay ranking loss; the NTHU thesis applied listwise losses to DTI top-k.

**Transferable elements.** RankNet (pairwise logistic) as the base few-shot objective; LambdaRank / 1-log(k+1) weighting to prioritize top-of-list (high-affinity) queries; ranking losses decouple CI from MSE (failure (3)) because their gradient ignores the mean; pair a ranking loss with a separate MSE level head so MSE shrinkage can't corrupt ordering (see Candidate 3).

## 10. 2024–2026 DTA work improving within-target ranking

**LigUnity / Hierarchical affinity landscape** — [bioRxiv 2025.02.17.638554](https://www.biorxiv.org/content/10.1101/2025.02.17.638554v2.full), [Cell Patterns 2025](https://www.sciencedirect.com/science/article/pii/S2666389925002193), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12546767/). Foundation model jointly optimizing virtual screening and hit-to-lead; **within-assay ranking loss with pair weights ∝ 1/log(k+1)** to prioritize top-ranked molecules; its pocket encoder needs pocket-level (3D) inputs — flag: the ranking-loss formulation transfers, the input modality does not.

**CrossLinker** — [JCIM 2025](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03216). Explicitly targets **cold-start and few-shot DTI** by aligning "relational" context (similar drugs/targets via a heterogeneous graph) with "sequential" context (SMILES + protein sequence) — a relational prior from similar known targets for cold targets.

**MetaDTA (2022)** — [KAIST](https://koasas.kaist.ac.kr/handle/10203/299261?mode=full) and **Contrastive Meta-Learning for DTA** — [IEEE 2022](https://ieeexplore.ieee.org/document/9995372): MAML/contrastive meta-learning variants with episodic cold-target evaluation; useful for task design, excluded by the no-inner-loop constraint. **AdaMBind** — [Nat. Commun. 2026](https://www.nature.com/articles/s41467-026-70554-5): meta-learning + task-adaptive module, few-shot DTA gains (abstract-level). **MolFeSCue** — [PMC 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10984949/): few-shot + contrastive for data-limited property prediction.

**Transferable elements.** Within-target ranking loss with top-of-list prioritization (LigUnity); relational/retrieval prior from similar targets (CrossLinker, kNN-DTA); episodic cold-target evaluation protocol (MetaDTA); task-conditioned modulation (AdaMBind).

---

## Top 5 concrete design candidates (ranked by fit to our failure modes)

**1. Within-target pairwise ranking (LambdaRank-weighted) as the primary episodic objective** — *from ActFound, PBCNet, LambdaRank, LigUnity.* Episode = one target; encode support+query ligands with GINE and the protein with pooled ESM-2; a small head scores each ligand. Loss = RankNet pairwise logistic over in-episode pairs (query–query and support–query), with cliff pairs oversampled and pair gradients weighted by NDCG-delta or 1/log(k+1) to prioritize top-ranked queries. k=0: pure ranking head — CI is the objective, a per-target constant is at chance, so the optimizer is forced to produce ligand-dependent variation (fixes failure (1)). k=1–5: support pairs join the same loss; the "correction" is inherent because support-involving pairs constrain the target's ordering — no scalar shift exists to learn (fixes failure (2)). Single-stage, one update, no inner loop. MSE is reported separately via Candidate 3's level head so ranking gradients never trade off against shrinkage (failure (3)).

**2. Support-attention per-query conditioning (Matching/ANP/FS-CAP/FEAT style)** — *from Matching Networks, ANP, FS-CAP, FEAT, Context-enriched molecules.* The query ligand embedding attends over support ligand embeddings (CNP-style: support tokens are (x, y) pairs so labels are visible without test-time gradients); the attended context is concatenated or FiLM-modulated into the prediction head. Each query gets a *distinct* context, so the k=1 correction is a function of the query's own structure — kills the constant-shift collapse (failure (2)) by construction. k=0 uses a learned default context or the pooled protein embedding as the "empty support"; the model still must produce within-target variation because training is episodic. One forward pass, standard backprop. Combine with #1 so the context is ranked by the pairwise loss (cliffs get different contexts than flat pairs).

**3. Level + residual two-head decomposition** — *from CNP (level from context) + ranking-trained residual (LambdaRank, ActFound).* Head L predicts the per-target mean (MSE; conditioned on pooled protein embedding + support labels via a set encoder); head R predicts each query's *within-target rank residual*, trained only by a pairwise ranking loss. Prediction = L + R·s, with s a global residual scale calibrated on validation. Because R's gradient comes from ranking, MSE can no longer "buy" gains by shrinking R to zero — CI and MSE are decoupled (failure (3)). k=0: L from the protein embedding alone, R from the ranking head. The two heads share the encoder but not the objective, so shrinkage to the target mean is confined to L, where it belongs.

**4. Cliff-aware hard-negative ranking and representation** — *from AC contrastive learning (ICML 2024, RS 2023) + ActFound pairwise.* Inside each episode, mine MMP/similar pairs (Tanimoto ≥ 0.6–0.7) with |ΔpK| ≥ 1.0 across support+query; (a) oversample these pairs in the pairwise ranking loss (weight ×(3–10)), and (b) add a margin contrastive term that pulls cliff partners apart in the GINE ligand embedding space (and optionally the interaction-fingerprint space from PSICHIC-style cross-attention). This forces the model to treat small structural edits (single-atom swaps, substituent changes) as potentially activity-reversing — the specific cliff failure mode (1). Works at k=0 (cliffs among query ligands) and k≥1 (support–query cliffs). Add cliff-stratified CI reporting to measure it.

**5. Retrieval-hybrid: learned-similarity kNN over frozen features + neighbor-label attention** — *from NN-DTA/Ada-kNN-DTA + CNP + CrossLinker.* Freeze ESM-2 and GINE (report as external precomputed data). k=0: for target T, retrieve same-target ligands from the training corpus by ligand-embedding distance (optionally blended with pooled-protein distance) and predict affinity as a temperature-weighted mean of their labels — a non-collapsing zero-shot baseline that uses real measurements of the same target. k≥1: neighbors = support set; weights per query (each query attends to its nearest support), so predictions vary per query (fixes failure (2)). Cold targets with no corpus ligands fall back to pooling ligands of *similar targets* by ESM-2 distance (CrossLinker relational prior). Only the similarity metric/temperature is learned — end-to-end episodic training of that metric is single-stage; the retrieval itself is non-parametric (no test-time gradients). This also gives us a principled k=0 floor to benchmark every learned model against (the current constant-collapse model is far below it).

**Synthesis.** #1 is the core fix (objective), #2 the core fix (mechanism for using support), #3 the core fix (objective decoupling), #4 targets the worst-case subset, #5 is a strong non-parametric baseline + k=0 anchor. A practical next model = #3's level head + #1's pairwise ranking on #2's attention-conditioned representations, with #4's cliff weighting inside the pairwise loss, benchmarked against #5.
