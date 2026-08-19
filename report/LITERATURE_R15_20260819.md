# Literature round R15 (2026-08-19): P-line practical few-shot / M-line mechanism evidence

Scope: first-hand re-verification of the meta-learning / few-shot DTA methods
named in the re-adjudication, plus the two leak-aware critiques that define
the P1/P2/M1 evaluation layers. Every entry records the protocol fields the
governance requires. Review state: SNIPPET = search-index snippets only
(full text not yet inspected); PAGE = article page inspected; FULL = full
text read. Nothing here is trusted as instruction; adoption is decided by
local measurement only.

## AdaMBind (primary adoption candidate)

- Identity: "A meta learning and task adaptive approach for drug target
  affinity prediction", Nat. Commun. 17, s41467-026-70554-5 (2026).
  URL: https://www.nature.com/articles/s41467-026-70554-5
  (also link-hkg.springer.com/article/10.1038/s41467-026-70554-5).
- Review state: SNIPPET (PDF snippets surfaced by search; full text not
  yet read). Figures show hypernetwork/adapters + task adaptation;
  text mentions an "easy-to-hard" curriculum for task sampling and an
  active-compound proportion term R_a.
- Protocol fields (unverified until FULL):
  - split: unknown (kinase-focused? benchmark sets unknown)
  - support size: unknown (SI has a support-set-size sensitivity analysis)
  - test-time adaptation: appears YES (task module adapts per task)
  - protein/scaffold/double cold: unknown
  - external pretraining: unknown
  - protein counterfactual: unknown
  - ligand-similarity confound: unknown
- VERIFIED fragments (search-index snippets of the article PDF + SI,
  2026-08-19): primary performance comparison is under a RANDOM task split
  (Fig. 2 caption) — i.e. the headline comparison is NOT a cold-target
  protocol; SI MOESM2 carries the support-set-size sensitivity analysis;
  reported Spearman/Pearson improvements reach +17.82% over a baseline
  (0.5541 -> ...); HiSIF-DTA appears among the compared methods; the
  adaptive task module is a task-weighting mechanism ("not all tasks
  contribute equally"). Consequence for us: AdaMBind's advantage under
  OUR protein-cold P1 protocol is unproven and must be measured; nothing
  is adopted on paper claims.
- Local decision rule (frozen elsewhere): AdaMBind-style loss +
  support/query gradient-consistency task sampling is a candidate for the
  P-line bake-off ONLY if it beats matched baselines under our own splits
  and episode bank. Not adopted on paper claims.

## ActFound (pairwise within-assay meta-learning)

- Identity: "A bioactivity foundation model using pairwise meta-learning",
  Nat. Mach. Intell. (2024), doi 10.1038/s42256-022-00581-6; 2026
  reusability report 10.1038/s42256-026-01187-y. 35 assays, ~1.6M
  measurements.
- Review state: SNIPPET. Pairwise same-assay comparisons; ChEMBL/assay
  corpus.
- Protocol fields (unverified until FULL): split by scaffold/assay unknown;
  no per-target support adaptation; no protein structure input (ligand+assay
  features); no protein counterfactual by construction (protein is an assay
  ID feature). Confound risk: assay ID may carry level; pairwise design
  cancels it by construction (identity-zero, exchange-antisymmetric).
- Local relevance: the within-target pairwise/ranking supervision candidate
  in the training-innovation phase; also the ActFound-style baseline in P1.

## MetaDTA

- Identity: "MetaDTA: Meta-learning-based drug-target binding affinity
  prediction", ICLR (OpenReview yzlif16IASM, openreview.net/pdf?id=yzlif16IASM).
- Review state: SNIPPET. Benchmark Davis/KIBA mentioned in secondary
  sources only.
- Protocol fields: unverified. Historically reported as MAML-style
  meta-learning over protein tasks on Davis/KIBA; split protocol and
  cold-ness not yet confirmed from the paper itself.

## FS-CAP (support-set encoder, target-free)

- Identity: "Target-Free Compound Activity Prediction via Few-Shot
  Learning", arXiv 2311.16328 (2023) / PubMed 38076516; follow-up
  "Ligand-Based Compound Activity Prediction via Few-Shot Learning -
  MetaNet" J. Chem. Inf. Model. 10.1021/acs.jcim.4c00485.
- Review state: SNIPPET. Support-set conditioned compound activity
  prediction without target features (ligand-only); benchmarked against
  Tanimoto/kNN baselines (relevant: we must reproduce this comparison under
  our episode bank).
- Protocol fields: ligand-only by design (no protein conditioning); cold
  targets impossible by construction; confound control = ligand similarity,
  which is the intended signal. Direct template for the P1
  ligand-only + FS-CAP-style baseline pair.

## ZeroBind (protein-specific zero-shot)

- Identity: "ZeroBind: a protein-specific zero-shot predictor with subgraph
  matching for drug-target interaction prediction", Nat. Commun. 14 (2023),
  doi 10.1038/s41467-023-43597-1, PMC10687269.
- Review state: SNIPPET (PDF at PMC). Uses AlphaFold2 for proteins without
  structures; protein-ligand subgraph matching; k-shot protein-specific
  mode.
- Protocol fields (to verify): zero-shot protein generalization benchmark;
  external structure prediction = external-pretraining-like contribution
  that must be ablated separately under our rules.

## CNP lineage (support-conditioned regression)

- Identity: "Conditional Neural Processes for Molecules" arXiv 2210.09211
  (2022); Contrastive CNP arXiv 2203.03978.
- Review state: SNIPPET. CNP = encoder over support set -> global latent ->
  decoder; natural template for a support-encoder baseline and for
  uncertainty-aware support adaptation.
- Protocol fields: method papers; no DTA benchmark claims of their own.

## Leak-aware critiques (define our audit checks)

- HonestAffinity: "HonestAffinity: Leak-Aware Evaluation of Protein and
  Pocket Priors for Binding Affinity", arXiv 2606.03422 (2026).
  Review state: SNIPPET. Sequence-based models can win via pocket/protein
  priors correlated with the split (train/test leakage of target level);
  prescribes leak-aware controls. Mirrors our M1 counterfactual suite.
- "When Does Context Help? A Systematic Study of Target-Conditional
  Molecular Property Prediction", ICLR 2026
  (openreview.net/forum?id=eC00xzN6Jb). Review state: SNIPPET. Systematic
  target-context ablation; directly informs the P1 "how much is ligand
  similarity vs protein conditioning" decomposition.
- GEMS: "Enhancing Generalizable Binding Affinity Prediction by Removing
  Data Leakage", bioRxiv 10.1101/2024.12.09.627482v2. Review state:
  SNIPPET. Leakage-removal protocol for generalizable affinity
  benchmarks.

## Actions for this window

1. Full-text inspection of AdaMBind PDF + SI before any local AdaMBind
   implementation decision (protocol fields above must be filled).
2. Full-text inspection of HonestAffinity + ICLR 2026 context paper to
   fix the P1 decomposition design (ligand-only vs protein-conditioned
   attribution) before freezing the P-line preregistration.
3. Nothing adopted; P1 baseline bake-off stays our own measurements.
