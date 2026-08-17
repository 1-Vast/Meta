# Stage D literature evidence (design evidence only; nothing here binds the design)

Review status per source: abstract-level (PubMed/PNAS eutils), full text not
fetched unless noted. These sources informed the Stage E design; the
preregistered hypotheses stand or fall on the D0 measurements and the
single-seed screen, not on these papers.

1. Luo et al., BatchDTA: implicit batch alignment enhances deep learning-based
   drug-target affinity estimation. Brief Bioinform 23(4):bbac260, 2022.
   doi:10.1093/bib/bbac260. PubMed 35794723. Inspected: abstract.
   Claim: measure metrics / assay information create batch effects that bias
   DTA labels; implicit alignment through learned compound orders across
   batches improves four DNN backbones (BindingDB, Davis, KIBA).
   Use: supports treating the level term as assay/batch-affected and training
   rank-stable (order) objectives; our D0_LEVEL_ANATOMY documents 70%
   in-fold document variance directly. Contradiction: none recorded.
   Uncertainty: abstract-level only; the published algorithm details were not
   inspected.

2. Gorantla et al., Learning Binding Affinities via Fine-Tuning of Protein and
   Ligand Language Models. J Chem Inf Model 65(22):12279-12291, 2025.
   doi:10.1021/acs.jcim.5c02063. PubMed 41171175. Inspected: abstract.
   Claim: BALM (protein/ligand LMs) generalizes to unseen drugs, scaffolds and
   targets on curated BindingDB; few-shot support beats ML and Vina on USP7/Mpro.
   Use: independent evidence that strict target-based splits are feasible and
   that LM features transfer; the local ESM-2 650M lane (D0) follows this.
   Uncertainty: their splits and metric definitions were not inspected in full.

3. Graber et al., Resolving data bias improves generalization in binding
   affinity prediction. Nat Mach Intell 7:1713-1725, 2025.
   doi:10.1038/s42256-025-01124-5. PubMed 41143208. Inspected: abstract.
   Claim: PDBbind/CASF leakage inflated prior models; CleanSplit retraining
   drops them substantially while a sparse-graph + LM-transfer model holds up.
   Use: independent support for leak-free protocols (this project's Stage B
   measured the same pattern: meta_val checkpoint selection worth ~0.62 pK^2).
   Uncertainty: abstract-level only.

4. Singh et al., Contrastive learning in protein language space predicts
   interactions between drugs and protein targets (ConPLex). PNAS
   120(24):e2220778120, 2023. doi:10.1073/pnas.2220778120. PubMed 37289807.
   Inspected: abstract.
   Claim: protein-anchored contrastive coembedding improves zero-shot DTI
   prediction and decoy specificity.
   Use: candidate future framework direction (contrastive coembedding); not
   used in Stage E because it changes two modules at once and the D0 evidence
   points at level/panel structure first.

5. Graph neural processes for molecules: an evaluation on docking scores and
   strategies to improve generalization. J Cheminform 16:94, 2024.
   doi:10.1186/s13321-024-00904-2. Inspected: discovered via search; abstract
   not inspected (redirect wall).
   Use: kept as a conditional-neural-process comparison lane; no claims built
   on it.

Design implications recorded (not evidence): (a) level/shape separation must
be explicit in the loss; (b) order-based objectives transfer across assays
better than absolute calibration; (c) every claim about generalization
requires a leakage-audited split.
