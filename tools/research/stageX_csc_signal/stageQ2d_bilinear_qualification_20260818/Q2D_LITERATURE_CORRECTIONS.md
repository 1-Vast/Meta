# Q2d literature corrections (per independent review, 2026-08-18)

These corrections supersede the earlier literature reading in the Q2c/Q3b
answer and are the operative evidence map for Stage Q2d.

| Claim (old) | Correction | Evidence |
|---|---|---|
| Merget et al. 2017 = 234-kinase x 380-inhibitor Kronecker PCM proving unseen-kinase transfer | WRONG. Merget et al. 2017 (J. Med. Chem.) trains ligand-side random-forest models per kinase; it is not a Kronecker PCM and does not demonstrate unseen-kinase transfer. Withdrawn as PCM evidence. | https://pubmed.ncbi.nlm.nih.gov/27966949/ |
| The literature uniformly uses delta-delta targets | WITHDRAWN. Many PCM works predict absolute activity/classification. Delta-delta is chosen HERE because it structurally zeroes ligand and protein main effects, not because all successful papers use it. | - |
| MdrDB = clean cold-protein benchmark (~100k measured WT-mutant dG records) | WRONG as stated. MdrDB mixes GDSC/DepMap-derived samples with structure-derived entries (PyMOL/AlphaFold); it is not 100k independent directly measured biochemical WT-mutant affinity pairs. It may serve as an EXTERNAL STRESS TEST only, after auditing label origin, structure replication, cell-line vs affinity mixing, protein/mutation/drug duplicates, and true component-isolated cold splits. | https://pmc.ncbi.nlm.nih.gov/articles/PMC10267113/ |
| Duong-Ly = same-study WT->mutant positive control | PARTIALLY WRONG. Duong-Ly WT values come from Anastassiadis (cross-study batch). B1 primary must be the Saifudeen same-study panel; Duong-Ly is B1-R cross-study replication with explicit study/batch covariates. | Q3 census; Duong-Ly supplement provenance |

Operative positive PCM evidence (first-hand, cited):
1. 317 WT/mutant kinases x 38 inhibitors x 12,046 Kd with explicit
   protein-ligand cross descriptors - the original kinome-wide PCM.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC2910025/
2. HIV variant PCM: mutation-site physicochemical descriptors x drug
   descriptors express drug-specific resistance direction.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC3002298/ ; external antivirogram
   validation https://pmc.ncbi.nlm.nih.gov/articles/PMC3578754/
3. Active-site sequence DTA on BindingDB: active-site sequence outperforms
   full sequence under unseen-kinase and stricter unseen-ligand+kinase
   splits. https://pmc.ncbi.nlm.nih.gov/articles/PMC9516689/

Consequence for Q2d: the 'local protein representation x ligand
representation' direction is literature-grounded; the unproven claim is
transfer across strict cold parent/pocket-component splits - which is
exactly what Q2d -> B1 -> C must test.
