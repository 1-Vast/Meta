# SAFSA-G0 preregistration: family-selectivity identifiability gate

Date: 2026-07-27  
Role: cheapest identifying gate before structure recovery or neural pretraining.  
Candidate provenance: this direction was supplied in the user's attached research opinion and does
not consume an agent-proposed candidate slot.

## Question and scope

Test whether the systematic, inactive-retaining Novartis Secondary Pharmacology Database (SPD)
contains a leakage-resistant **protein-family selectivity signal** that can be extracted from ligand
response patterns and transferred to a left-out gene without using that gene's labels to construct
the score.

This is not a dual-cold affinity benchmark and cannot overturn the RAMCI-S0 power verdict. SPD is
used only as a pretraining/anchor substrate. The historical FORT confirmation partition remains
permanently quarantined and the sealed test remains unconsumed.

## Disclosed feasibility inspection before freezing this file

Only schema and aggregate-count checks were performed:

- 1,948 drugs, 101 mapped genes and six broad target classes;
- 84 genes have at least five active and twenty tested-negative gene-level cells before the final
  controls;
- 91.3% of the 144 unmerged-assay rows are censored tested negatives;
- the source-count fields are sparse, but most rows have no DrugCentral or ChEMBL quantitative or
  single-concentration evidence;
- no family-selectivity score, AUPRC, bootstrap interval or gate contrast was calculated.

These counts may determine feasibility exclusions below but may not change the frozen effect-size
thresholds.

## Frozen source and target-class map

Input files are the CC-BY-4.0 SPD release recorded and hash-verified in
`dataset/public/spd_2023/manifest.json`.

Use only the 144 assay groups whose published `assay_group` is a single identifier. The 24
comma-joined assay groups are excluded from the primary gate because they are release-level merged
records, not independent evidence. Rows marked `Multiple_EGs == Yes` are excluded. Map assay group
to gene using `assay_group_vs_gene_map.txt`, then collapse to one drug--gene cell only after these
filters: active if any retained assay row has prefix `=`, tested-negative if every retained row has
prefix `>`.

The six broad classes follow the target-class taxonomy used by the SPD paper and
IUPHAR/ChEMBL:

- kinase: `EGFR, ERBB2, KDR`;
- ion channel: `CACNA1B, CACNA1C, CHRNA4, GABRA1, GRIN1, HTR3A, KCNH2, SCN5A`;
- transporter: `ABCB11, SLC18A2, SLC29A1, SLC6A2, SLC6A3, SLC6A4`;
- nuclear receptor: `AR, ESR1, ESR2, NR1H2, NR1H3, NR1H4, NR1I2, NR1I3, NR3C1, PGR, PPARA,
  PPARG, RARA, RARB, RARG`;
- enzyme: `ACHE, F2, MAOA, PDE3A, PDE4A, PDE4D, PDE6A, PTGS1, PTGS2`;
- GPCR: every other mapped SPD gene.

The runner must fail closed if the mapped gene set is not exactly partitioned once.

## Provenance and selection controls

The primary matrix is **public-evidence-clean**: exclude an entire drug--gene cell if any
contributing row has positive DrugCentral AC50 count, ChEMBL AC50 count, or ChEMBL
single-concentration count. Subscription evidence and rows without a public source count remain,
because the intended anchor is the systematic Novartis panel and its retained negatives.

A mandatory sensitivity matrix additionally removes every drug--gene cell marked as a DrugCentral
mechanism-of-action interaction. It tests whether the signal is merely rediscovery of known primary
targets.

Missing cells are never converted to negatives. For every query gene and drug, a score is evaluable
only when the query cell is tested, at least one other same-family gene is tested, at least five
outside-family genes are tested, and at least two outside families contribute.

## Frozen zero-query-label scores

For each query gene, remove all its cells before constructing ligand scores. For a set of genes
with `p` observed actives among `n` tested cells, use the Jeffreys-smoothed rate

`r = (p + 0.5) / (n + 1)`.

Construct:

1. `own_family`: rate over other genes in the query gene's class;
2. `family_cold`: pooled rate over all genes outside the query class;
3. `wrong_family`: equal-weight mean of separately smoothed rates from the contributing outside
   classes;
4. `global_promiscuity`: rate over all non-query genes;
5. `selection_coverage`: fraction of other same-family genes tested for that ligand.

No query-gene label is used in any score. A gene is evaluable only if the common score subset has at
least five actives and twenty tested negatives.

Primary per-gene metric is average precision (AUPRC), appropriate for the inactive-dominated panel.
AUROC is descriptive only. Aggregate by mean within gene, then equal-weight mean across target
families so the 60-GPCR class cannot dominate.

## Frozen uncertainty

Seed is 1729. Use 10,000 hierarchical bootstrap replicates: sample target families with replacement,
then sample evaluable genes with replacement within each sampled family. Report percentile 95%
intervals. No assay-row or drug-row bootstrap may substitute for this independent-family unit.

## Frozen gates

All must pass:

1. primary public-evidence-clean analysis has at least 24 evaluable genes across at least four
   target families;
2. family-macro `AUPRC(own_family) - AUPRC(family_cold) >= 0.03` and its 95% lower bound is above
   zero;
3. family-macro `AUPRC(own_family) - AUPRC(wrong_family) >= 0.03` and its 95% lower bound is above
   zero;
4. family-macro `AUPRC(own_family) - AUPRC(selection_coverage) >= 0.03` and its 95% lower bound is
   above zero;
5. after removing DrugCentral MOA cells, at least 20 genes across at least four families remain and
   `AUPRC(own_family) - AUPRC(family_cold) >= 0.02` with lower bound above zero;
6. tested negatives are at least 80% of primary cells, family-macro weighting is used, current-run
   FORT confirmation labels are unread, and the sealed test is unconsumed.

`own_family - global_promiscuity` is reported as a non-gating sentinel: a negative value means the
anchor is not adding selectivity beyond broad ligand activity even if other gates happen to pass.

Pass verdict:
`SAFSA_G0_PASS_AUTHORIZE_DATASET_BUILD_AND_PRETRAINING_PREREGISTRATION`.

Fail verdict:
`SAFSA_G0_FAMILY_SELECTIVITY_NOT_IDENTIFIABLE_STOP`.

On failure, do not recover structures, download protein models, tune smoothing, merge assay groups,
lower the family-level uncertainty unit, or train SAFSA. First issue an overall failure report as
requested. On pass, separately freeze the dataset firewall and neural objective before pretraining.

