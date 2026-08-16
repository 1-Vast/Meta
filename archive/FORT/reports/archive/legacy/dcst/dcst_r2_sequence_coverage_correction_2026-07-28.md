# DCST-R2 exact-sequence coverage correction

Date: 2026-07-28  
Status: frozen before fetching missing UniProt sequences  
Corrects only: the local-cache coverage statement in the R2 preregistration

## Trigger

The pre-retrieval audit established 2,344 firewalled rows with one
registry-consistent pocket UniProt (2,124 train, 220 development), but did not
check whether every accession was present in the historical local
`uniprot_sequences.json`. The first RCSB registry build exposed that omission:
the local cache covered 1,664 rows and 733 of the 863 exact accessions.

This is a coverage/accounting error, not a model result. RCSB retrieval itself
succeeded for all 1,366 entities requested for the locally covered subset.

## Correction

Before finalizing the R2 registry, fetch the 130 missing, pre-specified pocket
accessions from the official UniProt REST API using only accession and
sequence fields. Update the local sequence cache, then rebuild the same RCSB
entity list. No affinity, assay, split outcome, model score, or certificate is
used to choose accessions.

Report requested, recovered, and unresolved accessions separately. All R2
gates, steps, seed, feature model, mapping rule, and stop conditions remain
unchanged.

