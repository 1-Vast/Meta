# FIRE-DTA D0 Structural Registry

## Verdict

`FIRE_DTA_D0_PASS`

The label-blind PLINDER registry contains 9,264 training systems and
358 validation systems after structural-quality filters and all
registered target, protein, pocket, compound, scaffold, near-neighbor and edge firewalls. No affinity
column was read and the PLINDER test split was not selected.

## Gates

| Gate | Observed | Required | Pass |
| --- | ---: | ---: | --- |
| Training systems | 9,264 | 5,000 | True |
| Protein clusters | 1,375 | 500 | True |
| Pocket clusters | 2,696 | 500 | True |
| Official apo-or-pred links | 2,902 | 1,000 | True |
| Experimental apo identity certified | 1,828 | 1,000 | True |

The PLINDER metadata flag alone does not distinguish apo from predicted partners. The separately frozen
official apo link table and extracted structure inventory resolve that ambiguity and certify the reported
experimental apo count. Predicted structures remain separately typed and cannot satisfy the apo gate.

## Firewall

Protected axes include all ChEMBL target accessions, all canonical evaluation target accessions, their
observed PLINDER 70% protein-homology and pocket clusters, and all canonical evaluation compounds,
scaffolds and Morgan-neighborhood matches at Tanimoto >= 0.95. Edge overlap is zero by conservative
target-axis and compound-axis quarantine; no label-bearing edge registry was parsed.
