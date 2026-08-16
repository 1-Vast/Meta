# OpenMut `OMUT-X1` decision

**Verdict:** `OMUT_X1_DESCRIPTION_REGISTRY_INADEQUATE_STOP`

**Date:** 2026-07-28.
**Preregistration:** `reports/active/omut_x1_preregistration.md`,
SHA-256 `567e538829a677bac18ae498528488f0107f6dbd112d74cbf12373392faa69eb`.
**Runner:** `research/omut_x1.py`.
**Result:** `reports/active/omut_x1.json`.
**Predecessor:** `OMUT-X0`, verdict
`OMUT_X0_EVIDENCE_INADEQUATE_STOP`.

## 1. Execution and firewall

All execution, binding, cache, projector, candidate-disposition, exact-pair,
discovery-exclusion, and family-mapping gates pass. The sole failed gate is
`X1_PRIMARY_TOPOLOGY_ADEQUATE`.

The run reproduced the complete 119,801-row ChEMBL 37 variant projection,
completed all 143 deferred WT queries, resolved all 2,622 assay locators and
all 1,381 document locators, and dispositioned all 279 frozen k-capable
groups. It materialized zero numeric affinity, relation, censor, unit,
observed difference, empirical variance, reliability, or empirical MDE.
DAVIS-Complete and sealed test data were not read.

## 2. Exact components recovered

Four non-discovery components clear the exact per-ligand `k=4` rule:

| accession | mutation | endpoint | construct span | paired ligands | documents |
| --- | --- | --- | --- | ---: | ---: |
| `A3EZI9` | `R155K` | `Ki` | 1-631 | 6 | 1 |
| `O60674` | `V617F` | `Kd` | 536-812 | 4 | 1 |
| `P00533` | `L858R` | `Kd` | 669-1011 | 13 | 5 |
| `Q5S007` | `G2019S` | `Ki` | 970-2527 | 13 | 1 |

Every retained ligand has a source-native activity ID on each side, a
source-native assay ID on each side, one primary document identity, one
endpoint, the same explicit construct span, an exact normalized description
signature, and a mutation sequence consistent with the versioned UniProt
reference.

The pre-preregistered discovery accession `P15056` contributes zero primary
components. Its BRAF V600E `Kd` component recovers one exact ligand at span
429-741 and remains discovery-only below `k=4`.

## 3. Topology and decision

The primary topology is:

```text
components                 = 4
broad families             = 3
largest accession share    = 0.25
required components        = 25
required broad families    = 6
```

Twenty-five additional candidates recover one to three exact ligands but do
not clear `k=4`. The exact-description strategy is therefore a real recovery
mechanism, not a sufficient substrate.

`OMUT-I0`, representation fitting, and real affinity training remain blocked.
No threshold, fuzzy text similarity, canonical-WT inference, assay adjacency,
model capacity, or seed expansion may rescue X1.

## 4. Next admissible action

The next label-free action is an X2 source-accessibility audit restricted to
candidate ligand pairs that already match on molecule, target, endpoint,
primary document, mutation sequence, assay type, and exact normalized assay
context, and fail only because the descriptions do not state an explicit
construct span.

X2 may query DOI/PMID, open-full-text and supplement availability. It may not
read activity outcomes or use an LLM/fuzzy matcher to declare a construct.
Only a powered source-accessibility topology can authorize deterministic,
quote-bound full-text construct extraction.
