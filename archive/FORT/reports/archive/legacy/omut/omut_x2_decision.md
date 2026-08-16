# OpenMut `OMUT-X2` decision

**Verdict:** `OMUT_X2_FULLTEXT_RECOVERY_INSUFFICIENT_STOP`

**Date:** 2026-07-28.
**Preregistration:** `reports/active/omut_x2_preregistration.md`,
SHA-256 `fddbe2e4598653d867eb3526ba8eb644b0bee658e7c7c099f2ce888522c4a0d8`.
**Runner:** `research/omut_x2.py`.
**Result:** `reports/active/omut_x2.json`.
**Predecessor:** `OMUT-X1`.

## 1. Execution

All execution, firewall, candidate-disposition, near-pair, document, Europe
PMC, discovery-exclusion, and family gates pass. The only failure is
`X2_SOURCE_TOPOLOGY_ADEQUATE`.

X2 found 87 primary documents carrying ligand pairs that satisfy every frozen
non-construct condition. Europe PMC matched 83 identifiers, did not match two,
and dispositioned two as patents. Nine documents have an open-access EPMC
full-text record. No full-text body, activity value, relation, censor, unit,
observed difference, empirical variance, or empirical MDE was read.

## 2. Source-recoverable upper bound

Combining X1 exact ligands with near-exact ligands from the nine EPMC
open-full-text documents gives:

```text
source-recoverable components       = 16
broad families                      = 6
largest accession share             = 0.375
required components                 = 25
```

The family and concentration conditions pass, but the component floor fails.
Full-text availability is still only an optimistic bound: it does not prove
that a paper states the construct.

If all 87 near-pair documents were legally available, the post-run diagnostic
upper bound would be 110 non-BRAF components across 17 broad families, with a
largest accession share of 0.145. The route is therefore constrained by
licensed source discovery rather than candidate topology.

## 3. Decision

Europe PMC open-full-text records alone do not authorize a full extraction
campaign. `OMUT-I0` and affinity training remain blocked.

The next admissible label-free action is an exact DOI lookup in an independent
open-access index. Only accepted or published versions with an explicit open
license and a reproducible legal location may augment the nine EPMC records.
Preprints, snippets, social sharing sites, unlicensed PDFs, and paywall
landing pages do not count.
