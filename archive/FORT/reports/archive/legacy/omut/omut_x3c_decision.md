# OpenMut `OMUT-X3C` decision

**Date:** 2026-07-28.
**Verdict:** `OMUT_X3C_LICENSED_LINK_RECOVERY_FEASIBLE`.
**Scope:** label-free Crossref license/link metadata only.

## Decision

All seven frozen gates passed. Exact-DOI Crossref lookup dispositioned all 87
near-pair documents: 85 matched records and two documents had no DOI. Two
previously unavailable documents carried both an active admissible Creative
Commons license and a version-matched HTTPS full-text XML link. Union with
the nine EPMC records produced 11 optimistic accessible documents.

The resulting non-BRAF upper-bound topology is:

- 32 `k>=4` components;
- seven reported broad-family categories;
- ten accessions;
- largest accession share `0.50`.

This clears the frozen `25 / 6 / <=0.50` gate without changing a threshold.
No affinity outcome, abstract, linked full text, Davis value, or sealed-test
value was read by X3C.

## Concentration warning

The pass is deliberately narrow. `P61073` contributes 16 of the 32
components, all through one Crossref document
(`CHEMBL1156482`, DOI `10.1074/jbc.M704739200`) and exactly meets the maximum
share. The other new Crossref document is `CHEMBL4387756`, DOI
`10.1016/j.bmc.2019.05.040`. Both projected links are Elsevier text-mining XML
endpoints. Crossref's work-level license/link combination does not establish
unauthenticated transport, a shared WT construct, or usable labels.

Therefore this verdict authorizes only a preregistered transport and
construct-evidence audit. It does not authorize affinity training.

## Next allowed action

`OMUT-X4` may fetch only the 11 frozen source bodies, retain only
outcome-free construct evidence, and recompute topology using candidate-level
construct dispositions. `OMUT-I0` remains blocked unless X4 preserves the
source topology under exact evidence.
