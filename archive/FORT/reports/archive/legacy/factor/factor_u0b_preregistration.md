# FACTOR-U U0-B preregistration: mixed strict unlabeled-corpus audit

Date: 2026-07-26  
Route origin: user-supplied post-F0-C unlabeled-corpus expansion; no agent candidate slot consumed.  
Role: data-only audit. U0-B trains no encoder or affinity model.

## Why U0-B is allowed

PLINDER-only U0 retained 29,382 molecules and 17,003 scaffolds but failed because one scaffold
contributed 1.266%, above 0.5%. U0-B does not downsample that result or relax its threshold. It
changes the data source by forming a deduplicated union with the local ChEMBL 37 training
structures, testing whether an independently curated medicinal-chemistry collection removes the
structural concentration.

ChEMBL states that its data are provided under CC BY-SA 3.0. The local release is ChEMBL 37
(May 2026). Before redistribution, the database terms and share-alike obligations must be retained.

## Frozen sources and projections

1. PLINDER 2024-06/v2: exactly the seven structure/quality fields and molecular filter frozen in U0.
2. ChEMBL 37 dual-cold registry: predicate-pushdown `dual_cold_split == "train"` and project only
   `conn` and `scaffold`. Do not project affinity, endpoint, target, assay, document, accession or
   confirmation rows.

Canonicalize, deduplicate the union and retain 6--80 heavy atoms with elements in
`{B,C,N,O,F,Si,P,S,Cl,Se,Br,I}`.

## Global strict firewall

Delete every union molecule whose parent connectivity or nonempty Murcko scaffold occurs anywhere
in KIRHub2026, Reinecke2024 or Papyrus-Christmann2016. This is applied after source union and before
any statistic. The existing ChEMBL confirmation partition remains permanently quarantined; its
rows are not materialized.

## Frozen gates

All must pass:

1. zero connectivity and scaffold overlap with each evaluation source;
2. at least 50,000 retained unique molecules;
3. at least 25,000 retained nonempty scaffolds;
4. no single scaffold contributes more than 0.5% of retained unique molecules;
5. at least 10,000 retained molecules exclusive to each input source before union;
6. all evaluation element classes and pharmacophore roles supported;
7. retained heavy-atom q01 <= evaluation q05 and retained q99 >= evaluation q95;
8. both licenses and raw/registry SHA-256 values recorded;
9. current-run activity/affinity/protein columns read = false;
10. current-run confirmation labels read = false and sealed test consumed = false.

Pass: `FACTOR_U0B_PASS_AUTHORIZE_U1_PREREGISTRATION`.  
Fail: `FACTOR_U0B_CORPUS_INELIGIBLE_STOP`.

If U0-B passes, U1 still requires a separate fixed model/training preregistration. No F0-C1
checkpoint, external label or after-result corpus balancing may be reused.
