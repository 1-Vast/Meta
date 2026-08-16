# OpenMut `OMUT-X5` decision

**Date:** 2026-07-28.
**Verdict:** `OMUT_X5_SUPPLEMENT_CONSTRUCT_REGISTRY_INADEQUATE_STOP`.
**Bound result SHA-256:**
`8cfaf5dcbaed803f93abb3461f3982bc0d4cb1cb7929f79e5bbdeb9021770806`.

## Decision

All execution, member-policy, PDF-review, and outcome-firewall gates passed.
The nine official Europe PMC supplementary endpoints returned nine valid ZIP
archives containing:

- nine parsed PDF supplements;
- three inventoried but unread CSV tables;
- 168 GIF and 168 JPEG image assets inventoried but not read.

Twelve safe construct fragments survived the corrected firewall. One
candidate/member relation passed both the text rule and Poppler visual
review: `P08581 D1228V`, where WT and mutant are assigned the same
`His6-cMet(1038-1348)` protein construct. The strict topology therefore rises
from four to five components but remains at three reported family categories,
well below the `25 / 6` gate.

The EGFR `L858R` PDF hit was visually rejected because it described a
multi-substitution TMLR crystallography construct, not the required single
L858R WT-pair construct.

## Interpretation

The licensed main-text and supplement route is now source-exhausted for these
nine EPMC records under exact construct rules. Lowering the four-ligand
threshold does not rescue the registry: strict description evidence yields
only six components at `k>=2`.

No affinity result table, CSV member, Davis value, or sealed-test value was
read. Real-outcome training remains blocked.

## Next allowed action

`OMUT-X6` may inspect the already-bound ChEMBL assay descriptions for exact
supplier catalog numbers, product URLs, or other source-native reagent
locators. X6 is label-free and may only decide whether an official
reagent-registry lookup is sufficiently populated to justify a separate
external verification stage.
