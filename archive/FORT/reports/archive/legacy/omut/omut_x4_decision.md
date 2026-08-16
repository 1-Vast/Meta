# OpenMut `OMUT-X4` decision

**Date:** 2026-07-28.
**Verdict:** `OMUT_X4_CONSTRUCT_REGISTRY_INADEQUATE_STOP`.
**Bound result SHA-256:**
`5decf908ea6e708045792f1b4792c1f614c2ed62d6d56bcb87ba733e0b331ac6`.

## Decision

All execution and firewall gates passed, but the source topology gate failed.
The corrected transport dispositions are:

- nine complete Europe PMC article XML documents;
- two `invalid_document` Elsevier responses that contain only work metadata
  and no article blocks.

Only two outcome-free construct-related fragments were projected from the
nine full texts. Neither met a registered candidate-level common-construct
form. Zero near-exact candidate/document pairs were accepted. The remaining
registry is exactly the four X1 description-exact components across three
reported family categories, below the 25-component / six-family floor.

The initial 11-complete transport count was corrected under
`omut_x4_transport_correction.md`; the negative topology was unchanged.

## Interpretation

The X3C metadata pass did not survive source-body verification. Open article
availability alone is not sufficient: the WT-to-mutant construct relation is
missing from the main-text records at the specificity needed to bind ChEMBL
assays.

This stops X4 and keeps affinity training blocked. It does not establish that
the nine EPMC supplementary archives lack construct evidence, because X4
froze only `fullTextXML`.

## Next allowed action

`OMUT-X5` may query the official Europe PMC `supplementaryFiles` endpoint for
the same nine PMCIDs, inventory the returned archives, and project
candidate-level construct evidence from text-bearing supplements under the
same outcome firewall. No search-derived or unrelated source may enter.
