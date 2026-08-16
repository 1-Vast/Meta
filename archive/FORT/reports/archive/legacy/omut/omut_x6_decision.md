# OpenMut `OMUT-X6` decision

**Date:** 2026-07-28.
**Verdict:** `OMUT_X6_REAGENT_VERIFICATION_INSUFFICIENT_STOP`.
**Bound result SHA-256:**
`0de65b7efab34df6b659945853b865f9213385519f75327b65ff5def52b1c626`.

## Decision

All execution and firewall gates passed. X6 scanned the 145 distinct WT
assays used by frozen X2 near pairs and found:

- zero labeled catalog/product tokens;
- zero HTTPS product URLs;
- zero actionable supplier locators;
- zero candidate components eligible for an external supplier verification.

No supplier page or activity outcome was read. The reagent-registry route is
closed for these ChEMBL descriptions.

## Next allowed action

BindingDB remains the only local source with a materially larger
sequence-exact upper bound: 37 `Ki/Kd, k>=4` components. X0 excluded all of
them because the archive has no source-native assay identifier. A new stage
may test one conservative alternative only: an exact composite context
signature requiring same document, provenance fields, nonmissing pH and
temperature, and exact equality of those conditions between WT and mutant
rows. Missing conditions cannot pass the primary signature.
