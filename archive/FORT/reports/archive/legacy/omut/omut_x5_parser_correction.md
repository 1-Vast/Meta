# OpenMut `OMUT-X5` parser/firewall correction

**Date:** 2026-07-28.
**Scope:** implementation correction under frozen X5 evidence form 3 and the
frozen outcome firewall.

The first X5 execution used two patterns that were narrower than the written
preregistration:

1. an explicit construct span was recognized only after the words
   `residue(s)` or `amino acid(s)`, so
   `His6-cMet(1038-1348)` was missed;
2. `affinity` was blocked but its plural `affinities` was not.

Neither change alters a source, threshold, or semantic evidence form.
Projection schema 2 adds a protein-name parenthetical span pattern and blocks
both singular and plural outcome terms. The old projections cannot be reused.

The c-MET candidate page was rendered with Poppler and visually inspected
under the PDF skill. The relevant sentence assigns both WT and D1228V to the
same His6-cMet(1038-1348) construct. The visual disposition is bound in
`omut_x5_pdf_review.json`.
