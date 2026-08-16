# STATUS: `panel_gate_pd0.json` is VOID

`reports/active/panel_gate_pd0.json` was produced by a run that was launched before the audit
verdict `BLUEPRINT_REQUIRES_MAJOR_REVISION` arrived and that finished at 19:07 on 2026-07-25 before
the process could be stopped. The file is preserved **byte-identically** and is not edited, because
experimental artifacts are never modified after the fact.

It is void as a Gate PD0 result, for reasons that are structural rather than statistical. The
implementation predates the corrected contract in `reports/active/orrc_eb_blueprint_v2.md` and
differs from it in three ways that each change the estimand:

1. it used the unweighted projector `M_X` rather than the specified `M_X^W`;
2. its latent block used full-grid two-sided feature orthogonality, which section C.1 of the revised
   blueprint proves is **not** observed-space orthogonality on an incomplete panel — one missing cell
   out of nine already aligns the "orthogonal" latent block 29% with the interaction direction;
3. its candidate grid was filtered to effective rank `<= 8`, which reimposes a non-convex rank
   constraint through the selection rule.

Its statistics are therefore not cited, interpreted or carried forward anywhere. The only thing
retained from the run is a bug diagnosis: the fitted auxiliary-feature energy fell monotonically as
the latent penalty was relaxed, which is exactly the leakage defect (2) predicts.

The same applies to `reports/active/panel_pd0_preregistration.md` (superseded; see its amendment) and
to `research/panel_gate_pd0.py`, `research/panel_gate_pd1.py` and `research/orrc.py`, which carry
`SUPERSEDED` headers and must not be run until the corrected modules of blueprint section G exist.

`reports/active/panel_power_pd1.json` remains a valid arithmetic record of arm heterogeneity in the
Gate PB contrasts, but the gate it was written for (PD1 on the Gate PB development rows) is withdrawn
by blueprint section F, so its value is a reference number and not a live threshold.

No confirmation or sealed label was read at any point.
