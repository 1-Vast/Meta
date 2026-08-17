# Stage L report — support-gated assay-aware level head: REJECTED (G3)

Development evidence, single seed, meta_val read once after freezing;
meta_test sealed. Authorities: L_vs_T2.contrast.json, L_meta_val.rows.summary.json,
PREREGISTRATION.md, RESULT.json.

## Verdict

**Rejected by G3; nothing promoted.** The support gate did its preregistered
job for MSE — k=1/2/3/5 stay statistically indistinguishable from T2 (all
intervals cross zero, means slightly better) — but ordering degrades with
RESOLVED intervals at k=2/3/5 (Spearman -0.0865 / -0.0751 / -0.0596; CI
-0.0399 / -0.0345 / -0.0256) and k=0 CI -0.0300 [-0.0681, -0.0022].
Stop rule S2 fires.

The k=0 result is the best calibration in the project's record: MSE 2.0997
(-0.4964 vs T2, unresolved), level^2 1.2151 (-0.5163) — but the k=0
ranking means are lower (Spearman 0.019 vs 0.079), and the k>=1 ordering
regression shows the gate alone is not enough: training the level head on
k=0 episodes reshapes the SHARED trunk, which then orders worse at k>=1
even with the head switched off.

## What this closes

Three compositions of a learned zero-shot level head with this trunk have
now failed the ranking gate: ungated with residual supervision (Stage E),
ungated with assay covariates (Stage J), gated by support size (Stage L).
The coupling is through the shared encoder, not the head's output: the
zero-shot level objective and the within-target ordering objective are in
conflict on the same representation, and no routing rule tested separates
them without degrading one or the other with resolved intervals.

A fully separate (frozen-feature) level calibrator would avoid the shared
trunk conflict, but its measured ceiling is the D0c journal probe (1.62
level MSE -> k=0 MSE ~2.5), strictly worse than L's 2.10, and it would be
a two-stage composition that the governing contract discourages. The
protocol-level conclusion in report/BOUNDARY_20260817_NIGHT.md therefore
stands as the final state of this research programme.
