# RB-DR-QMAPD — Stage O0 feasibility / power audit (2026-07-26)

Model-free, cheapest gate first. Metz development panel only; no label read. Runner
`research/rb_dr_qmapd_o0.py`; machine-readable `reports/active/rb_dr_qmapd_power_audit.json`
(registry sha256 `94da6bb5…173e`).

## Construction

One label-blind permutation per development target; greedy scaffold-distinct support of size `K`;
fixed dual-cold query = ligands whose scaffold is disjoint from the support and whose max Tanimoto to
the support is `< 0.95`. Nested `S_4 ⊂ S_8 ⊂ S_16 ⊂ S_32`. The evaluated population for a pair `(k,K)`
is the set eligible at `K`. Dual disjointness at the mechanism gate: held homology component
(leave-component-out CV fold) + query scaffold disjoint from support.

## Result

| pair (k,K) | eligible targets | independent homology components | resolvable (≥30) |
|---|---:|---:|---|
| (4, 16) | 80 | 76 | yes |
| (4, 32) | 42 | 40 | yes |
| (8, 32) | 42 | 40 | yes |

Frozen reference MDE80 (panel_power_k4.json): paired k=4 component-macro Spearman MDE80 `0.0367`
(retrain-noise SD 0.1073/component). The empirical paired MDE80 for `Delta_info`/`Delta_total` is
computed inside O1 from a support-resampling null (same teacher, two independent support draws), not
assumed here.

Provenance limitation (honest): the panel is a **single document** (CHEMBL1201862) with **assay 1:1
with target**, so document overlap is trivially 1 and support/query necessarily share the target's
assay. This is a within-kinome, single-provenance substrate; it cannot test cross-document or
cross-assay generalisation. Recorded, not hidden.

## Verdict and selection

```
RB_DR_QMAPD_O0_FEASIBLE
```

Primary pair (frozen rule = most independent components → best power; uses only label-free counts):
**(k=4, K=16)**, 76 components. Secondary large-jump robustness arm (reported, non-gating):
**(k=4, K=32)**, 40 components. Proceeding to teacher preregistration and Stage O1.
